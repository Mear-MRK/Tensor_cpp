#pragma once
#include <iostream>
#include <stdexcept>
#include <string>
#include <cstring>
#include <cstdio>
#include <initializer_list>
#include <complex>
#include "Payload.hpp"

#ifdef DEBUG
#define LOG(s) std::cout << s << std::endl
#else
#define LOG(s)
#endif



inline void comp_index(size_t l_index[], size_t r_index[], size_t index[],
	unsigned l_rank, unsigned r_rank, unsigned l_i, unsigned r_i)
{
	unsigned r = 0;
	for (unsigned i = 0; i < l_rank; i++)
		if (i != l_i) {
			index[r++] = l_index[i];
		}
	for (unsigned i = 0; i < r_rank; i++)
		if (i != r_i) {
			index[r++] = r_index[i];
		}
}

inline void decomp_index(size_t index[], size_t l_index[], size_t r_index[],
	unsigned l_rank, unsigned r_rank, unsigned l_i, unsigned r_i)
{
	unsigned i = 0;
	while (i < l_i)
	{
		l_index[i] = index[i];
		i++;
	}
	unsigned j = i++;
	while (i < l_rank)
	{
		l_index[i++] = index[j++];
	}
	i = 0;
	while (i < r_i)
	{
		r_index[i++] = index[j++];
	}
	i++;
	while (i < r_rank)
	{
		r_index[i++] = index[j++];
	}
}

template<typename F>
inline std::string to_str(std::complex<F> c) {
	return std::to_string(c.real()) + "+" + std::to_string(c.imag()) + "i";
}

template<typename F>
inline std::string to_str(F f) {
	return std::to_string(f);
}



template<typename F = float>
class Tensor {

private:
	char meta[64] = "";
	unsigned rank = 0;
	size_t size = 0;
	size_t off = 0;
	size_t* dim = nullptr;
	size_t* jump = nullptr;
	Payload<F> data;


	inline void reset() {
		rank = 0;
		size = 0;
		off = 0;
		delete[] dim;
		dim = nullptr;
		delete[] jump;
		jump = nullptr;
		data.release();
	}

	inline void init(const size_t* t_dim, unsigned t_rank)
	{
		if (!t_rank) {
			return;
		}
		for (unsigned r = 0; r < t_rank; r++)
			if (t_dim[r] == 0)
				return;

		rank = t_rank;
		dim = new size_t[rank];
		memcpy(dim, t_dim, t_rank * sizeof(size_t));
		jump = new size_t[rank];
		size = 1;
		jump[rank - 1] = 1;
		for (unsigned i = rank - 1; i > 0; i--) {
			size *= dim[i];
			jump[i - 1] = size;
		}
		size *= dim[0];
	}

	inline void init(const std::initializer_list<size_t>& dim_l) {
		init(dim_l.begin(), dim_l.size());
	}

public:

	Tensor()
	{
		LOG("tns cns, dfl " << this << " : " << meta);
	}

	Tensor(const std::string& meta_str)
	{
		set_meta(meta_str);
		LOG("tns cns, mta " << this << " : " << meta);
	}

	Tensor(size_t size) {
		if (size) {
			size_t dim[] = { size };
			init(dim, 1);
			data = Payload<F>(size);
		}
		LOG("tns cns, siz " << this << " : " << meta);
	}

	~Tensor()
	{
		data.release();
		delete[] jump;
		delete[] dim;
		LOG("tns dst,     " << this << " : " << meta);
	}

	// copy constructor
	Tensor(const Tensor<F>& rhs) {
		if (rhs.size) {
			rank = rhs.rank;
			size = rhs.size;
			off = rhs.off;
			dim = new size_t[rank];
			memcpy(dim, rhs.dim, rank * sizeof(size_t));
			jump = new size_t[rank];
			memcpy(jump, rhs.jump, rank * sizeof(size_t));
			data = rhs.data;
		}
		set_meta(rhs.get_meta() + "Cc;");
		LOG("tns cns, cop " << &rhs << " -> " << this << " : " << meta);
	}
	
	// copy assignment
	Tensor<F>& operator=(const Tensor<F>& rhs) {
		if (this == &rhs || !rhs.size)
			return *this;
		if (!size) {
			rank = rhs.rank;
			size = rhs.size;
			off = rhs.off;
			dim = new size_t[rank];
			memcpy(dim, rhs.dim, rank * sizeof(size_t));
			jump = new size_t[rank];
			memcpy(jump, rhs.jump, rank * sizeof(size_t));
			data = rhs.data;
		}
		if (size) {
			data.assign(rhs.data, rhs.off, size, off, rhs.size);
		}

		this->append_meta("]" + rhs.get_meta() + "C=,");

		LOG("tns asn, cop " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	// move constructor
	Tensor(Tensor<F>&& rhs) noexcept {
		rank = rhs.rank;
		size = rhs.size;
		off = rhs.off;
		dim = rhs.dim;
		jump = rhs.jump;
		data = std::move(rhs.data);

		rhs.rank = 0;
		rhs.size = 0;
		rhs.off = 0;
		rhs.dim = nullptr;
		rhs.jump = nullptr;

		set_meta(rhs.get_meta() + "Mc;");

		LOG("tns cns, mov " << &rhs << " -> " << this << " : " << meta);

	}
	
	// move assignement
	Tensor<F>& operator=(Tensor<F>&& rhs) noexcept {
		if (this == &rhs || !rhs.size)
			return *this;

		if (!size) {
			rank = rhs.rank;
			size = rhs.size;
			off = rhs.off;
			dim = rhs.dim;
			jump = rhs.jump;
			data = std::move(rhs.data);
		}
		if (size) {
			data.assign(rhs.data, rhs.off, size, off, rhs.size);
			rhs.data.release();
		}
		rhs.rank = 0;
		rhs.size = 0;
		rhs.off = 0;
		rhs.dim = nullptr;
		rhs.jump = nullptr;

		append_meta("]" + rhs.get_meta() + "M=,");

		LOG("tns asn, mov " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	Tensor(const Payload<F>& data) : data(data) {
		if (data.get_size()) {
			size_t dim[] = { data.get_size() };
			init(dim, 1);
		}
		LOG("tns cns, pyl " << this << " : " << meta);
	}

	Tensor(const size_t dim[], unsigned rank) {
		init(dim, rank);
		if (size) {
			data = Payload<F>(size, 0);
		}
		set_meta("0;");
		LOG("tns cns, dr0 " << this << " : " << meta);
	}

	Tensor(const size_t dim[], unsigned rank, const Payload<F>& t_data) {
		init(dim, rank);
		if (size <= t_data.get_size()) {
			data = t_data;
		}
		else {
			reset();
		}
		set_meta("drP;");
		LOG("tns cns, drP " << this << " : " << meta);
	}

	Tensor(const size_t dim[], unsigned rank, F value) {
		init(dim, rank);
		if (size) {
			data = Payload<F>(size, value);
		}

		set_meta(std::to_string(value) + ";");
		LOG("tns cns, drv " << this << " : " << meta);

	}

	Tensor(const size_t dim[], unsigned rank, const F* array) {
		init(dim, rank);
		if (size) {
			data = Payload<F>(size);
			std::copy(array, array + size, data.get_arr());
		}
		set_meta("array;");

		LOG("tns cns, dra " << this << " : " << meta);

	}

	Tensor(const size_t dim[], unsigned rank, F(*rnd)(void)) {
		init(dim, rank);
		if (size) {
			data = Payload<F>(size);
			for (size_t i = 0; i < size; i++)
				data[i] = rnd();
		}
		set_meta("rnd;");

		LOG("tns cns, drr " << this << " : " << meta);

	}

	Tensor(const size_t dim[], unsigned rank, F start, F end) {
		init(dim, rank);
		if (size) {
			data = Payload<F>(size, start, end);
		}
		set_meta("lin;");
		LOG("tns cns, drl " << this << " : " << meta);
	}

	Tensor(const std::initializer_list<size_t>& dim_l) : Tensor(dim_l.begin(), dim_l.size()) {}

	Tensor(const std::initializer_list<size_t>& dim_l, const Payload<F>& data) : Tensor(dim_l.begin(), dim_l.size(), data) {
	}

	Tensor(const std::initializer_list<size_t>& dim_l, F(*rnd)(void)) : Tensor(dim_l.begin(), dim_l.size(), rnd) {
	}

	Tensor(const std::initializer_list<size_t>& dim_l, const F arr[]) : Tensor(dim_l.begin(), dim_l.size(), arr) {
	}

	Tensor(const std::initializer_list<size_t>& dim_l, F start, F end) : Tensor(dim_l.begin(), dim_l.size(), start, end)
	{}

	Tensor<F> copy() {
		Tensor<F> out;
		if (size) {
			out.rank = rank;
			out.size = size;
			out.off = 0;
			out.dim = new size_t[rank];
			memcpy(out.dim, dim, rank * sizeof(size_t));
			out.jump = new size_t[rank];
			memcpy(out.jump, jump, rank * sizeof(size_t));
			out.data = data.copy(off, out.size);
		}
		return out;
	}

	Tensor<F>& reshape(const size_t new_dim[], unsigned new_rank) {
		if (!new_rank)
			if (!rank)
				return *this;
			else
				throw std::invalid_argument("reshape : mismatch sizes.");
		size_t* new_cum_dim = new size_t[new_rank];
		new_cum_dim[new_rank - 1] = 1;
		for (unsigned i = 1; i < new_rank; i++) {
			new_cum_dim[new_rank - i - 1] = new_cum_dim[new_rank - i] * new_dim[new_rank - i];
		}
		size_t new_size = new_cum_dim[0] * new_dim[0];
		if (new_size != size)
			throw std::invalid_argument("reshape : incompatible dimensions: mismatch size.");
		if (!new_size) {
			delete[] new_cum_dim;
			return *this;
		}
		delete[] dim;
		delete[] jump;
		dim = new size_t[new_rank];
		jump = new_cum_dim;
		rank = new_rank;
		memcpy(dim, new_dim, new_rank * sizeof(size_t));
		append_meta("rshp,");
		return *this;
	}

	Tensor<F>& reshape(const std::initializer_list<size_t>& dim_l) {
		return reshape(dim_l.begin(), dim_l.size());
	}

	Tensor<F> operator+(const Tensor<F>& rhs) const {
		if ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t)))
			throw std::invalid_argument(" + : incompatible dimensions.");
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] + rhs.data[i + rhs.off];
		out.append_meta("(" + get_meta() + "+" + rhs.get_meta() + "),");
		return out;
	}

	Tensor<F> operator-(const Tensor<F>& rhs) const {
		if ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t)))
			throw std::invalid_argument(" - : incompatible dimensions.");
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] - rhs.data[i + rhs.off];
		return out;
	}

	Tensor<F> operator*(const Tensor<F>& rhs) const {
		if ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t)))
			throw std::invalid_argument(" * : incompatible dimensions.");
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] * rhs.data[i + rhs.off];
		return out;
	}

	Tensor<F> operator/(const Tensor<F>& rhs) const {
		if ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t)))
			throw std::invalid_argument(" / : incompatible dimensions.");
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] / rhs.data[i + rhs.off];
		return out;
	}

	Tensor<F>& operator+=(Tensor<F>& rhs) {
		if (this != &rhs && ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t))))
			throw std::invalid_argument(" += : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] += rhs.data[i + rhs.off];

		return *this;
	}

	Tensor<F>& operator-=(Tensor<F>& rhs) {
		if (this != &rhs && ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t))))
			throw std::invalid_argument(" -= : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] -= rhs.data[i + rhs.off];

		return *this;

	}

	Tensor<F>& operator*=(Tensor<F>& rhs) {
		if (this != &rhs && ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t))))
			throw std::invalid_argument(" *= : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] *= rhs.data[i + rhs.off];

		return *this;
	}

	Tensor<F> operator*(F rhs) const {
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] * rhs;
		return out;
	}

	Tensor<F> operator+(F rhs) const {
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] + rhs;
		return out;
	}

	Tensor<F> operator/(F rhs) const {
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] / rhs;
		return out;
	}

	Tensor<F> operator-(F rhs) const {
		Tensor<F> out{ dim, rank };
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] - rhs;
		return out;
	}

	Tensor<F>& operator*=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] *= value;

		return *this;
	}

	Tensor<F>& operator+=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] += value;

		return *this;
	}

	Tensor<F>& operator/=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] /= value;

		return *this;
	}

	Tensor<F>& operator-=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] -= value;

		return *this;
	}

	Tensor<F>& operator=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] = value;

		return *this;
	}


	bool is_close_to(const Tensor<F>& rhs, F eps) const {
		if (this == &rhs)
			return true;
		if ((rank != rhs.rank) || memcmp(dim, rhs.dim, rank * sizeof(size_t)))
			throw std::invalid_argument("is_close_to : incompatible dimensions.");
		for (size_t i = 0; i < size; i++) {
			F dif = data[i + off] - rhs.data[i + rhs.off];
			if (std::abs(dif) > std::abs(eps))
				return false;
		}
		return true;
	}


	std::string to_string() const {
		std::string out = "Rank " + std::to_string(rank) + " Tensor, dims: ";
		if (!rank) {
			out += "-\n";
		}
		else {
			for (unsigned r = 0; r < rank; r++)
				out += std::to_string(dim[r]) + ((r != rank - 1) ? "," : "");
		}
		out += "\nmeta: " + std::string(meta) + "\n";

		size_t* index = new size_t[rank]{ 0 };
		size_t i = 0;
		out += (rank == 1) ? "\n\t" : "";
		do {
			if ((rank > 1) && i % dim[rank - 1] == 0) {
				out += "\n";
				if ((rank > 2) && i % (dim[rank - 1] * dim[rank - 2]) == 0) {
					out += ">";
					for (unsigned r = 0; r < rank - 2; r++)
						out += std::to_string(index[r]) + ",";
					out += "\n";
				}
				out += std::to_string(index[rank - 2]) + ":\t";
			}
			out += to_str(data.at(i + off)) + " ";
			i++;
		} while (i < size && inc_index(index));
		out += "\n";
		delete[] index;
		return out;
	}

	unsigned get_rank(void) const {
		return rank;
	}

	size_t get_size(void) const {
		return size;
	}

	void get_dim(size_t cap[]) const {
		memcpy(cap, this->dim, rank * sizeof(size_t));
	}

	Payload<F> get_data() const {
		return data;
	}

	Tensor<F>& apply(F(*func)(F)) {
		for (size_t i = 0; i < size; i++)
			data[i + off] = func(data[i + off]);
		return *this;
	}

	Tensor<F> operator[](size_t i) const {
		if (i >= size || !rank)
			return Tensor<F>();
		if (rank == 1) {
			Tensor<F> out;
			out.rank = 1;
			out.size = 1;
			out.dim = new size_t{ 1 };
			out.jump = new size_t{ 1 };
			out.off = i + off;
			out.data = data;
			return out;
		}

		Tensor<F> out;
		out.rank = rank - 1;
		out.dim = new size_t[out.rank];
		memcpy(out.dim, dim + 1, out.rank * sizeof(size_t));
		out.jump = new size_t[out.rank];
		memcpy(out.jump, jump + 1, out.rank * sizeof(size_t));
		out.size = size / dim[0];
		out.off = off + i * jump[0];
		out.data = data;
		return out;

	}

	Tensor<F> operator[](const std::initializer_list<size_t>& index) {
		if (rank < index.size() || !rank)
			return Tensor<F>();
		for (unsigned r = 0; r < index.size(); r++)
			if (index.begin()[r] >= dim[r])
				return Tensor<F>();

		if (rank == index.size()) {
			size_t pos = off;
			for (unsigned i = 0; i < rank; i++)
				pos += index.begin()[i] * jump[i];
			Tensor<F> out;
			out.rank = 1;
			out.size = 1;
			out.dim = new size_t{ 1 };
			out.jump = new size_t{ 1 };
			out.off = pos;
			out.data = data;
			return out;
		}
		Tensor<F> out;
		out.rank = rank - index.size();
		out.dim = new size_t[out.rank];
		memcpy(out.dim, dim + index.size(), out.rank * sizeof(size_t));
		out.jump = new size_t[out.rank];
		memcpy(out.jump, jump + index.size(), out.rank * sizeof(size_t));
		out.off = off;
		out.size = size;
		for (unsigned i = 0; i < index.size(); i++) {
			out.off += index.begin()[i] * jump[i];
			out.size /= dim[i];
		}
		out.data = data;
		return out;
	}

private:

	inline F& elm(size_t i) {
		throw std::logic_error("elm : not implemented.");
	}

	inline bool check_index(size_t index[]) {
		if (!rank)
			return false;
		for (unsigned r = 0; r < rank; r++)
			if (index[r] >= dim[r]) {
				//throw std::invalid_argument("index: " + std::to_string(r) + " is out of the bound.");
				return false;
			}
		return true;
	}

	inline bool inc_index(size_t index[]) const {
		if (!rank || !size)
			return false;
		unsigned i = rank - 1;
		while (++index[i] == dim[i]) {
			if (i == 0)
				return false;
			index[i--] = 0;
		}
		return true;
	}

	inline size_t index_to_pos(size_t index[]) const {
		size_t pos = off;
		for (unsigned r = 0; r < rank; r++) {
			pos += index[r] * jump[r];	// out = index[r] + cap[r] * out;
		}
		return pos;
	}

	inline void pos_to_index(size_t pos, size_t index[]) const {
		//if (i >= size)
		//	throw std::invalid_argument("Out of bound index.");

		size_t m, rm = (pos - off);
		for (unsigned r = 0; r < rank; r++) {
			m = rm / jump[r];
			rm = rm % jump[r];
			index[r] = m;
		}
	}

public:
	// Tensor contraction
	static Tensor<F> cnt(const Tensor<F>& lt, const Tensor<F>& rt, unsigned li, unsigned ri) {

		Tensor<F> out("cnt(" + lt.get_meta() + ":" + rt.get_meta() + ");");
		if (lt.rank == 0 || rt.rank == 0) {
			return out;
		}
		if (li >= lt.rank || ri >= rt.rank || lt.dim[li] != rt.dim[ri])
			throw std::invalid_argument("cnt : out of rank index or mismatch dimension.");
		size_t out_rank = lt.rank + rt.rank - 2;
		if (out_rank == 0) {
			F val = 0;
			for (size_t k = 0; k < lt.dim[0]; k++)
				val += lt.data[k + lt.off] * rt.data[k + rt.off];
			out.init({ 1 });
			out.data = Payload<F>(1, val);
			return out;
		}

		size_t* out_dim;
		out_dim = new size_t[out_rank];
		comp_index(lt.dim, rt.dim, out_dim, lt.rank, rt.rank, li, ri);
		out.init(out_dim, out_rank);

		out.data = Payload<F>(out.size);

		size_t l_step = lt.jump[li];
		size_t r_step = rt.jump[ri];
		size_t l_jump = (li > 0) ? lt.jump[li - 1] : lt.size;
		size_t r_jump = (ri > 0) ? rt.jump[ri - 1] : rt.size;
		size_t pos = 0;
		for (size_t l_i_off = lt.off; l_i_off < lt.size + lt.off; l_i_off += l_jump) {
			for (size_t l_off = l_i_off; l_off < l_i_off + lt.jump[li]; l_off++) {
				for (size_t r_i_off = rt.off; r_i_off < rt.size + rt.off; r_i_off += r_jump) {
					for (size_t r_off = r_i_off; r_off < r_i_off + rt.jump[ri]; r_off++) {
						out.data[pos] = 0;
						for (size_t l_pos = l_off, r_pos = r_off; l_pos < l_off + l_jump; l_pos += l_step, r_pos += r_step)
							out.data[pos] += lt.data[l_pos] * rt.data[r_pos];
						pos++;
					}
				}
			}
		}
		return out;
	}

	Tensor<F> operator%(const Tensor<F>& rhs) {
		return cnt(*this, rhs, rank - 1, 0);
	}

	void set_meta(const char* meta_str) {
		strncpy(meta, meta_str, 63);
		meta[63] = 0;
	}

	void set_meta(const std::string& meta_str) {
		set_meta(meta_str.c_str());
	}

	std::string get_meta(void) const {
		return std::string(meta);
	}

	void append_meta(const std::string& mt_str) {
		set_meta(get_meta() + mt_str);
	}


	static Tensor<F> range(F start, F end, F jump = 1) {
		if (!jump || (end <= start && jump > 0) || (end >= start && jump < 0))
			return Tensor<F>();
		size_t size = (size_t)ceil((double)(end - start) / jump);
		size_t dim[1]{ size };
		Tensor<F> out(dim, 1);
		for (size_t i = 0; i < out.size; i++)
			out.data[i] = start + i * jump;
		return out;
	}

	static Tensor<F> linspace(F start, F end, size_t size) {
		if (!size || (size == 1 && start != end))
			return Tensor<F>();
		size_t dim[1]{ size };
		return Tensor<F>(dim, 1, Payload<F>(size, start, end));
	}

};