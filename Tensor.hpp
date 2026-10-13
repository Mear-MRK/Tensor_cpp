#pragma once
#include <algorithm>
#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cstring>
#include <cstdio>
#include <initializer_list>
#include <complex>
#include <limits>
#include <utility>
#include <vector>
#include "Payload.hpp"

#ifdef TENSOR_USE_BLAS
#if __has_include(<cblas.h>)
#include <cblas.h>
#elif __has_include(<Accelerate/Accelerate.h>)
#include <Accelerate/Accelerate.h>
#else
#error TENSOR_USE_BLAS requires a CBLAS-compatible header.
#endif
#endif

#ifdef DEBUG
#define LOG(s) std::cout << s << std::endl
#else
#define LOG(s)
#endif

// Combines the uncontracted left and right indices into one output index.
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

// Splits an output index into left and right indices around contracted axes.
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

// Formats a complex value for tensor display.
template<typename F>
inline std::string to_str(std::complex<F> c) {
	return std::to_string(c.real()) + "+" + std::to_string(c.imag()) + "i";
}

// Formats a scalar value for tensor display.
template<typename F>
inline std::string to_str(F f) {
	return std::to_string(f);
}

// Reports whether matrix dimensions fit the integer BLAS interface.
inline bool blas_size_supported(size_t m, size_t k, size_t n) {
	size_t max_size = static_cast<size_t>(std::numeric_limits<int>::max());
	return m <= max_size && k <= max_size && n <= max_size;
}

// Declines BLAS dispatch for unsupported element types.
template<typename F>
inline bool blas_matrix_product(const F*, const F*, F*, size_t, size_t, size_t) {
	return false;
}

#ifdef TENSOR_USE_BLAS
// Multiplies float matrices through the configured CBLAS implementation.
inline bool blas_matrix_product(const float* left, const float* right, float* out,
	size_t m, size_t k, size_t n) {
	if (!blas_size_supported(m, k, n))
		return false;
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		static_cast<int>(m), static_cast<int>(n), static_cast<int>(k), 1.0f,
		left, static_cast<int>(k), right, static_cast<int>(n), 0.0f,
		out, static_cast<int>(n));
	return true;
}

// Multiplies double matrices through the configured CBLAS implementation.
inline bool blas_matrix_product(const double* left, const double* right, double* out,
	size_t m, size_t k, size_t n) {
	if (!blas_size_supported(m, k, n))
		return false;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		static_cast<int>(m), static_cast<int>(n), static_cast<int>(k), 1.0,
		left, static_cast<int>(k), right, static_cast<int>(n), 0.0,
		out, static_cast<int>(n));
	return true;
}

// Multiplies complex-float matrices through the configured CBLAS implementation.
inline bool blas_matrix_product(const std::complex<float>* left,
	const std::complex<float>* right, std::complex<float>* out,
	size_t m, size_t k, size_t n) {
	if (!blas_size_supported(m, k, n))
		return false;
	std::complex<float> alpha{ 1.0f, 0.0f };
	std::complex<float> beta{ 0.0f, 0.0f };
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		static_cast<int>(m), static_cast<int>(n), static_cast<int>(k), &alpha,
		left, static_cast<int>(k), right, static_cast<int>(n), &beta,
		out, static_cast<int>(n));
	return true;
}

// Multiplies complex-double matrices through the configured CBLAS implementation.
inline bool blas_matrix_product(const std::complex<double>* left,
	const std::complex<double>* right, std::complex<double>* out,
	size_t m, size_t k, size_t n) {
	if (!blas_size_supported(m, k, n))
		return false;
	std::complex<double> alpha{ 1.0, 0.0 };
	std::complex<double> beta{ 0.0, 0.0 };
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		static_cast<int>(m), static_cast<int>(n), static_cast<int>(k), &alpha,
		left, static_cast<int>(k), right, static_cast<int>(n), &beta,
		out, static_cast<int>(n));
	return true;
}
#endif

// Stores immutable tensor dimensions and strides with thread-safe shared ownership.
class TensorShape {
private:
	// Owns immutable vectors of dimensions and row-major strides.
	struct Storage {
		std::atomic<size_t> ref_count{ 1 };
		std::vector<size_t> dim;
		std::vector<size_t> jump;

		// Builds checked row-major dimensions and strides.
		Storage(const size_t dimensions[], unsigned rank)
			: dim(dimensions, dimensions + rank), jump(rank) {
			size_t size = 1;
			for (unsigned i = rank; i > 0; i--) {
				jump[i - 1] = size;
				if (dim[i - 1] > std::numeric_limits<size_t>::max() / size)
					throw std::overflow_error("tensor dimensions overflow size_t.");
				size *= dim[i - 1];
			}
		}
	};

	Storage* storage = nullptr;
	unsigned first = 0;
	unsigned rank = 0;
	bool scalar = false;

	// Increments the shared shape reference count when present.
	void inc_ref_count() {
		if (storage)
			storage->ref_count.fetch_add(1, std::memory_order_relaxed);
	}

	// Decrements the shape reference count and frees unreferenced storage.
	void dec_ref_count() {
		if (storage && storage->ref_count.fetch_sub(1, std::memory_order_acq_rel) == 1)
			delete storage;
	}

public:
	// Creates an empty shape.
	TensorShape() = default;

	// Creates a checked row-major shape from dimensions and rank.
	TensorShape(const size_t dimensions[], unsigned t_rank) {
		if (!t_rank)
			return;
		if (!dimensions)
			throw std::invalid_argument("tensor dimensions cannot be null.");
		for (unsigned i = 0; i < t_rank; i++)
			if (!dimensions[i])
				return;
		storage = new Storage(dimensions, t_rank);
		rank = t_rank;
	}

	// Shares immutable shape storage with another shape.
	TensorShape(const TensorShape& oth)
		: storage(oth.storage), first(oth.first), rank(oth.rank), scalar(oth.scalar) {
		inc_ref_count();
	}

	// Replaces this shape with shared immutable shape storage.
	TensorShape& operator=(const TensorShape& rhs) {
		if (this == &rhs)
			return *this;
		dec_ref_count();
		storage = rhs.storage;
		first = rhs.first;
		rank = rhs.rank;
		scalar = rhs.scalar;
		inc_ref_count();
		return *this;
	}

	// Transfers immutable shape storage from another shape.
	TensorShape(TensorShape&& oth) noexcept
		: storage(oth.storage), first(oth.first), rank(oth.rank), scalar(oth.scalar) {
		oth.storage = nullptr;
		oth.first = 0;
		oth.rank = 0;
		oth.scalar = false;
	}

	// Replaces this shape by moving immutable storage from another shape.
	TensorShape& operator=(TensorShape&& rhs) noexcept {
		if (this == &rhs)
			return *this;
		dec_ref_count();
		storage = rhs.storage;
		first = rhs.first;
		rank = rhs.rank;
		scalar = rhs.scalar;
		rhs.storage = nullptr;
		rhs.first = 0;
		rhs.rank = 0;
		rhs.scalar = false;
		return *this;
	}

	// Releases the shared shape storage.
	~TensorShape() {
		dec_ref_count();
	}

	// Returns the number of represented dimensions.
	unsigned get_rank() const {
		return rank;
	}

	// Returns the represented element count.
	size_t get_size() const {
		if (!storage || !rank)
			return 0;
		if (scalar)
			return 1;
		size_t size = 1;
		for (unsigned i = 0; i < rank; i++) {
			if (dimension(i) > std::numeric_limits<size_t>::max() / size)
				throw std::overflow_error("tensor dimensions overflow size_t.");
			size *= dimension(i);
		}
		return size;
	}

	// Returns one represented dimension.
	size_t dimension(unsigned i) const {
		return scalar ? 1 : storage->dim[first + i];
	}

	// Returns one represented row-major stride.
	size_t stride(unsigned i) const {
		return scalar ? 1 : storage->jump[first + i];
	}

	// Reports whether two represented shapes have equal dimensions.
	bool same_dimensions(const TensorShape& rhs) const {
		if (rank != rhs.rank)
			return false;
		for (unsigned i = 0; i < rank; i++)
			if (dimension(i) != rhs.dimension(i))
				return false;
		return true;
	}

	// Copies the represented dimensions into caller-provided storage.
	void copy_dimensions(size_t out[]) const {
		for (unsigned i = 0; i < rank; i++)
			out[i] = dimension(i);
	}

	// Returns a shared shape with leading dimensions removed.
	TensorShape slice(unsigned count) const {
		if (count > rank)
			return TensorShape();
		TensorShape out(*this);
		out.first += count;
		out.rank -= count;
		return out;
	}

	// Returns a shared synthetic rank-one scalar shape.
	TensorShape scalar_shape() const {
		TensorShape out(*this);
		out.rank = 1;
		out.scalar = true;
		return out;
	}
};

// Represents an owned tensor or a shared view over contiguous tensor data.
template<typename F = float>
class Tensor {

private:
	char meta[64] = "";
	TensorShape shape;
	unsigned rank = 0;
	size_t size = 0;
	size_t off = 0;
	Payload<F> data;
	bool view = false;

	// Selects construction with allocated but uninitialized element storage.
	struct UninitializedTag {};

	// Creates a checked rank-one shape for a requested size.
	static TensorShape shape_from_size(size_t t_size) {
		if (!t_size)
			return TensorShape();
		size_t dimensions[]{ t_size };
		return TensorShape(dimensions, 1);
	}

	// Returns shared view storage or an independent owned copy as appropriate.
	static Payload<F> copy_payload(const Tensor<F>& rhs) {
		if (rhs.view)
			return rhs.data;
		return rhs.data.copy(rhs.off, rhs.size);
	}

	// Creates an owned tensor with uninitialized storage for a prepared shape.
	Tensor(TensorShape t_shape, UninitializedTag)
		: shape(std::move(t_shape)), rank(shape.get_rank()), size(shape.get_size()),
		data(size) {}

	// Creates a shared mutable view without allocating shape or payload storage.
	Tensor(const Tensor<F>& source, TensorShape t_shape, size_t t_size, size_t t_off)
		: shape(std::move(t_shape)), rank(shape.get_rank()), size(t_size), off(t_off),
		data(source.data), view(true) {}

	// Exchanges all tensor state without allocation.
	void swap(Tensor<F>& oth) noexcept {
		for (size_t i = 0; i < sizeof(meta); i++)
			std::swap(meta[i], oth.meta[i]);
		std::swap(shape, oth.shape);
		std::swap(rank, oth.rank);
		std::swap(size, oth.size);
		std::swap(off, oth.off);
		std::swap(data, oth.data);
		std::swap(view, oth.view);
	}


	// Returns this tensor to its empty state and releases its resources.
	inline void reset() {
		rank = 0;
		size = 0;
		off = 0;
		view = false;
		shape = TensorShape();
		data.release();
	}

public:
	// Provides read-only access to a shared tensor slice.
	class ConstView {
	private:
		TensorShape shape;
		size_t size = 0;
		size_t off = 0;
		Payload<F> data;

		friend class Tensor<F>;

		// Creates a read-only view over shared tensor storage.
		ConstView(const Payload<F>& t_data, TensorShape t_shape, size_t t_size, size_t t_off)
			: shape(std::move(t_shape)), size(t_size), off(t_off), data(t_data) {}

	public:
		// Creates an empty read-only view.
		ConstView() = default;

		// Returns the number of view dimensions.
		unsigned get_rank() const {
			return shape.get_rank();
		}

		// Returns the number of visible elements.
		size_t get_size() const {
			return size;
		}

		// Copies the view dimensions into caller-provided storage.
		void get_dim(size_t out[]) const {
			shape.copy_dimensions(out);
		}

		// Returns an immutable pointer to the first visible element.
		const F* get_data() const {
			return data.get_arr() ? data.get_arr() + off : nullptr;
		}

		// Returns a checked element selected by a complete index.
		const F& at(const std::initializer_list<size_t>& index) const {
			if (index.size() != shape.get_rank())
				throw std::out_of_range("const tensor view index rank mismatch.");
			size_t pos = off;
			for (unsigned i = 0; i < index.size(); i++) {
				if (index.begin()[i] >= shape.dimension(i))
					throw std::out_of_range("const tensor view index is out of bounds.");
				pos += index.begin()[i] * shape.stride(i);
			}
			return data[pos];
		}

		// Returns a read-only view selected by the first-axis index.
		ConstView operator[](size_t i) const {
			unsigned rank = shape.get_rank();
			if (!rank || i >= shape.dimension(0))
				return ConstView();
			if (rank == 1)
				return ConstView(data, shape.scalar_shape(), 1, off + i);
			return ConstView(data, shape.slice(1), size / shape.dimension(0),
				off + i * shape.stride(0));
		}

		// Returns a read-only view selected by a leading index list.
		ConstView operator[](const std::initializer_list<size_t>& index) const {
			unsigned rank = shape.get_rank();
			if (!rank || index.size() > rank)
				return ConstView();
			for (unsigned i = 0; i < index.size(); i++)
				if (index.begin()[i] >= shape.dimension(i))
					return ConstView();
			size_t out_off = off;
			for (unsigned i = 0; i < index.size(); i++)
				out_off += index.begin()[i] * shape.stride(i);
			if (index.size() == rank)
				return ConstView(data, shape.scalar_shape(), 1, out_off);
			size_t out_size = size;
			for (unsigned i = 0; i < index.size(); i++)
				out_size /= shape.dimension(i);
			return ConstView(data, shape.slice(index.size()), out_size, out_off);
		}
	};

	// Creates an empty tensor.
	Tensor()
	{
		LOG("tns cns, dfl " << this << " : " << meta);
	}

	// Creates an empty tensor with metadata text.
	Tensor(const std::string& meta_str)
	{
		set_meta(meta_str);
		LOG("tns cns, mta " << this << " : " << meta);
	}

	// Creates an uninitialized rank-one tensor of the requested size.
	Tensor(size_t t_size)
		: shape(shape_from_size(t_size)), rank(shape.get_rank()), size(shape.get_size()),
		data(size) {
		LOG("tns cns, siz " << this << " : " << meta);
	}

	// Releases tensor shape metadata and its payload reference.
	~Tensor()
	{
		LOG("tns dst,     " << this << " : " << meta);
	}

	// Copies owned tensors deeply while preserving shared view semantics.
	Tensor(const Tensor<F>& rhs)
		: shape(rhs.shape), rank(rhs.rank), size(rhs.size), off(rhs.view ? rhs.off : 0),
		data(copy_payload(rhs)), view(rhs.view) {
		set_meta(rhs.get_meta() + "Cc;");
		LOG("tns cns, cop " << &rhs << " -> " << this << " : " << meta);
	}
	
	// Copies values into a view or conventionally replaces an owned tensor.
	Tensor<F>& operator=(const Tensor<F>& rhs) {
		if (this == &rhs)
			return *this;
		if (view) {
			if (!shape.same_dimensions(rhs.shape))
				throw std::invalid_argument("copy assignment to view: incompatible dimensions.");
			data.assign(rhs.data, rhs.off, size, off, rhs.size);
			append_meta("]" + rhs.get_meta() + "C=,");
		}
		else {
			Tensor<F> replacement(rhs);
			swap(replacement);
		}

		LOG("tns asn, cop " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	// Transfers tensor shape and payload resources from another tensor.
	Tensor(Tensor<F>&& rhs) noexcept
		: shape(std::move(rhs.shape)), rank(rhs.rank), size(rhs.size), off(rhs.off),
		data(std::move(rhs.data)), view(rhs.view) {
		set_meta(rhs.meta);
		size_t meta_length = strlen(meta);
		strncat(meta, "Mc;", sizeof(meta) - meta_length - 1);

		rhs.rank = 0;
		rhs.size = 0;
		rhs.off = 0;
		rhs.view = false;

		LOG("tns cns, mov " << &rhs << " -> " << this << " : " << meta);

	}
	
	// Copies values into a view or transfers ownership into an owned tensor.
	Tensor<F>& operator=(Tensor<F>&& rhs) {
		if (this == &rhs)
			return *this;
		if (view) {
			if (!shape.same_dimensions(rhs.shape))
				throw std::invalid_argument("move assignment to view: incompatible dimensions.");
			data.assign(rhs.data, rhs.off, size, off, rhs.size);
			rhs.reset();
			append_meta("]M=,");
		}
		else {
			Tensor<F> replacement(std::move(rhs));
			swap(replacement);
		}

		LOG("tns asn, mov " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	// Creates a rank-one tensor that shares an existing payload.
	Tensor(const Payload<F>& t_data)
		: shape(shape_from_size(t_data.get_size())), rank(shape.get_rank()),
		size(shape.get_size()), data(t_data) {
		LOG("tns cns, pyl " << this << " : " << meta);
	}

	// Creates a zero-filled tensor with the requested shape.
	Tensor(const size_t dim[], unsigned t_rank)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size, F{}) {
		set_meta("0;");
		LOG("tns cns, dr0 " << this << " : " << meta);
	}

	// Creates a shaped tensor backed by an existing payload.
	Tensor(const size_t dim[], unsigned t_rank, const Payload<F>& t_data)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(t_data) {
		if (size > t_data.get_size())
			reset();
		set_meta("drP;");
		LOG("tns cns, drP " << this << " : " << meta);
	}

	// Creates a tensor filled with a constant value.
	Tensor(const size_t dim[], unsigned t_rank, F value)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size, value) {
		set_meta(to_str(value) + ";");
		LOG("tns cns, drv " << this << " : " << meta);

	}

	// Creates a tensor by copying values from an array.
	Tensor(const size_t dim[], unsigned t_rank, const F* array)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()),
		data(array, size) {
		set_meta("array;");

		LOG("tns cns, dra " << this << " : " << meta);

	}

	// Creates a tensor filled by repeatedly calling a generator function.
	Tensor(const size_t dim[], unsigned t_rank, F(*rnd)(void))
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size) {
		if (size) {
			for (size_t i = 0; i < size; i++)
				data[i] = rnd();
		}
		set_meta("rnd;");

		LOG("tns cns, drr " << this << " : " << meta);

	}

	// Creates a tensor filled with evenly spaced values.
	Tensor(const size_t dim[], unsigned t_rank, F start, F end)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()),
		data(size, start, end) {
		set_meta("lin;");
		LOG("tns cns, drl " << this << " : " << meta);
	}

	// Creates a zero-filled tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l) : Tensor(dim_l.begin(), dim_l.size()) {}

	// Creates a constant-filled tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l, F value) : Tensor(dim_l.begin(), dim_l.size(), value) {
	}

	// Creates a tensor from a dimension list and existing payload.
	Tensor(const std::initializer_list<size_t>& dim_l, const Payload<F>& data) : Tensor(dim_l.begin(), dim_l.size(), data) {
	}

	// Creates a generated tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l, F(*rnd)(void)) : Tensor(dim_l.begin(), dim_l.size(), rnd) {
	}

	// Creates a tensor from a dimension list and copied array values.
	Tensor(const std::initializer_list<size_t>& dim_l, const F arr[]) : Tensor(dim_l.begin(), dim_l.size(), arr) {
	}

	// Creates a linearly spaced tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l, F start, F end) : Tensor(dim_l.begin(), dim_l.size(), start, end)
	{}

	// Returns an independent deep copy of this tensor.
	Tensor<F> copy() const {
		Tensor<F> out(shape, UninitializedTag{});
		if (size)
			std::copy_n(data.get_arr() + off, size, out.data.get_arr());
		return out;
	}

	// Changes the tensor shape while preserving its element count.
	Tensor<F>& reshape(const size_t new_dim[], unsigned new_rank) & {
		if (!new_rank) {
			if (!rank)
				return *this;
			else
				throw std::invalid_argument("reshape : mismatch sizes.");
		}
		TensorShape new_shape(new_dim, new_rank);
		size_t new_size = new_shape.get_size();
		if (new_size != size) {
			throw std::invalid_argument("reshape : incompatible dimensions: mismatch size.");
		}
		if (!new_size)
			return *this;
		std::string reshaped_meta = get_meta() + "rshp,";
		shape = std::move(new_shape);
		rank = new_rank;
		set_meta(reshaped_meta);
		return *this;
	}

	// Reshapes a temporary tensor and returns it without copying its payload.
	Tensor<F> reshape(const size_t new_dim[], unsigned new_rank) && {
		static_cast<Tensor<F>&>(*this).reshape(new_dim, new_rank);
		return std::move(*this);
	}

	// Changes the tensor shape using a dimension list.
	Tensor<F>& reshape(const std::initializer_list<size_t>& dim_l) & {
		return reshape(dim_l.begin(), dim_l.size());
	}

	// Reshapes a temporary tensor from a dimension list without copying its payload.
	Tensor<F> reshape(const std::initializer_list<size_t>& dim_l) && {
		return std::move(*this).reshape(dim_l.begin(), dim_l.size());
	}

	// Returns the element-wise sum of equally shaped tensors.
	Tensor<F> operator+(const Tensor<F>& rhs) const {
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" + : incompatible dimensions.");
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] + rhs.data[i + rhs.off];
		out.append_meta("(" + get_meta() + "+" + rhs.get_meta() + "),");
		return out;
	}

	// Returns the element-wise difference of equally shaped tensors.
	Tensor<F> operator-(const Tensor<F>& rhs) const {
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" - : incompatible dimensions.");
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] - rhs.data[i + rhs.off];
		return out;
	}

	// Returns the element-wise product of equally shaped tensors.
	Tensor<F> operator*(const Tensor<F>& rhs) const {
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" * : incompatible dimensions.");
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] * rhs.data[i + rhs.off];
		return out;
	}

	// Returns the element-wise quotient of equally shaped tensors.
	Tensor<F> operator/(const Tensor<F>& rhs) const {
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" / : incompatible dimensions.");
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] / rhs.data[i + rhs.off];
		return out;
	}

	// Adds an equally shaped tensor to this tensor in place.
	Tensor<F>& operator+=(const Tensor<F>& rhs) {
		if (this != &rhs && !shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" += : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] += rhs.data[i + rhs.off];

		return *this;
	}

	// Subtracts an equally shaped tensor from this tensor in place.
	Tensor<F>& operator-=(const Tensor<F>& rhs) {
		if (this != &rhs && !shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" -= : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] -= rhs.data[i + rhs.off];

		return *this;

	}

	// Multiplies this tensor by an equally shaped tensor in place.
	Tensor<F>& operator*=(const Tensor<F>& rhs) {
		if (this != &rhs && !shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(" *= : incompatible dimensions.");

		for (size_t i = 0; i < size; i++)
			data[i + off] *= rhs.data[i + rhs.off];

		return *this;
	}

	// Returns a tensor with every element multiplied by a scalar.
	Tensor<F> operator*(F rhs) const {
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] * rhs;
		return out;
	}

	// Returns a tensor with a scalar added to every element.
	Tensor<F> operator+(F rhs) const {
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] + rhs;
		return out;
	}

	// Returns a tensor with every element divided by a scalar.
	Tensor<F> operator/(F rhs) const {
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] / rhs;
		return out;
	}

	// Returns a tensor with a scalar subtracted from every element.
	Tensor<F> operator-(F rhs) const {
		Tensor<F> out(shape, UninitializedTag{});
		for (size_t i = 0; i < size; i++)
			out.data[i] = data[i + off] - rhs;
		return out;
	}

	// Multiplies every element by a scalar in place.
	Tensor<F>& operator*=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] *= value;

		return *this;
	}

	// Adds a scalar to every element in place.
	Tensor<F>& operator+=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] += value;

		return *this;
	}

	// Divides every element by a scalar in place.
	Tensor<F>& operator/=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] /= value;

		return *this;
	}

	// Subtracts a scalar from every element in place.
	Tensor<F>& operator-=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] -= value;

		return *this;
	}

	// Assigns one scalar value to every tensor element.
	Tensor<F>& operator=(F value) {
		for (size_t i = 0; i < size; i++)
			data[i + off] = value;

		return *this;
	}


	// Reports whether equally shaped tensors match within a tolerance.
	bool is_close_to(const Tensor<F>& rhs, F eps) const {
		if (this == &rhs)
			return true;
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument("is_close_to : incompatible dimensions.");
		for (size_t i = 0; i < size; i++) {
			F dif = data[i + off] - rhs.data[i + rhs.off];
			if (std::abs(dif) > std::abs(eps))
				return false;
		}
		return true;
	}


	// Formats tensor shape, metadata, and values as text.
	std::string to_string() const {
		std::string out = "Rank " + std::to_string(rank) + " Tensor, dims: ";
		if (!rank) {
			out += "-\n";
		}
		else {
			for (unsigned r = 0; r < rank; r++)
				out += std::to_string(shape.dimension(r)) + ((r != rank - 1) ? "," : "");
		}
		out += "\nmeta: " + std::string(meta) + "\n";
		if (!size)
			return out;

		std::vector<size_t> index(rank, 0);
		size_t i = 0;
		out += (rank == 1) ? "\n\t" : "";
		do {
			if ((rank > 1) && i % shape.dimension(rank - 1) == 0) {
				out += "\n";
				if ((rank > 2) && i % (shape.dimension(rank - 1) * shape.dimension(rank - 2)) == 0) {
					out += ">";
					for (unsigned r = 0; r < rank - 2; r++)
						out += std::to_string(index[r]) + ",";
					out += "\n";
				}
				out += std::to_string(index[rank - 2]) + ":\t";
			}
			out += to_str(data[i + off]) + " ";
			i++;
		} while (i < size && inc_index(index.data()));
		out += "\n";
		return out;
	}

	// Returns the number of tensor dimensions.
	unsigned get_rank(void) const {
		return rank;
	}

	// Returns the total number of tensor elements.
	size_t get_size(void) const {
		return size;
	}

	// Copies the tensor dimensions into a caller-provided array.
	void get_dim(size_t cap[]) const {
		shape.copy_dimensions(cap);
	}

	// Returns a shared mutable reference to the underlying payload.
	Payload<F> get_data() {
		return data;
	}

	// Returns an immutable pointer to the first visible element.
	const F* get_data() const {
		return data.get_arr() ? data.get_arr() + off : nullptr;
	}

	// Applies a unary function to every tensor element in place.
	Tensor<F>& apply(F(*func)(F)) {
		for (size_t i = 0; i < size; i++)
			data[i + off] = func(data[i + off]);
		return *this;
	}

	// Returns a shared view selected by the first-axis index.
	Tensor<F> operator[](size_t i) {
		if (!rank || i >= shape.dimension(0))
			return Tensor<F>();
		if (rank == 1)
			return Tensor<F>(*this, shape.scalar_shape(), 1, off + i);
		return Tensor<F>(*this, shape.slice(1), size / shape.dimension(0),
			off + i * shape.stride(0));
	}

	// Returns a read-only shared view selected from a const tensor.
	ConstView operator[](size_t i) const {
		return ConstView(data, shape, size, off)[i];
	}

	// Returns a shared view selected by a leading index list.
	Tensor<F> operator[](const std::initializer_list<size_t>& index) {
		if (rank < index.size() || !rank)
			return Tensor<F>();
		for (unsigned r = 0; r < index.size(); r++)
			if (index.begin()[r] >= shape.dimension(r))
				return Tensor<F>();

		if (rank == index.size()) {
			size_t pos = off;
			for (unsigned i = 0; i < rank; i++)
				pos += index.begin()[i] * shape.stride(i);
			return Tensor<F>(*this, shape.scalar_shape(), 1, pos);
		}
		size_t out_off = off;
		size_t out_size = size;
		for (unsigned i = 0; i < index.size(); i++) {
			out_off += index.begin()[i] * shape.stride(i);
			out_size /= shape.dimension(i);
		}
		return Tensor<F>(*this, shape.slice(index.size()), out_size, out_off);
	}

	// Returns a read-only shared view selected from a const tensor.
	ConstView operator[](const std::initializer_list<size_t>& index) const {
		return ConstView(data, shape, size, off)[index];
	}

private:

	// Rejects direct element lookup because it is not implemented.
	inline F& elm(size_t) {
		throw std::logic_error("elm : not implemented.");
	}

	// Reports whether every coordinate lies within the tensor shape.
	inline bool check_index(size_t index[]) {
		if (!rank)
			return false;
		for (unsigned r = 0; r < rank; r++)
			if (index[r] >= shape.dimension(r)) {
				//throw std::invalid_argument("index: " + std::to_string(r) + " is out of the bound.");
				return false;
			}
		return true;
	}

	// Advances a multidimensional index in row-major order.
	inline bool inc_index(size_t index[]) const {
		if (!rank || !size)
			return false;
		unsigned i = rank - 1;
		while (++index[i] == shape.dimension(i)) {
			if (i == 0)
				return false;
			index[i--] = 0;
		}
		return true;
	}

	// Converts a multidimensional index to its payload position.
	inline size_t index_to_pos(size_t index[]) const {
		size_t pos = off;
		for (unsigned r = 0; r < rank; r++) {
			pos += index[r] * shape.stride(r);
		}
		return pos;
	}

	// Converts a payload position to a multidimensional index.
	inline void pos_to_index(size_t pos, size_t index[]) const {
		//if (i >= size)
		//	throw std::invalid_argument("Out of bound index.");

		size_t m, rm = (pos - off);
		for (unsigned r = 0; r < rank; r++) {
			m = rm / shape.stride(r);
			rm = rm % shape.stride(r);
			index[r] = m;
		}
	}

public:
	// Contracts one axis from each tensor through packed matrix multiplication.
	static Tensor<F> cnt(const Tensor<F>& lt, const Tensor<F>& rt, unsigned li, unsigned ri) {
		std::string out_meta = "cnt(" + lt.get_meta() + ":" + rt.get_meta() + ");";
		if (lt.rank == 0 || rt.rank == 0)
			return Tensor<F>(out_meta);
		if (li >= lt.rank || ri >= rt.rank
			|| lt.shape.dimension(li) != rt.shape.dimension(ri))
			throw std::invalid_argument("cnt : out of rank index or mismatch dimension.");
		unsigned out_rank = lt.rank + rt.rank - 2;
		if (out_rank == 0) {
			F val{};
			for (size_t k = 0; k < lt.shape.dimension(0); k++)
				val += lt.data[k + lt.off] * rt.data[k + rt.off];
			Tensor<F> out({ 1 }, val);
			out.set_meta(out_meta);
			return out;
		}

		std::vector<size_t> out_dim;
		out_dim.reserve(out_rank);
		for (unsigned i = 0; i < lt.rank; i++)
			if (i != li)
				out_dim.push_back(lt.shape.dimension(i));
		for (unsigned i = 0; i < rt.rank; i++)
			if (i != ri)
				out_dim.push_back(rt.shape.dimension(i));
		Tensor<F> out(TensorShape(out_dim.data(), out_rank), UninitializedTag{});
		out.set_meta(out_meta);

		size_t k_size = lt.shape.dimension(li);
		size_t m_size = lt.size / k_size;
		size_t n_size = rt.size / k_size;

		Payload<F> l_matrix;
		const F* l_arr = lt.data.get_arr() + lt.off;
		if (li != lt.rank - 1) {
			// Move the contracted left index to the last matrix index.
			l_matrix = Payload<F>(lt.size);
			size_t suffix_size = lt.shape.stride(li);
			for (size_t pos = 0; pos < lt.size; pos++) {
				size_t k = (pos / suffix_size) % k_size;
				size_t row = (pos / (k_size * suffix_size)) * suffix_size + pos % suffix_size;
				l_matrix[row * k_size + k] = lt.data[pos + lt.off];
			}
			l_arr = l_matrix.get_arr();
		}

		Payload<F> r_matrix;
		const F* r_arr = rt.data.get_arr() + rt.off;
		if (ri != 0) {
			// Move the contracted right index to the first matrix index.
			r_matrix = Payload<F>(rt.size);
			size_t suffix_size = rt.shape.stride(ri);
			for (size_t pos = 0; pos < rt.size; pos++) {
				size_t k = (pos / suffix_size) % k_size;
				size_t col = (pos / (k_size * suffix_size)) * suffix_size + pos % suffix_size;
				r_matrix[k * n_size + col] = rt.data[pos + rt.off];
			}
			r_arr = r_matrix.get_arr();
		}

		F* out_arr = out.data.get_arr();
		if (!blas_matrix_product(l_arr, r_arr, out_arr, m_size, k_size, n_size)) {
			std::fill_n(out_arr, out.size, F{});
			for (size_t m = 0; m < m_size; m++) {
				for (size_t k = 0; k < k_size; k++) {
					F l_value = l_arr[m * k_size + k];
					for (size_t n = 0; n < n_size; n++)
						out_arr[m * n_size + n] += l_value * r_arr[k * n_size + n];
				}
			}
		}
		return out;
	}

	// Contracts this tensor's last axis with another tensor's first axis.
	Tensor<F> operator%(const Tensor<F>& rhs) const {
		return cnt(*this, rhs, rank - 1, 0);
	}

	// Replaces tensor metadata from a C string.
	void set_meta(const char* meta_str) {
		strncpy(meta, meta_str, 63);
		meta[63] = 0;
	}

	// Replaces tensor metadata from a string.
	void set_meta(const std::string& meta_str) {
		set_meta(meta_str.c_str());
	}

	// Returns the tensor metadata text.
	std::string get_meta(void) const {
		return std::string(meta);
	}

	// Appends text to the tensor metadata.
	void append_meta(const std::string& mt_str) {
		set_meta(get_meta() + mt_str);
	}


	// Creates a rank-one tensor containing a half-open arithmetic sequence.
	static Tensor<F> range(F start, F end, F jump = 1) {
		if (!jump || (end <= start && jump > 0) || (end >= start && jump < 0))
			return Tensor<F>();
		long double count = std::ceil((static_cast<long double>(end)
			- static_cast<long double>(start)) / static_cast<long double>(jump));
		long double size_limit = std::ldexp(1.0L, std::numeric_limits<size_t>::digits);
		if (!std::isfinite(count) || count >= size_limit)
			throw std::overflow_error("range size overflows size_t.");
		size_t size = static_cast<size_t>(count);
		Tensor<F> out(size);
		out.set_meta("0;");
		for (size_t i = 0; i < out.size; i++)
			out.data[i] = start + i * jump;
		return out;
	}

	// Creates a rank-one tensor of evenly spaced endpoint-inclusive values.
	static Tensor<F> linspace(F start, F end, size_t size) {
		if (!size || (size == 1 && start != end))
			return Tensor<F>();
		size_t dim[1]{ size };
		return Tensor<F>(dim, 1, start, end);
	}

};
