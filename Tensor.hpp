#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <initializer_list>
#include <complex>
#include <limits>
#include <memory>
#include <type_traits>
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

// Stores immutable tensor dimensions and strides with standard shared ownership.
class TensorShape {
private:
	// Owns immutable vectors of dimensions and row-major strides.
	struct Storage {
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

	std::shared_ptr<const Storage> storage;
	unsigned first = 0;
	unsigned rank = 0;
	bool scalar = false;

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
		storage = std::make_shared<Storage>(dimensions, t_rank);
		rank = t_rank;
	}

	// Shares immutable shape storage with another shape.
	TensorShape(const TensorShape&) = default;

	// Replaces this shape with shared immutable shape storage.
	TensorShape& operator=(const TensorShape&) = default;

	// Transfers immutable shape storage from another shape.
	TensorShape(TensorShape&& oth) noexcept
		: storage(std::move(oth.storage)), first(oth.first), rank(oth.rank),
		scalar(oth.scalar) {
		oth.first = 0;
		oth.rank = 0;
		oth.scalar = false;
	}

	// Replaces this shape by moving immutable storage from another shape.
	TensorShape& operator=(TensorShape&& rhs) noexcept {
		if (this == &rhs)
			return *this;
		storage = std::move(rhs.storage);
		first = rhs.first;
		rank = rhs.rank;
		scalar = rhs.scalar;
		rhs.first = 0;
		rhs.rank = 0;
		rhs.scalar = false;
		return *this;
	}

	// Releases the standard shared shape reference.
	~TensorShape() = default;

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

	// Returns the represented dimensions in a capacity-safe container.
	std::vector<size_t> dimensions() const {
		std::vector<size_t> out;
		out.reserve(rank);
		for (unsigned i = 0; i < rank; i++)
			out.push_back(dimension(i));
		return out;
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

// Represents an independently owned contiguous tensor.
template<typename F = float>
class Tensor {

private:
	std::string meta;
	TensorShape shape;
	unsigned rank = 0;
	size_t size = 0;
	Payload<F> data;

	// Selects construction with allocated but uninitialized element storage.
	struct UninitializedTag {};

	// Creates a checked rank-one shape for a requested size.
	static TensorShape shape_from_size(size_t t_size) {
		if (!t_size)
			return TensorShape();
		size_t dimensions[]{ t_size };
		return TensorShape(dimensions, 1);
	}

	// Creates an owned tensor with uninitialized storage for a prepared shape.
	Tensor(TensorShape t_shape, UninitializedTag)
		: shape(std::move(t_shape)), rank(shape.get_rank()), size(shape.get_size()),
		data(size) {}

	// Exchanges all tensor state without allocation.
	void swap(Tensor<F>& oth) noexcept {
		std::swap(meta, oth.meta);
		std::swap(shape, oth.shape);
		std::swap(rank, oth.rank);
		std::swap(size, oth.size);
		std::swap(data, oth.data);
	}


	// Returns this tensor to its empty state and releases its resources.
	inline void reset() {
		rank = 0;
		size = 0;
		meta.clear();
		shape = TensorShape();
		data.release();
	}

	// Validates that an operand has the same tensor dimensions.
	void require_same_shape(const Tensor<F>& rhs, const char* operation) const {
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument(std::string(operation) + " : incompatible dimensions.");
	}

	// Creates an element-wise result from two equally shaped tensors.
	template<typename Operation>
	Tensor<F> binary_transform(const Tensor<F>& rhs, const char* operation,
		Operation&& transform) const {
		require_same_shape(rhs, operation);
		Tensor<F> out(shape, UninitializedTag{});
		F* out_data = out.data.get_arr();
		const F* left_data = data.get_arr();
		const F* right_data = rhs.data.get_arr();
		for (size_t i = 0; i < size; i++)
			out_data[i] = std::invoke(transform, left_data[i], right_data[i]);
		return out;
	}

	// Creates an element-wise result from this tensor and a scalar.
	template<typename Operation>
	Tensor<F> scalar_transform(F value, Operation&& transform) const {
		Tensor<F> out(shape, UninitializedTag{});
		F* out_data = out.data.get_arr();
		const F* source_data = data.get_arr();
		for (size_t i = 0; i < size; i++)
			out_data[i] = std::invoke(transform, source_data[i], value);
		return out;
	}

	// Applies an in-place operation with an equally shaped tensor.
	template<typename Operation>
	Tensor<F>& binary_compound(const Tensor<F>& rhs, const char* operation,
		Operation&& transform) {
		require_same_shape(rhs, operation);
		F* left_data = data.get_arr();
		const F* right_data = rhs.data.get_arr();
		for (size_t i = 0; i < size; i++)
			std::invoke(transform, left_data[i], right_data[i]);
		return *this;
	}

	// Applies an in-place operation with a scalar.
	template<typename Operation>
	Tensor<F>& scalar_compound(F value, Operation&& transform) {
		F* values = data.get_arr();
		for (size_t i = 0; i < size; i++)
			std::invoke(transform, values[i], value);
		return *this;
	}

public:
	class View;

	// Provides read-only access to a shared tensor slice.
	class ConstView {
	private:
		TensorShape shape;
		size_t size = 0;
		size_t off = 0;
		Payload<F> data;

		friend class Tensor<F>;
		friend class View;

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

		// Returns all view dimensions in a capacity-safe container.
		std::vector<size_t> dimensions() const {
			return shape.dimensions();
		}

		// Returns one checked view dimension.
		size_t dimension(unsigned i) const {
			if (i >= shape.get_rank())
				throw std::out_of_range("const tensor view dimension is out of bounds.");
			return shape.dimension(i);
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

	// Provides mutable access to a shared tensor slice without owning tensor semantics.
	class View {
	private:
		TensorShape shape;
		size_t size = 0;
		size_t off = 0;
		Payload<F> data;

		friend class Tensor<F>;

		// Creates a mutable view over shared tensor storage.
		View(const Payload<F>& t_data, TensorShape t_shape, size_t t_size, size_t t_off)
			: shape(std::move(t_shape)), size(t_size), off(t_off), data(t_data) {}

		// Copies values from a compatible payload range into this view.
		View& assign_values(const Payload<F>& source, const TensorShape& source_shape,
			size_t source_size, size_t source_off) {
			if (!shape.same_dimensions(source_shape))
				throw std::invalid_argument("assignment to tensor view: incompatible dimensions.");
			data.assign(source, source_off, size, off, source_size);
			return *this;
		}

	public:
		// Creates an empty mutable view.
		View() = default;

		// Creates another handle to the same mutable view.
		View(const View&) = default;

		// Creates another handle to the same mutable view from a temporary handle.
		View(View&&) noexcept = default;

		// Copies values from another view while preserving this view's binding.
		View& operator=(const View& rhs) {
			if (this == &rhs)
				return *this;
			return assign_values(rhs.data, rhs.shape, rhs.size, rhs.off);
		}

		// Copies values from a temporary view while preserving this view's binding.
		View& operator=(View&& rhs) {
			return assign_values(rhs.data, rhs.shape, rhs.size, rhs.off);
		}

		// Copies values from an owned tensor into this view.
		View& operator=(const Tensor<F>& rhs) {
			return assign_values(rhs.data, rhs.shape, rhs.size, 0);
		}

		// Copies values from a read-only view into this view.
		View& operator=(const ConstView& rhs) {
			return assign_values(rhs.data, rhs.shape, rhs.size, rhs.off);
		}

		// Assigns one scalar value to every visible element.
		View& operator=(F value) {
			F* values = data.get_arr();
			if (values)
				values += off;
			for (size_t i = 0; i < size; i++)
				values[i] = value;
			return *this;
		}

		// Returns the number of view dimensions.
		unsigned get_rank() const {
			return shape.get_rank();
		}

		// Returns the number of visible elements.
		size_t get_size() const {
			return size;
		}

		// Returns all view dimensions in a capacity-safe container.
		std::vector<size_t> dimensions() const {
			return shape.dimensions();
		}

		// Returns one checked view dimension.
		size_t dimension(unsigned i) const {
			if (i >= shape.get_rank())
				throw std::out_of_range("tensor view dimension is out of bounds.");
			return shape.dimension(i);
		}

		// Returns a mutable pointer to the first visible element.
		F* get_data() {
			return data.get_arr() ? data.get_arr() + off : nullptr;
		}

		// Returns an immutable pointer to the first visible element.
		const F* get_data() const {
			return data.get_arr() ? data.get_arr() + off : nullptr;
		}

		// Returns a checked mutable element selected by a complete index.
		F& at(const std::initializer_list<size_t>& index) {
			if (index.size() != shape.get_rank())
				throw std::out_of_range("tensor view index rank mismatch.");
			size_t pos = off;
			for (unsigned i = 0; i < index.size(); i++) {
				if (index.begin()[i] >= shape.dimension(i))
					throw std::out_of_range("tensor view index is out of bounds.");
				pos += index.begin()[i] * shape.stride(i);
			}
			return data[pos];
		}

		// Returns a checked immutable element selected by a complete index.
		const F& at(const std::initializer_list<size_t>& index) const {
			return ConstView(data, shape, size, off).at(index);
		}

		// Returns a mutable subview selected by the first-axis index.
		View operator[](size_t i) {
			unsigned rank = shape.get_rank();
			if (!rank || i >= shape.dimension(0))
				return View();
			if (rank == 1)
				return View(data, shape.scalar_shape(), 1, off + i);
			return View(data, shape.slice(1), size / shape.dimension(0),
				off + i * shape.stride(0));
		}

		// Returns a read-only subview selected from a const view handle.
		ConstView operator[](size_t i) const {
			return ConstView(data, shape, size, off)[i];
		}

		// Returns a mutable subview selected by a leading index list.
		View operator[](const std::initializer_list<size_t>& index) {
			unsigned rank = shape.get_rank();
			if (!rank || index.size() > rank)
				return View();
			for (unsigned i = 0; i < index.size(); i++)
				if (index.begin()[i] >= shape.dimension(i))
					return View();
			size_t out_off = off;
			for (unsigned i = 0; i < index.size(); i++)
				out_off += index.begin()[i] * shape.stride(i);
			if (index.size() == rank)
				return View(data, shape.scalar_shape(), 1, out_off);
			size_t out_size = size;
			for (unsigned i = 0; i < index.size(); i++)
				out_size /= shape.dimension(i);
			return View(data, shape.slice(index.size()), out_size, out_off);
		}

		// Returns a read-only subview selected from a const view handle.
		ConstView operator[](const std::initializer_list<size_t>& index) const {
			return ConstView(data, shape, size, off)[index];
		}

		// Applies any compatible callable to every visible element in place.
		template<typename Function,
			typename std::enable_if<std::is_invocable_r<F, Function&, F>::value, int>::type = 0>
		View& apply(Function&& function) {
			F* values = data.get_arr();
			if (values)
				values += off;
			for (size_t i = 0; i < size; i++)
				values[i] = std::invoke(function, values[i]);
			return *this;
		}

		// Returns an independent owned copy of the visible values.
		Tensor<F> copy() const {
			return Tensor<F>(*this);
		}

		// Formats the visible tensor shape and values as text.
		std::string to_string() const {
			return copy().to_string();
		}
	};

	// Creates an empty tensor.
	Tensor()
	{
		TENSOR_LOG("tns cns, dfl " << this << " : " << meta);
	}

	// Creates an empty tensor with metadata text.
	Tensor(const std::string& meta_str)
	{
		set_meta(meta_str);
		TENSOR_LOG("tns cns, mta " << this << " : " << meta);
	}

	// Creates an uninitialized rank-one tensor of the requested size.
	Tensor(size_t t_size)
		: shape(shape_from_size(t_size)), rank(shape.get_rank()), size(shape.get_size()),
		data(size) {
		TENSOR_LOG("tns cns, siz " << this << " : " << meta);
	}

	// Releases tensor shape metadata and its payload reference.
	~Tensor()
	{
		TENSOR_LOG("tns dst,     " << this << " : " << meta);
	}

	// Creates an independent deep copy of an owned tensor.
	Tensor(const Tensor<F>& rhs)
		: meta(rhs.meta + "Cc;"), shape(rhs.shape), rank(rhs.rank), size(rhs.size),
		data(rhs.data.copy(0, rhs.size)) {
		TENSOR_LOG("tns cns, cop " << &rhs << " -> " << this << " : " << meta);
	}

	// Creates an independent owned tensor from a mutable view.
	Tensor(const View& source)
		: shape(source.shape), rank(shape.get_rank()), size(source.size),
		data(source.data.copy(source.off, source.size)) {
		TENSOR_LOG("tns cns, viw " << &source << " -> " << this << " : " << meta);
	}

	// Creates an independent owned tensor from a read-only view.
	Tensor(const ConstView& source)
		: shape(source.shape), rank(shape.get_rank()), size(source.size),
		data(source.data.copy(source.off, source.size)) {
		TENSOR_LOG("tns cns, cvw " << &source << " -> " << this << " : " << meta);
	}
	
	// Replaces this owned tensor with an independent copy.
	Tensor<F>& operator=(const Tensor<F>& rhs) {
		if (this == &rhs)
			return *this;
		Tensor<F> replacement(rhs);
		swap(replacement);

		TENSOR_LOG("tns asn, cop " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	// Transfers tensor shape and payload resources from another tensor.
	Tensor(Tensor<F>&& rhs) noexcept
		: meta(std::move(rhs.meta)), shape(std::move(rhs.shape)), rank(rhs.rank),
		size(rhs.size), data(std::move(rhs.data)) {

		rhs.rank = 0;
		rhs.size = 0;

		TENSOR_LOG("tns cns, mov " << &rhs << " -> " << this << " : " << meta);

	}
	
	// Transfers complete tensor ownership into this tensor.
	Tensor<F>& operator=(Tensor<F>&& rhs) noexcept {
		if (this == &rhs)
			return *this;
		Tensor<F> replacement(std::move(rhs));
		swap(replacement);

		TENSOR_LOG("tns asn, mov " << &rhs << " -> " << this << " : " << meta);

		return *this;
	}

	// Creates a rank-one tensor that shares an existing payload.
	Tensor(const Payload<F>& t_data)
		: shape(shape_from_size(t_data.get_size())), rank(shape.get_rank()),
		size(shape.get_size()), data(t_data) {
		TENSOR_LOG("tns cns, pyl " << this << " : " << meta);
	}

	// Creates a zero-filled tensor with the requested shape.
	Tensor(const size_t dim[], unsigned t_rank)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size, F{}) {
		set_meta("0;");
		TENSOR_LOG("tns cns, dr0 " << this << " : " << meta);
	}

	// Creates a shaped tensor backed by an existing payload.
	Tensor(const size_t dim[], unsigned t_rank, const Payload<F>& t_data)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(t_data) {
		if (size > t_data.get_size())
			reset();
		set_meta("drP;");
		TENSOR_LOG("tns cns, drP " << this << " : " << meta);
	}

	// Creates a tensor filled with a constant value.
	Tensor(const size_t dim[], unsigned t_rank, F value)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size, value) {
		set_meta(to_str(value) + ";");
		TENSOR_LOG("tns cns, drv " << this << " : " << meta);

	}

	// Creates a tensor by copying values from an array.
	Tensor(const size_t dim[], unsigned t_rank, const F* array)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()),
		data(array, size) {
		set_meta("array;");

		TENSOR_LOG("tns cns, dra " << this << " : " << meta);

	}

	// Creates a tensor filled by repeatedly invoking any compatible generator.
	template<typename Generator,
		typename std::enable_if<std::is_invocable_r<F, Generator&>::value, int>::type = 0>
	Tensor(const size_t dim[], unsigned t_rank, Generator&& generator)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()), data(size) {
		F* values = data.get_arr();
		for (size_t i = 0; i < size; i++)
			values[i] = std::invoke(generator);
		set_meta("rnd;");

		TENSOR_LOG("tns cns, drr " << this << " : " << meta);

	}

	// Creates a tensor filled with evenly spaced values.
	Tensor(const size_t dim[], unsigned t_rank, F start, F end)
		: shape(dim, t_rank), rank(shape.get_rank()), size(shape.get_size()),
		data(size, start, end) {
		set_meta("lin;");
		TENSOR_LOG("tns cns, drl " << this << " : " << meta);
	}

	// Creates a zero-filled tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l) : Tensor(dim_l.begin(), dim_l.size()) {}

	// Creates a constant-filled tensor from a dimension list.
	Tensor(const std::initializer_list<size_t>& dim_l, F value) : Tensor(dim_l.begin(), dim_l.size(), value) {
	}

	// Creates a tensor from a dimension list and existing payload.
	Tensor(const std::initializer_list<size_t>& dim_l, const Payload<F>& data) : Tensor(dim_l.begin(), dim_l.size(), data) {
	}

	// Creates a generated tensor from a dimension list and compatible callable.
	template<typename Generator,
		typename std::enable_if<std::is_invocable_r<F, Generator&>::value, int>::type = 0>
	Tensor(const std::initializer_list<size_t>& dim_l, Generator&& generator)
		: Tensor(dim_l.begin(), dim_l.size(), std::forward<Generator>(generator)) {
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
			std::copy_n(data.get_arr(), size, out.data.get_arr());
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
		Tensor<F> out = binary_transform(rhs, "+",
			[](const F& left, const F& right) { return left + right; });
		out.append_meta("(" + get_meta() + "+" + rhs.get_meta() + "),");
		return out;
	}

	// Returns the element-wise difference of equally shaped tensors.
	Tensor<F> operator-(const Tensor<F>& rhs) const {
		return binary_transform(rhs, "-",
			[](const F& left, const F& right) { return left - right; });
	}

	// Returns the element-wise product of equally shaped tensors.
	Tensor<F> operator*(const Tensor<F>& rhs) const {
		return binary_transform(rhs, "*",
			[](const F& left, const F& right) { return left * right; });
	}

	// Returns the element-wise quotient of equally shaped tensors.
	Tensor<F> operator/(const Tensor<F>& rhs) const {
		return binary_transform(rhs, "/",
			[](const F& left, const F& right) { return left / right; });
	}

	// Adds an equally shaped tensor to this tensor in place.
	Tensor<F>& operator+=(const Tensor<F>& rhs) {
		return binary_compound(rhs, "+=",
			[](F& left, const F& right) { left += right; });
	}

	// Subtracts an equally shaped tensor from this tensor in place.
	Tensor<F>& operator-=(const Tensor<F>& rhs) {
		return binary_compound(rhs, "-=",
			[](F& left, const F& right) { left -= right; });
	}

	// Multiplies this tensor by an equally shaped tensor in place.
	Tensor<F>& operator*=(const Tensor<F>& rhs) {
		return binary_compound(rhs, "*=",
			[](F& left, const F& right) { left *= right; });
	}

	// Divides this tensor by an equally shaped tensor in place.
	Tensor<F>& operator/=(const Tensor<F>& rhs) {
		return binary_compound(rhs, "/=",
			[](F& left, const F& right) { left /= right; });
	}

	// Returns a tensor with every element multiplied by a scalar.
	Tensor<F> operator*(F rhs) const {
		return scalar_transform(rhs,
			[](const F& left, const F& right) { return left * right; });
	}

	// Returns a tensor with a scalar added to every element.
	Tensor<F> operator+(F rhs) const {
		return scalar_transform(rhs,
			[](const F& left, const F& right) { return left + right; });
	}

	// Returns a tensor with every element divided by a scalar.
	Tensor<F> operator/(F rhs) const {
		return scalar_transform(rhs,
			[](const F& left, const F& right) { return left / right; });
	}

	// Returns a tensor with a scalar subtracted from every element.
	Tensor<F> operator-(F rhs) const {
		return scalar_transform(rhs,
			[](const F& left, const F& right) { return left - right; });
	}

	// Multiplies every element by a scalar in place.
	Tensor<F>& operator*=(F value) {
		return scalar_compound(value,
			[](F& left, const F& right) { left *= right; });
	}

	// Adds a scalar to every element in place.
	Tensor<F>& operator+=(F value) {
		return scalar_compound(value,
			[](F& left, const F& right) { left += right; });
	}

	// Divides every element by a scalar in place.
	Tensor<F>& operator/=(F value) {
		return scalar_compound(value,
			[](F& left, const F& right) { left /= right; });
	}

	// Subtracts a scalar from every element in place.
	Tensor<F>& operator-=(F value) {
		return scalar_compound(value,
			[](F& left, const F& right) { left -= right; });
	}

	// Assigns one scalar value to every tensor element.
	Tensor<F>& operator=(F value) {
		return scalar_compound(value,
			[](F& left, const F& right) { left = right; });
	}


	// Reports whether equally shaped tensors match within a tolerance.
	bool is_close_to(const Tensor<F>& rhs, F eps) const {
		if (this == &rhs)
			return true;
		if (!shape.same_dimensions(rhs.shape))
			throw std::invalid_argument("is_close_to : incompatible dimensions.");
		for (size_t i = 0; i < size; i++) {
			F dif = data[i] - rhs.data[i];
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
			out += to_str(data[i]) + " ";
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

	// Returns all tensor dimensions in a capacity-safe container.
	std::vector<size_t> dimensions() const {
		return shape.dimensions();
	}

	// Returns one checked tensor dimension.
	size_t dimension(unsigned i) const {
		if (i >= rank)
			throw std::out_of_range("tensor dimension is out of bounds.");
		return shape.dimension(i);
	}

	// Returns a shared mutable reference to the underlying payload.
	Payload<F> get_data() {
		return data;
	}

	// Returns an immutable pointer to the first visible element.
	const F* get_data() const {
		return data.get_arr();
	}

	// Applies any compatible callable to every tensor element in place.
	template<typename Function,
		typename std::enable_if<std::is_invocable_r<F, Function&, F>::value, int>::type = 0>
	Tensor<F>& apply(Function&& function) {
		F* values = data.get_arr();
		for (size_t i = 0; i < size; i++)
			values[i] = std::invoke(function, values[i]);
		return *this;
	}

	// Returns a mutable view selected by the first-axis index.
	View operator[](size_t i) {
		if (!rank || i >= shape.dimension(0))
			return View();
		if (rank == 1)
			return View(data, shape.scalar_shape(), 1, i);
		return View(data, shape.slice(1), size / shape.dimension(0),
			i * shape.stride(0));
	}

	// Returns a read-only shared view selected from a const tensor.
	ConstView operator[](size_t i) const {
		return ConstView(data, shape, size, 0)[i];
	}

	// Returns a mutable view selected by a leading index list.
	View operator[](const std::initializer_list<size_t>& index) {
		if (rank < index.size() || !rank)
			return View();
		for (unsigned r = 0; r < index.size(); r++)
			if (index.begin()[r] >= shape.dimension(r))
				return View();

		if (rank == index.size()) {
			size_t pos = 0;
			for (unsigned i = 0; i < rank; i++)
				pos += index.begin()[i] * shape.stride(i);
			return View(data, shape.scalar_shape(), 1, pos);
		}
		size_t out_off = 0;
		size_t out_size = size;
		for (unsigned i = 0; i < index.size(); i++) {
			out_off += index.begin()[i] * shape.stride(i);
			out_size /= shape.dimension(i);
		}
		return View(data, shape.slice(index.size()), out_size, out_off);
	}

	// Returns a read-only shared view selected from a const tensor.
	ConstView operator[](const std::initializer_list<size_t>& index) const {
		return ConstView(data, shape, size, 0)[index];
	}

private:

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

	// Contracts tensor-like operands through a shared packed-matrix implementation.
	template<typename LeftTensor, typename RightTensor>
	static Tensor<F> contract(const LeftTensor& lt, const RightTensor& rt,
		unsigned li, unsigned ri, const std::string& out_meta) {
		unsigned left_rank = lt.get_rank();
		unsigned right_rank = rt.get_rank();
		if (left_rank == 0 || right_rank == 0)
			return Tensor<F>(out_meta);
		if (li >= left_rank || ri >= right_rank
			|| lt.dimension(li) != rt.dimension(ri))
			throw std::invalid_argument("cnt : out of rank index or mismatch dimension.");
		unsigned out_rank = left_rank + right_rank - 2;
		const F* left_data = lt.get_data();
		const F* right_data = rt.get_data();
		if (out_rank == 0) {
			F val{};
			for (size_t k = 0; k < lt.dimension(0); k++)
				val += left_data[k] * right_data[k];
			Tensor<F> out({ 1 }, val);
			out.set_meta(out_meta);
			return out;
		}

		std::vector<size_t> out_dim;
		out_dim.reserve(out_rank);
		for (unsigned i = 0; i < left_rank; i++)
			if (i != li)
				out_dim.push_back(lt.dimension(i));
		for (unsigned i = 0; i < right_rank; i++)
			if (i != ri)
				out_dim.push_back(rt.dimension(i));
		Tensor<F> out(TensorShape(out_dim.data(), out_rank), UninitializedTag{});
		out.set_meta(out_meta);

		size_t k_size = lt.dimension(li);
		size_t m_size = lt.get_size() / k_size;
		size_t n_size = rt.get_size() / k_size;

		Payload<F> l_matrix;
		const F* l_arr = left_data;
		if (li != left_rank - 1) {
			l_matrix = Payload<F>(lt.get_size());
			F* packed_left = l_matrix.get_arr();
			size_t suffix_size = 1;
			for (unsigned i = li + 1; i < left_rank; i++)
				suffix_size *= lt.dimension(i);
			for (size_t pos = 0; pos < lt.get_size(); pos++) {
				size_t k = (pos / suffix_size) % k_size;
				size_t row = (pos / (k_size * suffix_size)) * suffix_size
					+ pos % suffix_size;
				packed_left[row * k_size + k] = left_data[pos];
			}
			l_arr = l_matrix.get_arr();
		}

		Payload<F> r_matrix;
		const F* r_arr = right_data;
		if (ri != 0) {
			r_matrix = Payload<F>(rt.get_size());
			F* packed_right = r_matrix.get_arr();
			size_t suffix_size = 1;
			for (unsigned i = ri + 1; i < right_rank; i++)
				suffix_size *= rt.dimension(i);
			for (size_t pos = 0; pos < rt.get_size(); pos++) {
				size_t k = (pos / suffix_size) % k_size;
				size_t col = (pos / (k_size * suffix_size)) * suffix_size
					+ pos % suffix_size;
				packed_right[k * n_size + col] = right_data[pos];
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

public:
	// Contracts one axis from each owned tensor through packed matrix multiplication.
	static Tensor<F> cnt(const Tensor<F>& lt, const Tensor<F>& rt, unsigned li, unsigned ri) {
		return contract(lt, rt, li, ri,
			"cnt(" + lt.get_meta() + ":" + rt.get_meta() + ");");
	}

	// Contracts tensor owners or views without materializing view payloads.
	template<typename LeftTensor, typename RightTensor>
	static Tensor<F> cnt(const LeftTensor& lt, const RightTensor& rt, unsigned li, unsigned ri) {
		return contract(lt, rt, li, ri, "cnt(view);");
	}

	// Contracts this tensor's last axis with another tensor's first axis.
	Tensor<F> operator%(const Tensor<F>& rhs) const {
		return cnt(*this, rhs, rank - 1, 0);
	}

	// Replaces tensor metadata from a C string.
	void set_meta(const char* meta_str) {
		meta = meta_str ? meta_str : "";
	}

	// Replaces tensor metadata from a string.
	void set_meta(const std::string& meta_str) {
		set_meta(meta_str.c_str());
	}

	// Returns the tensor metadata text without copying it.
	const std::string& get_meta(void) const {
		return meta;
	}

	// Appends text to the tensor metadata.
	void append_meta(const std::string& mt_str) {
		meta += mt_str;
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
		F* values = out.data.get_arr();
		for (size_t i = 0; i < out.size; i++)
			values[i] = start + i * jump;
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
