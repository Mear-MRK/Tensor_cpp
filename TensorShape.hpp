#pragma once
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

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
