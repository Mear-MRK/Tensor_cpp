#pragma once
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include "TensorDebug.hpp"

// Manages a reference-counted contiguous array of values.
template <typename F>
class Payload
{
private:
	size_t size = 0;
	std::shared_ptr<F[]> arr;

	// Allocates an array managed by a standard shared control block.
	static std::shared_ptr<F[]> allocate(size_t t_size)
	{
		return t_size ? std::shared_ptr<F[]>(new F[t_size]) : nullptr;
	}

public:
	// Creates an empty payload.
	Payload() { TENSOR_LOG("pyl cns, dfl " << this); }

	// Allocates an uninitialized payload with the requested size.
	explicit Payload(size_t t_size) : size(t_size), arr(allocate(t_size))
	{
		TENSOR_LOG("pyl cns, siz " << this);
	}

	// Copies an array into a newly allocated payload.
	Payload(const F array[], size_t t_size) : size(t_size), arr(allocate(t_size))
	{
		if (size && !array)
			throw std::invalid_argument("payload input array cannot be null.");
		if (size)
			std::copy_n(array, size, arr.get());
		TENSOR_LOG("pyl cns, arr " << this);
	}

	// Shares the referenced storage with another payload.
	Payload(const Payload<F> &oth) : size(oth.size), arr(oth.arr)
	{
		TENSOR_LOG("pyl cns, cop " << &oth << "-c->" << this);
	}

	// Replaces this payload with shared storage from another payload.
	Payload<F> &operator=(const Payload<F> &rhs)
	{
		if (this == &rhs)
			return *this;
		size = rhs.size;
		arr = rhs.arr;
		TENSOR_LOG("pyl asn, cop " << &rhs << "-c=>" << this);
		return *this;
	}

	// Transfers storage ownership from another payload.
	Payload(Payload<F> &&oth) noexcept : size(oth.size), arr(std::move(oth.arr))
	{
		oth.size = 0;
		TENSOR_LOG("pyl cns, mov " << &oth << "-m->" << this);
	}

	// Replaces this payload by moving storage from another payload.
	Payload<F> &operator=(Payload<F> &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;
		size = rhs.size;
		arr = std::move(rhs.arr);
		rhs.size = 0;
		TENSOR_LOG("pyl asn, mov " << &rhs << "-m=>" << this);
		return *this;
	}

	// Allocates a payload and fills every element with a value.
	Payload(size_t t_size, F value) : size(t_size), arr(allocate(t_size))
	{
		if (size)
			std::fill_n(arr.get(), size, value);
		TENSOR_LOG("pyl cns, szv " << this);
	}

	// Allocates a payload filled with evenly spaced values.
	Payload(size_t t_size, F start, F end) : size(t_size), arr(allocate(t_size))
	{
		if (size > 1)
		{
			F step = (end - start) / (size - 1);
			for (size_t i = 0; i < size; i++)
				arr[i] = start + i * step;
		}
		else if (size)
		{
			arr[0] = start;
		}
		TENSOR_LOG("pyl cns, lin " << this);
	}

	// Releases this payload's reference to its storage.
	void release()
	{
		arr.reset();
		size = 0;
	}

	// Returns an independent copy of a selected data range.
	Payload<F> copy(size_t offset = 0, size_t new_size = ~0ull) const
	{
		if (offset >= get_size())
			return Payload<F>();
		size_t m_size = std::min<size_t>(new_size, get_size() - offset);
		Payload<F> out = Payload<F>(m_size);
		if (m_size)
			std::copy_n(get_arr() + offset, m_size, out.get_arr());
		return out;
	}

	// Copies a bounded range from another payload into this payload.
	Payload<F> &assign(const Payload<F> &oth, size_t oth_off = 0, size_t cp_size = ~0ull, size_t cp_off = 0, size_t oth_size = ~0ull)
	{
		if (oth_off >= oth.get_size() || !cp_size || !oth_size || !oth.get_arr())
			return *this;
		if (oth_size == ~0ull || oth_size > oth.get_size() - oth_off)
			oth_size = oth.get_size() - oth_off;
		if (cp_size > oth_size)
			cp_size = oth_size;
		if (!get_size())
		{
			if (cp_size > std::numeric_limits<size_t>::max() - cp_off)
				throw std::overflow_error("payload assignment size overflow.");
			init(cp_off + cp_size);
		}
		if (cp_off >= get_size())
			return *this;
		if (cp_size > get_size() - cp_off)
			cp_size = get_size() - cp_off;

		F *destination = get_arr() + cp_off;
		const F *source = oth.get_arr() + oth_off;
		if (destination == source)
			return *this;
		if constexpr (std::is_trivially_copyable<F>::value)
		{
			std::memmove(destination, source, cp_size * sizeof(F));
		}
		else if (arr == oth.arr && destination > source && destination < source + cp_size)
		{
			std::copy_backward(source, source + cp_size, destination + cp_size);
		}
		else
		{
			std::copy_n(source, cp_size, destination);
		}
		return *this;
	}

	// Initializes empty storage with the requested size.
	bool init(size_t t_size)
	{
		if (arr)
		{
			TENSOR_LOG("pyl init, already initialized " << this);
			return false;
		}
		if (t_size)
		{
			std::shared_ptr<F[]> new_arr = allocate(t_size);
			size = t_size;
			arr = std::move(new_arr);
			return true;
		}
		return false;
	}

	// Initializes empty storage by copying an input array.
	bool init(const F array[], size_t t_size)
	{
		if (arr)
		{
			TENSOR_LOG("pyl init, already initialized " << this);
			return false;
		}
		if (t_size)
		{
			if (!array)
				throw std::invalid_argument("payload input array cannot be null.");
			std::shared_ptr<F[]> new_arr = allocate(t_size);
			std::copy_n(array, t_size, new_arr.get());
			size = t_size;
			arr = std::move(new_arr);
			return true;
		}
		return false;
	}

	// Releases the payload's standard shared array reference.
	~Payload()
	{
		TENSOR_LOG("pyl dst,     " << this);
	}

	// Returns the current storage reference count or -1 when empty.
	int get_ref_count() const
	{
		if (arr)
			return static_cast<int>(arr.use_count());
		return -1;
	}

	// Returns checked mutable access to an element.
	inline F &operator[](size_t i)
	{
		if (!arr || i >= size)
			throw std::out_of_range("payload index is out of bounds.");
		return arr[i];
	}

	// Returns checked immutable access to an element.
	inline const F &operator[](size_t i) const
	{
		if (!arr || i >= size)
			throw std::out_of_range("payload index is out of bounds.");
		return arr[i];
	}


	// Returns the mutable raw storage pointer.
	inline F *get_arr()
	{
		return arr.get();
	}

	// Returns the immutable raw storage pointer.
	inline const F *get_arr() const
	{
		return arr.get();
	}

	// Returns the number of stored elements.
	inline size_t get_size() const
	{
		return size;
	}
};
