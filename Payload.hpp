#pragma once
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <stdexcept>

#ifdef DEBUG
#define LOG(s) std::cout << s << std::endl
#else
#define LOG(s)
#endif

// Manages a reference-counted contiguous array of values.
template <typename F>
class Payload
{
private:
	int *ref_count;
	size_t size;
	F *arr;

public:
	// Creates an empty payload.
	Payload() : ref_count(nullptr), size(0), arr(nullptr) { LOG("pyl cns, dfl " << this); }

	// Allocates an uninitialized payload with the requested size.
	explicit Payload(size_t size) : size(size)
	{
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
		}
		LOG("pyl cns, siz " << this);
	}

	// Copies an array into a newly allocated payload.
	Payload(const F array[], size_t size) : size(size)
	{
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
			std::copy_n(array, size, arr);
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
		}
		LOG("pyl cns, arr " << this);
	}

	// Shares the referenced storage with another payload.
	Payload(const Payload<F> &oth) : ref_count(oth.ref_count), size(oth.size), arr(oth.arr)
	{
		inc_ref_count();
		LOG("pyl cns, cop " << &oth << "-c->" << this);
	}

	// Replaces this payload with shared storage from another payload.
	Payload<F> &operator=(const Payload<F> &rhs)
	{
		if (this == &rhs)
			return *this;
		dec_ref_count();
		ref_count = rhs.ref_count;
		size = rhs.size;
		arr = rhs.arr;
		inc_ref_count();
		LOG("pyl asn, cop " << &rhs << "-c=>" << this);
		return *this;
	}

	// Transfers storage ownership from another payload.
	Payload(Payload<F> &&oth) noexcept : ref_count(oth.ref_count), size(oth.size), arr(oth.arr)
	{
		oth.ref_count = nullptr;
		oth.size = 0;
		oth.arr = nullptr;
		LOG("pyl cns, mov " << &oth << "-m->" << this);
	}

	// Replaces this payload by moving storage from another payload.
	Payload<F> &operator=(Payload<F> &&rhs) noexcept
	{
		if (this == &rhs)
			return *this;
		dec_ref_count();
		ref_count = rhs.ref_count;
		size = rhs.size;
		arr = rhs.arr;
		rhs.ref_count = nullptr;
		rhs.size = 0;
		rhs.arr = nullptr;
		LOG("pyl asn, mov " << &rhs << "-m=>" << this);
		return *this;
	}

	// Allocates a payload and fills every element with a value.
	Payload(size_t size, F value) : size(size)
	{
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
			std::fill_n(arr, size, value);
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
		}
		LOG("pyl cns, szv " << this);
	}

	// Allocates a payload filled with evenly spaced values.
	Payload(size_t size, F start, F end) : size(size)
	{
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
			if (size > 1)
			{
				F step = (end - start) / (size - 1);
				for (size_t i = 0; i < size; i++)
				{
					arr[i] = start + i * step;
				}
			}
			else
			{
				arr[0] = start;
			}
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
		}
		LOG("pyl cns, lin " << this);
	}

	// Releases this payload's reference to its storage.
	void release()
	{
		dec_ref_count();
		size = 0;
		ref_count = nullptr;
		arr = nullptr;
	}

	// Returns an independent copy of a selected data range.
	Payload<F> copy(size_t offset = 0, size_t new_size = ~0ull) const
	{
		if (offset >= size)
			return Payload<F>();
		size_t m_size = std::min<size_t>(new_size, size);
		m_size = (m_size + offset <= size) ? m_size : size - offset;
		Payload<F> out = Payload<F>(m_size);
		if (m_size)
		{
			memcpy(out.arr, arr + offset, m_size * sizeof(F));
		}
		return out;
	}

	// Copies a bounded range from another payload into this payload.
	Payload<F> &assign(const Payload<F> &oth, size_t oth_off = 0, size_t cp_size = ~0ull, size_t cp_off = 0, size_t oth_size = ~0ull)
	{
		if (this == &oth || oth_off >= oth.size || !cp_size || !oth_size || !oth.arr)
			return *this;
		if (oth_size == ~0ull || oth_off + oth_size > oth.size)
			oth_size = oth.size - oth_off;
		if (cp_size > oth_size)
			cp_size = oth_size;
		if (!size)
		{
			init(cp_off + cp_size);
		}
		if (cp_off >= size)
			return *this;
		if (cp_size + cp_off > size)
			cp_size = size - cp_off;

		memcpy(arr + cp_off, oth.arr + oth_off, cp_size * sizeof(F));
		return *this;
	}

	// Initializes empty storage with the requested size.
	bool init(size_t size)
	{
		if (ref_count)
		{
			LOG("pyl init, already initialized " << this);
			return false;
		}
		this->size = size;
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
			return true;
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
			return false;
		}
	}

	// Initializes empty storage by copying an input array.
	bool init(const F array[], size_t size)
	{
		if (ref_count)
		{
			LOG("pyl init, already initialized " << this);
			return false;
		}
		this->size = size;
		if (size)
		{
			ref_count = new int(1);
			arr = new F[size];
			std::copy_n(array, size, arr);
			return true;
		}
		else
		{
			ref_count = nullptr;
			arr = nullptr;
			return false;
		}
	}

	// Releases the payload's reference and destroys unshared storage.
	~Payload()
	{
		dec_ref_count();
		LOG("pyl dst,     " << this);
	}

	// Returns the current storage reference count or -1 when empty.
	int get_ref_count() const
	{
		if (ref_count != nullptr)
			return *ref_count;
		return -1;
	}

	// Returns mutable unchecked access to an element.
	inline F &operator[](size_t i)
	{
		// Doesn't check if arr is null or index out of bound
		return arr[i];
	}

	// Returns immutable unchecked access to an element.
	inline const F &operator[](size_t i) const
	{
		// Doesn't check if arr is null or index out of bound
		return arr[i];
	}


	// Returns the mutable raw storage pointer.
	inline F *get_arr()
	{
		return arr;
	}

	// Returns the immutable raw storage pointer.
	inline const F *get_arr() const
	{
		return arr;
	}

	// Returns the number of stored elements.
	inline size_t get_size() const
	{
		return size;
	}

private:
	// Increments the shared storage reference count when present.
	void inline inc_ref_count()
	{
		if (ref_count != nullptr)
			(*ref_count)++;
	}

	// Decrements the reference count and frees unreferenced storage.
	void inline dec_ref_count()
	{
		if (ref_count != nullptr)
		{
			(*ref_count)--;
			if (*ref_count == 0)
			{
				delete[] arr;
				delete ref_count;
			}
		}
	}
};
