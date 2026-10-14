#pragma once
#include <complex>
#include <cstddef>
#include <limits>

#ifdef TENSOR_USE_BLAS
#if __has_include(<cblas.h>)
#include <cblas.h>
#elif __has_include(<Accelerate/Accelerate.h>)
#include <Accelerate/Accelerate.h>
#else
#error TENSOR_USE_BLAS requires a CBLAS-compatible header.
#endif
#endif

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
