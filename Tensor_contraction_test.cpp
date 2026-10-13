#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Tensor.hpp"

size_t passed_tests = 0;

// Throws with a diagnostic message when a test condition fails.
void require(bool condition, const std::string& message) {
	if (!condition)
		throw std::runtime_error(message);
}

// Returns the element count represented by a tensor shape.
size_t tensor_size(const std::vector<size_t>& dims) {
	size_t size = 1;
	for (size_t dim : dims)
		size *= dim;
	return size;
}

// Decodes a row-major position into a multidimensional index.
std::vector<size_t> decode_index(size_t pos, const std::vector<size_t>& dims) {
	std::vector<size_t> index(dims.size());
	for (size_t i = dims.size(); i > 0; i--) {
		index[i - 1] = pos % dims[i - 1];
		pos /= dims[i - 1];
	}
	return index;
}

// Encodes a multidimensional index as a row-major position.
size_t encode_index(const std::vector<size_t>& index, const std::vector<size_t>& dims) {
	size_t pos = 0;
	for (size_t i = 0; i < dims.size(); i++)
		pos = pos * dims[i] + index[i];
	return pos;
}

// Returns exact-comparison tolerance for integral element types.
template<typename F>
double tolerance() {
	return 0.0;
}

// Returns comparison tolerance for float elements.
template<>
double tolerance<float>() {
	return 1e-5;
}

// Returns comparison tolerance for double elements.
template<>
double tolerance<double>() {
	return 1e-12;
}

// Returns comparison tolerance for complex-float elements.
template<>
double tolerance<std::complex<float>>() {
	return 1e-5;
}

// Returns comparison tolerance for complex-double elements.
template<>
double tolerance<std::complex<double>>() {
	return 1e-12;
}

// Compares two values using type-appropriate relative tolerance.
template<typename F>
bool is_close(F left, F right) {
	double dif = static_cast<double>(std::abs(left - right));
	double scale = 1.0 + std::max(static_cast<double>(std::abs(left)),
		static_cast<double>(std::abs(right)));
	return dif <= tolerance<F>() * scale;
}

// Produces a deterministic scalar test value for an index and seed.
template<typename F>
F test_value(size_t i, int seed) {
	int value = static_cast<int>((i * 7 + seed * 3) % 17) - 8;
	return static_cast<F>(value) / static_cast<F>(4);
}

// Produces a deterministic integer test value for an index and seed.
template<>
int test_value<int>(size_t i, int seed) {
	return static_cast<int>((i * 7 + seed * 3) % 17) - 8;
}

// Produces a deterministic complex-float test value.
template<>
std::complex<float> test_value<std::complex<float>>(size_t i, int seed) {
	float real = static_cast<float>(static_cast<int>((i * 5 + seed) % 11) - 5) / 3.0f;
	float imag = static_cast<float>(static_cast<int>((i * 3 + seed) % 7) - 3) / 4.0f;
	return { real, imag };
}

// Produces a deterministic complex-double test value.
template<>
std::complex<double> test_value<std::complex<double>>(size_t i, int seed) {
	double real = static_cast<double>(static_cast<int>((i * 5 + seed) % 11) - 5) / 3.0;
	double imag = static_cast<double>(static_cast<int>((i * 3 + seed) % 7) - 3) / 4.0;
	return { real, imag };
}

// Builds deterministic test data of the requested size.
template<typename F>
std::vector<F> make_values(size_t size, int seed) {
	std::vector<F> values(size);
	for (size_t i = 0; i < size; i++)
		values[i] = test_value<F>(i, seed);
	return values;
}

// Compares optimized contraction output with a direct reference calculation.
template<typename F, typename LeftTensor, typename RightTensor>
void verify_contraction(const LeftTensor& left, const std::vector<size_t>& left_dims,
	const std::vector<F>& left_values, unsigned li, const RightTensor& right,
	const std::vector<size_t>& right_dims, const std::vector<F>& right_values,
	unsigned ri, const std::string& name) {
	Tensor<F> result = Tensor<F>::cnt(left, right, li, ri);
	std::vector<size_t> out_dims;
	for (size_t i = 0; i < left_dims.size(); i++)
		if (i != li)
			out_dims.push_back(left_dims[i]);
	for (size_t i = 0; i < right_dims.size(); i++)
		if (i != ri)
			out_dims.push_back(right_dims[i]);
	if (out_dims.empty())
		out_dims.push_back(1);

	require(result.get_size() == tensor_size(out_dims), name + ": output size mismatch");
	require(result.dimensions() == out_dims, name + ": output dimensions mismatch");

	Payload<F> result_data = result.get_data();
	for (size_t pos = 0; pos < result.get_size(); pos++) {
		std::vector<size_t> out_index = decode_index(pos, out_dims);
		std::vector<size_t> left_index(left_dims.size());
		std::vector<size_t> right_index(right_dims.size());
		size_t out_i = 0;
		for (size_t i = 0; i < left_dims.size(); i++)
			if (i != li)
				left_index[i] = out_index[out_i++];
		for (size_t i = 0; i < right_dims.size(); i++)
			if (i != ri)
				right_index[i] = out_index[out_i++];

		F expected = F{};
		for (size_t k = 0; k < left_dims[li]; k++) {
			left_index[li] = k;
			right_index[ri] = k;
			expected += left_values[encode_index(left_index, left_dims)]
				* right_values[encode_index(right_index, right_dims)];
		}
		require(is_close(result_data[pos], expected), name + ": output value mismatch");
	}
	passed_tests++;
}

// Tests every compatible contraction-axis pair for two tensor ranks.
template<typename F>
void test_axis_pairs(const std::vector<size_t>& left_dims,
	const std::vector<size_t>& right_base_dims, const std::string& name) {
	for (unsigned li = 0; li < left_dims.size(); li++) {
		for (unsigned ri = 0; ri < right_base_dims.size(); ri++) {
			std::vector<size_t> right_dims = right_base_dims;
			right_dims[ri] = left_dims[li];
			std::vector<F> left_values = make_values<F>(tensor_size(left_dims), li + 1);
			std::vector<F> right_values = make_values<F>(tensor_size(right_dims), ri + 7);
			Tensor<F> left(left_dims.data(), static_cast<unsigned>(left_dims.size()), left_values.data());
			Tensor<F> right(right_dims.data(), static_cast<unsigned>(right_dims.size()), right_values.data());
			verify_contraction(left, left_dims, left_values, li,
				right, right_dims, right_values, ri, name);
		}
	}
}

// Verifies contractions performed on shared tensor views.
template<typename F>
void test_views(const std::string& name) {
	std::vector<size_t> left_base_dims{ 2, 2, 3, 4 };
	std::vector<size_t> right_base_dims{ 2, 5, 3, 6 };
	std::vector<F> left_base_values = make_values<F>(tensor_size(left_base_dims), 11);
	std::vector<F> right_base_values = make_values<F>(tensor_size(right_base_dims), 13);
	Tensor<F> left_base(left_base_dims.data(), static_cast<unsigned>(left_base_dims.size()), left_base_values.data());
	Tensor<F> right_base(right_base_dims.data(), static_cast<unsigned>(right_base_dims.size()), right_base_values.data());
	auto left_view = left_base[1];
	auto right_view = right_base[1];
	std::vector<F> left_values(left_base_values.begin() + 24, left_base_values.end());
	std::vector<F> right_values(right_base_values.begin() + 90, right_base_values.end());
	verify_contraction(left_view, { 2, 3, 4 }, left_values, 1,
		right_view, { 5, 3, 6 }, right_values, 1, name);
}

// Verifies rank-one contraction as a dot product.
template<typename F>
void test_dot_product(const std::string& name) {
	std::vector<size_t> dims{ 7 };
	std::vector<F> left_values = make_values<F>(7, 17);
	std::vector<F> right_values = make_values<F>(7, 19);
	Tensor<F> left(dims.data(), 1, left_values.data());
	Tensor<F> right(dims.data(), 1, right_values.data());
	verify_contraction(left, dims, left_values, 0,
		right, dims, right_values, 0, name);
}

// Verifies dimension, axis, and empty-input validation behavior.
void test_invalid_inputs() {
	Tensor<int> left({ 2, 3 }, 1);
	Tensor<int> right({ 4, 2 }, 1);
	bool threw = false;
	try {
		Tensor<int>::cnt(left, right, 1, 0);
	}
	catch (const std::invalid_argument&) {
		threw = true;
	}
	require(threw, "mismatched contraction dimensions must throw");

	threw = false;
	try {
		Tensor<int>::cnt(left, right, 2, 0);
	}
	catch (const std::invalid_argument&) {
		threw = true;
	}
	require(threw, "out-of-range contraction index must throw");

	Tensor<int> empty;
	require(Tensor<int>::cnt(empty, right, 0, 0).get_size() == 0,
		"empty contraction input must return an empty tensor");
	passed_tests += 3;
}

// Runs the complete generic or BLAS-backed contraction test suite.
int main() {
	try {
		test_axis_pairs<int>({ 2, 3, 4 }, { 5, 6, 7 }, "rank-3 int");
		test_axis_pairs<float>({ 2, 3, 4 }, { 5, 6, 7 }, "rank-3 float");
		test_axis_pairs<double>({ 2, 3, 4 }, { 5, 6, 7 }, "rank-3 double");
		test_axis_pairs<std::complex<float>>({ 2, 3, 4 }, { 5, 6, 7 }, "rank-3 complex float");
		test_axis_pairs<std::complex<double>>({ 2, 3, 4 }, { 5, 6, 7 }, "rank-3 complex double");

		test_axis_pairs<int>({ 2, 3, 2, 4 }, { 3, 2, 5, 2 }, "rank-4 int");
		test_axis_pairs<float>({ 2, 3, 2, 4 }, { 3, 2, 5, 2 }, "rank-4 float");

		test_views<int>("int views");
		test_views<float>("float views");
		test_views<std::complex<double>>("complex views");

		test_dot_product<int>("int dot product");
		test_dot_product<float>("float dot product");
		test_dot_product<double>("double dot product");
		test_dot_product<std::complex<float>>("complex float dot product");
		test_dot_product<std::complex<double>>("complex double dot product");

		test_invalid_inputs();

		std::cout << "Tensor contraction tests passed: " << passed_tests;
#ifdef TENSOR_USE_BLAS
		std::cout << " (BLAS enabled)";
#else
		std::cout << " (generic matrix product)";
#endif
		std::cout << std::endl;
	}
	catch (const std::exception& error) {
		std::cerr << "Tensor contraction test failed: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}
