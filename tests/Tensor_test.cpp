#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../Tensor.hpp"

size_t passed_tensor_tests = 0;

// Throws with a diagnostic message when a tensor test fails.
void require_tensor(bool condition, const std::string& message) {
	if (!condition)
		throw std::runtime_error(message);
}

// Verifies that a tensor contains the expected contiguous values.
template<typename F>
void require_values(const Tensor<F>& tensor, const std::vector<F>& expected,
	const std::string& message) {
	require_tensor(tensor.get_size() == expected.size(), message + ": size mismatch");
	const F* values = tensor.get_data();
	for (size_t i = 0; i < expected.size(); i++)
		require_tensor(values[i] == expected[i], message + ": value mismatch");
}

// Creates a mutable view whose original tensor has already been destroyed.
Tensor<int>::View make_persistent_view() {
	Tensor<int> source = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	return source[1];
}

// Verifies constructors, generated values, ranges, and reshaping.
void test_construction() {
	Tensor<double> filled({ 2, 3, 4 }, 1.0);
	require_tensor(filled.dimensions() == std::vector<size_t>({ 2, 3, 4 }),
		"constant tensor dimensions mismatch");
	require_tensor(filled.get_size() == 24, "constant tensor size mismatch");
	for (size_t i = 0; i < filled.get_size(); i++)
		require_tensor(filled.get_data()[i] == 1.0, "constant tensor value mismatch");

	int next = 0;
	Tensor<int> generated({ 2, 3 }, [&next]() { return ++next; });
	require_values(generated, { 1, 2, 3, 4, 5, 6 }, "capturing generator");

	Tensor<int> range = Tensor<int>::range(-2, 7, 2);
	require_values(range, { -2, 0, 2, 4, 6 }, "range");
	range.reshape({ 5, 1 });
	require_tensor(range.dimensions() == std::vector<size_t>({ 5, 1 }),
		"reshape dimensions mismatch");

	Tensor<double> line = Tensor<double>::linspace(1.0, 2.0, 3);
	require_values(line, { 1.0, 1.5, 2.0 }, "linspace");
	passed_tensor_tests += 7;
}

// Verifies binary, scalar, compound, and callable arithmetic operations.
void test_arithmetic() {
	Tensor<int> left = Tensor<int>::range(1, 5).reshape({ 2, 2 });
	Tensor<int> right({ 2, 2 }, 2);
	require_values(left + right, { 3, 4, 5, 6 }, "tensor addition");
	require_values(left - right, { -1, 0, 1, 2 }, "tensor subtraction");
	require_values(left * right, { 2, 4, 6, 8 }, "tensor multiplication");
	require_values((left * right) / right, { 1, 2, 3, 4 }, "tensor division");
	require_values(left + 3, { 4, 5, 6, 7 }, "scalar addition");
	require_values(left - 1, { 0, 1, 2, 3 }, "scalar subtraction");
	require_values(left * 3, { 3, 6, 9, 12 }, "scalar multiplication");
	require_values((left * 4) / 2, { 2, 4, 6, 8 }, "scalar division");

	const Tensor<int> constant({ 2, 2 }, 1);
	Tensor<int> compound({ 2, 2 }, 2);
	compound += constant;
	compound -= Tensor<int>({ 2, 2 }, 1);
	compound *= constant + 1;
	compound /= constant + 1;
	require_values(compound, { 2, 2, 2, 2 }, "compound tensor operands");
	compound += 2;
	compound -= 1;
	compound *= 3;
	compound /= 3;
	require_values(compound, { 3, 3, 3, 3 }, "compound scalar operands");

	int captured = 4;
	compound.apply([captured](int value) { return value + captured; });
	require_values(compound, { 7, 7, 7, 7 }, "capturing apply callable");
	passed_tensor_tests += 11;
}

// Verifies explicit mutable and immutable tensor-view behavior.
void test_views() {
	Tensor<int> tensor = Tensor<int>::range(1, 13).reshape({ 2, 2, 3 });
	auto plane = tensor[1];
	static_assert(std::is_same<decltype(plane), Tensor<int>::View>::value,
		"mutable slicing must return Tensor::View");
	require_tensor(plane.dimensions() == std::vector<size_t>({ 2, 3 }),
		"mutable view dimensions mismatch");

	auto row = tensor[{ 1, 0 }];
	row = Tensor<int>::range(20, 23);
	require_tensor(tensor.get_data()[6] == 20 && tensor.get_data()[8] == 22,
		"view assignment did not update owner");
	auto destination = tensor[{ 0, 1 }];
	destination = row;
	require_tensor(tensor.get_data()[3] == 20 && tensor.get_data()[5] == 22,
		"view-to-view assignment failed");

	Tensor<int> owned_copy = row;
	row = 9;
	require_values(owned_copy, { 20, 21, 22 }, "owned copy from view");
	require_tensor(tensor.get_data()[6] == 9 && tensor.get_data()[8] == 9,
		"scalar view assignment failed");

	const Tensor<int>& const_tensor = tensor;
	auto const_view = const_tensor[0];
	static_assert(std::is_same<decltype(const_view), Tensor<int>::ConstView>::value,
		"const slicing must return Tensor::ConstView");
	static_assert(std::is_same<decltype(const_view.get_data()), const int*>::value,
		"const view data must be immutable");
	require_tensor(const_view.at({ 1, 2 }) == 22, "const view lookup failed");

	auto persistent = make_persistent_view();
	require_tensor(persistent.at({ 0 }) == 4 && persistent.at({ 2 }) == 6,
		"mutable view did not retain storage lifetime");
	persistent.apply([](int value) { return value * 2; });
	require_tensor(persistent.at({ 0 }) == 8 && persistent.at({ 2 }) == 12,
		"callable application to view failed");
	passed_tensor_tests += 9;
}

// Verifies checked dimensions, full metadata, copies, and exception paths.
void test_api_and_errors() {
	Tensor<int> tensor = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	require_tensor(tensor.dimension(0) == 2 && tensor.dimension(1) == 3,
		"checked dimension access failed");
	bool threw = false;
	try {
		tensor.dimension(2);
	}
	catch (const std::out_of_range&) {
		threw = true;
	}
	require_tensor(threw, "out-of-range dimension access did not throw");

	std::string long_meta(160, 'm');
	tensor.set_meta(long_meta);
	tensor.append_meta("-tail");
	require_tensor(tensor.get_meta() == long_meta + "-tail", "metadata was truncated");

	Tensor<int> independent = tensor.copy();
	tensor[0] = 0;
	require_tensor(independent.get_data()[0] == 1,
		"explicit tensor copy shared mutable storage");

	threw = false;
	try {
		Tensor<int> mismatch({ 3, 2 }, 1);
		static_cast<void>(tensor + mismatch);
	}
	catch (const std::invalid_argument&) {
		threw = true;
	}
	require_tensor(threw, "mismatched arithmetic did not throw");

	threw = false;
	try {
		tensor.reshape({ 7 });
	}
	catch (const std::invalid_argument&) {
		threw = true;
	}
	require_tensor(threw, "mismatched reshape did not throw");
	require_tensor(!Tensor<int>().to_string().empty(), "empty tensor formatting failed");

	Payload<int> empty_payload;
	threw = false;
	try {
		static_cast<void>(empty_payload[0]);
	}
	catch (const std::out_of_range&) {
		threw = true;
	}
	require_tensor(threw, "empty payload access did not throw");
	passed_tensor_tests += 8;
}

// Runs the automated general tensor test suite.
int main() {
	try {
		test_construction();
		test_arithmetic();
		test_views();
		test_api_and_errors();
		std::cout << "Tensor tests passed: " << passed_tensor_tests << std::endl;
	}
	catch (const std::exception& error) {
		std::cerr << "Tensor test failed: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}
