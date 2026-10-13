#include <atomic>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

#include "Tensor.hpp"

size_t passed_refactor_tests = 0;

// Throws with a diagnostic message when a refactoring test fails.
void require_refactor(bool condition, const std::string& message) {
	if (!condition)
		throw std::runtime_error(message);
}

// Verifies that owned copy and move assignment replace complete tensor state.
void test_owned_assignment() {
	Tensor<int> source = Tensor<int>::range(1, 5);
	Tensor<int> target({ 2 }, 9);
	target = source;
	require_refactor(target.get_size() == 4, "copy assignment did not replace size");
	size_t dim[1]{};
	target.get_dim(dim);
	require_refactor(dim[0] == 4, "copy assignment did not replace shape");
	source[0] = 99;
	require_refactor(target.get_data()[0] == 1, "owned copy assignment shared payload");

	Tensor<int> move_source = Tensor<int>::range(3, 8);
	int* original_data = move_source.get_data().get_arr();
	Tensor<int> move_target({ 2 }, 0);
	move_target = std::move(move_source);
	require_refactor(move_target.get_data().get_arr() == original_data,
		"owned move assignment did not transfer payload");
	passed_refactor_tests += 4;
}

// Verifies matching view assignment and overlapping shared-storage assignment.
void test_view_assignment() {
	Tensor<int> base = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	Tensor<int> values({ 3 }, 8);
	Tensor<int> row = base[1];
	row = values;
	Payload<int> data = base.get_data();
	require_refactor(data[3] == 8 && data[4] == 8 && data[5] == 8,
		"view assignment did not update its source tensor");

	base = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	base[1] = base[0];
	data = base.get_data();
	require_refactor(data[3] == 1 && data[4] == 2 && data[5] == 3,
		"overlapping shared-storage assignment failed");
	base[0] = base[0];

	bool threw = false;
	try {
		Tensor<int> mismatch({ 2 }, 1);
		row = mismatch;
	}
	catch (const std::invalid_argument&) {
		threw = true;
	}
	require_refactor(threw, "mismatched view assignment did not throw");
	passed_refactor_tests += 3;
}

// Returns a const view whose source tensor has already been destroyed.
Tensor<int>::ConstView make_persistent_const_view() {
	Tensor<int> source = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	const Tensor<int>& const_source = source;
	return const_source[1];
}

// Verifies compile-time read-only access and shared const-view lifetime.
void test_const_views() {
	Tensor<int> source = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	const Tensor<int>& const_source = source;
	auto view = const_source[1];
	static_assert(std::is_same<decltype(view.get_data()), const int*>::value,
		"const view data must be immutable");
	static_assert(!std::is_assignable<decltype(view)&, int>::value,
		"const views must reject scalar assignment");
	require_refactor(view.get_rank() == 1 && view.get_size() == 3,
		"const view shape mismatch");
	require_refactor(view.at({ 0 }) == 4 && view.at({ 2 }) == 6,
		"const view values mismatch");
	auto persistent = make_persistent_const_view();
	require_refactor(persistent.at({ 0 }) == 4 && persistent.at({ 2 }) == 6,
		"const view did not retain shared storage");
	passed_refactor_tests += 3;
}

// Verifies safe copying for non-trivial and overlapping payload ranges.
void test_generic_payload_copying() {
	std::string values[]{ "a", "b", "c", "d" };
	Payload<std::string> payload(values, 4);
	Payload<std::string> alias = payload;
	alias.assign(payload, 0, 3, 1);
	require_refactor(payload[0] == "a" && payload[1] == "a"
		&& payload[2] == "b" && payload[3] == "c",
		"backward overlapping payload assignment failed");
	alias.assign(payload, 1, 3, 0);
	require_refactor(payload[0] == "a" && payload[1] == "b"
		&& payload[2] == "c" && payload[3] == "c",
		"forward overlapping payload assignment failed");
	Payload<std::string> independent = payload.copy();
	independent[0] = "changed";
	require_refactor(payload[0] == "a", "payload copy was not independent");
	passed_refactor_tests += 3;
}

// Verifies atomic payload reference counting across independent threads.
void test_thread_safe_reference_counting() {
	Payload<int> payload(32, 7);
	Tensor<int> tensor = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	const Tensor<int>& const_tensor = tensor;
	auto view = const_tensor[1];
	std::atomic<bool> valid{ true };
	std::vector<std::thread> threads;
	for (unsigned thread = 0; thread < 8; thread++) {
		threads.emplace_back([&payload, &view, &valid]() {
			for (unsigned i = 0; i < 10000; i++) {
				Payload<int> copy = payload;
				auto view_copy = view;
				if (copy[0] != 7 || copy[31] != 7
					|| view_copy.at({ 0 }) != 4 || view_copy.at({ 2 }) != 6)
					valid.store(false, std::memory_order_relaxed);
			}
		});
	}
	for (std::thread& thread : threads)
		thread.join();
	require_refactor(valid.load(std::memory_order_relaxed),
		"concurrent payload copies observed invalid data");
	require_refactor(payload.get_ref_count() == 1,
		"atomic payload reference count did not return to one");
	passed_refactor_tests += 2;
}

// Verifies dimension overflow detection and strong reshape exception safety.
void test_shape_safety() {
	size_t overflowing[]{ std::numeric_limits<size_t>::max(), 2 };
	bool threw = false;
	try {
		Tensor<int> invalid(overflowing, 2);
	}
	catch (const std::overflow_error&) {
		threw = true;
	}
	require_refactor(threw, "overflowing tensor dimensions did not throw");

	Tensor<int> tensor = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	Payload<int> before = tensor.get_data().copy();
	threw = false;
	try {
		tensor.reshape(overflowing, 2);
	}
	catch (const std::overflow_error&) {
		threw = true;
	}
	require_refactor(threw, "overflowing reshape dimensions did not throw");
	size_t dim[2]{};
	tensor.get_dim(dim);
	require_refactor(dim[0] == 2 && dim[1] == 3,
		"failed reshape changed tensor dimensions");
	Payload<int> after = tensor.get_data();
	for (size_t i = 0; i < tensor.get_size(); i++)
		require_refactor(after[i] == before[i], "failed reshape changed tensor data");
	passed_refactor_tests += 4;
}

// Verifies shared shape lifetime when an owner is reshaped after slicing.
void test_shared_view_shape() {
	Tensor<int> tensor = Tensor<int>::range(1, 7).reshape({ 2, 3 });
	Tensor<int> view = tensor[1];
	tensor.reshape({ 3, 2 });
	size_t dim[1]{};
	view.get_dim(dim);
	require_refactor(view.get_rank() == 1 && dim[0] == 3,
		"view shape changed with owner reshape");
	require_refactor(view.get_data()[3] == 4, "view payload changed with owner reshape");
	passed_refactor_tests += 2;
}

// Runs the focused ownership, safety, and view regression tests.
int main() {
	try {
		test_owned_assignment();
		test_view_assignment();
		test_const_views();
		test_generic_payload_copying();
		test_thread_safe_reference_counting();
		test_shape_safety();
		test_shared_view_shape();
		std::cout << "Tensor refactor tests passed: " << passed_refactor_tests << std::endl;
	}
	catch (const std::exception& error) {
		std::cerr << "Tensor refactor test failed: " << error.what() << std::endl;
		return 1;
	}
	return 0;
}
