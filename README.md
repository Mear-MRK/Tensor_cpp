# C++ Tensor Library

A lightweight and flexible C++ library for multi-dimensional tensor operations, designed with a focus on efficient memory management and a clean API.

## Features

-   **Multi-dimensional Tensors:** Create and manipulate tensors of any rank.
-   **Efficient Memory Management:** Utilizes reference counting and copy-on-write semantics via a `Payload` class to minimize memory overhead.
-   **Rich Set of Operations:** Supports a wide range of operations including:
    -   Element-wise arithmetic (+, -, *, /)
    -   Scalar operations
    -   Tensor contraction and matrix product (`%`)
-   **Flexible Initialization:** Tensors can be created from dimensions, initializer lists, existing arrays, or generated with ranges and linearly spaced values.
-   **Slicing and Views:** Easily create views of tensors using indexing without copying data.
-   **Type Safety:** Templated `Tensor` class to support various data types (e.g., `float`, `int`, `std::complex`).
-   **Debugging Support:** In-built logging that can be enabled with the `DEBUG` preprocessor definition.

## Getting Started

### Prerequisites

-   A C++ compiler that supports C++11 or later (e.g., g++, clang++).

### Building and Running Tests

The project includes a test suite to verify the functionality of the `Tensor` library. To compile and run the tests:

```bash
g++ -std=c++11 -DDEBUG Tensor_test.cpp -o test_runner
./test_runner
```

This will compile the test file and run the executable. The `-DDEBUG` flag enables detailed logging from the library.

## Usage

Here are some examples of how to use the `Tensor` library.

### Creating Tensors

```cpp
#include "Tensor.hpp"
#include <iostream>

// Create a 2x3 tensor of doubles, initialized to 1.0
Tensor<double> tensor1({2, 3}, 1.0);

// Create a rank-1 tensor (vector) of 5 complex numbers
Tensor<std::complex<float>> tensor2(5);

// Create a 3x3 tensor of random integers
Tensor<int> tensor3({3, 3}, rand);

// Create a tensor with a range of values [1-24] and reshape it
Tensor<> tensor4 = Tensor<int>::range(1, 25).reshape({2, 2, 3, 2});

// Create a tensor with 12 linearly spaced values between 1 and 5
Tensor<> tensor5 = Tensor<float>::linspace(1, 5, 12);
```

### Operations

```cpp
// Create two 2x2 tensors
Tensor<float> a({2, 2}, 1.0f);
Tensor<float> b({2, 2}, 2.0f);

// Element-wise addition
Tensor<float> c = a + b;

// Scalar multiplication
c *= 2.0f;

// Tensor contraction (matrix product in this case)
auto res = a % b;

// Print the result
std::cout << res.to_string() << std::endl;
```

### Slicing and Views

```cpp
Tensor<int> t({3, 3, 3}, 1, 27);

// Get a view of the second slice along the first dimension
Tensor<int> view = t[1];

// Modify the view (this will modify the original tensor)
view[0] = 99;

std::cout << t.to_string() << std::endl;
```

## Code Structure

-   `Tensor.hpp`: The main header file containing the `Tensor` class implementation.
-   `Payload.hpp`: Contains the `Payload` class, which manages the underlying data array with reference counting.
-   `Tensor_test.cpp`: The test suite for the library.
-   `README.md`: This file.

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue.