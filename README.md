# Tensor Library in C++

This project provides an implementation of a Tensor library in C++. It includes handling of multi-dimensional arrays with support for various operations and memory management techniques. The library also includes reference counting for efficient memory use and copy-on-write semantics to optimize performance.


## Features

- **Tensor Class**: Supports multi-dimensional arrays with various mathematical operations.
- **Payload Class**: Manages the underlying data with reference counting and efficient memory operations.
- **Memory Management**: Includes copy constructors, move constructors, and assignment operators for efficient memory handling.
- **Debugging**: Debugging support with logging enabled via the `DEBUG` preprocessor definition.


## Usage

You can create a Tensor with various constructors provided:

```cpp
#include <iostream>
#include "Tensor.hpp"

int main() {
    Tensor<double> tensor1({2, 3, 4}); // Create a 2x3x4 tensor (rank-3)
    Tensor<float> tensor2(5); // Create a one-dimensional tensor of size 5
    Tensor<> tensor3 = Tensor<int>::range(1, 25).reshape({ 2,2,3,2 }); // Create an int tensor with a range of values [1-24] then reshape it to 2x2x3x2.
    Tensor<> tensor4 = Tensor<float>::linspace(1, 5, 12);  // Create a one-dimensional tensor containing 12 linearly spaced numbers between 1 and 5 (inclusive).

    std::cout << tensor3.to_string() << std::endl;
    return 0;
}
```

The library supports various operations like arithmetic operations and **tensor contraction**:

```cpp
Tensor<float> tensor1({2, 2}, 1.0); // 2x2 tensor initialized with 1.0
Tensor<float> tensor2({2, 2}, -2.0); // 2x2 tensor initialized with 2.0

Tensor<float> result = tensor1 + tensor2; // Addition
result *= 2.0; // Multiplication by scalar

result = Tensor<float>::cnt(tensor1, tensor2, 1, 0) // Contraction of 2nd index of tensor1 with 1st index of tensor2

result.apply(func)  // apply element-wise float func(float) to result
```

The library includes support for copy and move semantics to optimize performance:

```cpp
Tensor<float> tensor1({2, 2}, 1.0);
Tensor<float> tensor2 = tensor1; // Copy constructor
Tensor<float> tensor3 = std::move(tensor1); // Move constructor
```

Enable debugging by defining the `DEBUG` preprocessor directive:

```cpp
#define DEBUG 1
#include "Tensor.hpp"
```
