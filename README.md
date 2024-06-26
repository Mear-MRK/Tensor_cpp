# Tensor Library in C++

This project provides an implementation of a Tensor library in C++. It includes handling of multi-dimensional arrays with support for various operations and memory management techniques. The library also includes reference counting for efficient memory use and copy-on-write semantics to optimize performance.


## Features

- **Tensor Class**: Supports multi-dimensional arrays with various mathematical operations.
- **Payload Class**: Manages the underlying data with reference counting and efficient memory operations.
- **Memory Management**: Includes copy constructors, move constructors, and assignment operators for efficient memory handling.
- **Debugging**: Debugging support with logging enabled via the `DEBUG` preprocessor definition.

## Important Note
The provided Tensor library, while functional, is not *yet* optimized to utilize CPU intrinsics and multi-thread parallelism. Additionally, it does not currently employ techniques such as rearranging of the elements and loop optimizations for efficient utilization of CPU cache. Performance improvements in these areas require further development.

## Usage

You can create a Tensor with various constructors provided:

```cpp
#include <iostream>
#include <cstdlib>
#include <complex>
#include "Tensor.hpp"

int main() {

    Tensor<double> tensor1({2, 3, 4}, 1); // Creates a 2x3x4 tensor (rank-3) of doubles initialized with 1.

    Tensor<std::complex<float>> tensor2(5); // Creates a one-dimensional tensor of size 5 of complex numbers.

    Tensor<int> tensor3({3, 3}, rand); // Creates a 3x3 tensor of random integers generated by 'int rand(void)' function.

    Tensor<> tensor4 = Tensor<int>::range(1, 25).reshape({2, 2, 3, 2}); // Creates a tensor with a range of values [1-24] and reshapes it to 2x2x3x2.

    Tensor<> tensor5 = Tensor<float>::linspace(1, 5, 12); // Creates a one-dimensional tensor containing 12 linearly spaced numbers between 1 and 5 (inclusive).

    std::cout << tensor3.to_string() << std::endl; // Prints tensor3.

    return 0;
}
```

The library supports various operations like arithmetic operations and **tensor contraction**:

```cpp
Tensor<float> tensor1({2, 2}, 1); // 2x2 tensor initialized with 1
Tensor<float> tensor2({2, 2}, -2); // 2x2 tensor initialized with -2

tensor1[{1}] = 3; // Sets the 2nd row of tensor1 to 3

Tensor<float> result = tensor1 + tensor2; // Addition
result *= 2.0; // Multiplication by scalar

auto res = Tensor<float>::cnt(tensor1, tensor2, 1, 0); // Contraction of 2nd index of tensor1 with 1st index of tensor2
auto res2 = tensor1 % tensor2; // Matrix product
bool b = res.is_close_to(res2, 1e-7); // True

auto res3 = result.copy(); // Cloning result

res3.apply(func);  // Applying element-wise function 'float func(float)' to res3
```

Enable debugging by defining the `DEBUG` preprocessor directive:

```cpp
#define DEBUG 1
#include "Tensor.hpp"
```
