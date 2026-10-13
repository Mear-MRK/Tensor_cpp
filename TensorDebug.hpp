#pragma once

#ifdef TENSOR_DEBUG
#include <iostream>
#define TENSOR_LOG(message) do { std::cout << message << std::endl; } while (false)
#else
#define TENSOR_LOG(message) do {} while (false)
#endif
