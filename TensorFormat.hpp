#pragma once
#include <complex>
#include <string>

// Formats a complex value for tensor display.
template<typename F>
inline std::string to_str(std::complex<F> c) {
	return std::to_string(c.real()) + "+" + std::to_string(c.imag()) + "i";
}

// Formats a scalar value for tensor display.
template<typename F>
inline std::string to_str(F f) {
	return std::to_string(f);
}
