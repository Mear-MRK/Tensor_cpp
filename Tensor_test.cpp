#include <iostream>
#include <random>
#include <functional>
#include <complex>
#include <cstdlib>

#include "tensor.hpp"


std::default_random_engine rnd_gen;
std::uniform_real_distribution<float> flt_dist(-1.0f, 1.0f);

float rnd(void) {
	return flt_dist(rnd_gen);
}

template<typename T>
void print_arr(T arr[], size_t size) {
	for (size_t i = 0; i < size; i++)
		std::cout << arr[i] << " ";
	std::cout << std::endl;
}

float sqr(float f) {
	return f * f;
}

void test_basics() {
	size_t dim[] = { 2, 3 };
	unsigned rank = 2;
	Tensor<float> tn{ dim, rank, 1.11111f };
	tn *= 3.24f;
	Tensor<> r_tn(dim, rank, rnd);
	std::cout << tn.to_string() << std::endl;
	std::cout << r_tn.to_string() << std::endl;
	Tensor<> res{ {2,3} };
	res = tn + r_tn;
	res.reshape({ 3, 2 }).apply(sqr);
	std::cout << res.to_string() << std::endl;
	Tensor<std::complex<float>> ct(3);
	ct[1] = std::complex<float>{-1, 1};
	std::cout << ct.to_string() << std::endl;
	ct += std::complex<float>{0,3};
	std::cout << ct.to_string() << std::endl;
	ct /= std::complex<float>{1,1};
	std::cout << ct.to_string() << std::endl;

}

void test_cnt() {

	Tensor<int> t1("t1");
	t1 = Tensor<int>::range(1, 25).reshape({ 2,2,3,2 });
	std::cout << "t1:\n" << t1.to_string() << std::endl;
	int v[] = { 11,12,21,22 };
	Tensor<int> t2({ 2, 2 }, v); t2.set_meta("t2");
	std::cout << "t2:\n" << t2.to_string() << std::endl;
	std::cout << "cnt t1,t2:\n";
	Tensor<int> t = Tensor<int>::cnt(t1, t2, 1, 1);
	std::cout << t.to_string() << '\n' << std::endl;

	Tensor<float> l_tf({ 2, 2 }, 1.0f, 4.0f); 
	std::cout << "l_tf:\n" << l_tf.to_string() << std::endl;
	Tensor<float> r_tf = std::move(Tensor<float>::range(0.0f, 4.0f).reshape({ 2, 2 }));
	std::cout << "r_tf:\n" << r_tf.to_string() << std::endl;
	std::cout << "cnt l,r:\n" << Tensor<float>::cnt(l_tf, r_tf, 1, 0).to_string() << '\n' << std::endl;

	Tensor<float> l_v({ 3 }, 1.0f, 3.0f);
	std::cout << l_v.to_string() << std::endl;
	Tensor<float> r_v = Tensor<float>::range(-2, 1).reshape({ 3 });
	std::cout << r_v.to_string() << std::endl;
	std::cout << Tensor<float>::cnt(l_v, r_v, 0, 0).to_string() << std::endl;
	std::cout << (l_v % r_v).to_string() << (r_v % l_v).to_string() << std::endl;

	std::cout << "--End cnt--" << std::endl;
}

void test_index() {
	size_t i1[]{ 1,2,3,4 }; print_arr<size_t>(i1, 4);
	size_t i2[]{ 5,6,7 };   print_arr<size_t>(i2, 3);
	size_t i[5]{ 0 };
	unsigned li = 1;
	unsigned ri = 2;
	std::cout << li << ", " << ri << std::endl;
	comp_index(i1, i2, i, 4, 3, li, ri);
	print_arr<size_t>(i, 5);
	decomp_index(i, i1, i2, 4, 3, 0, 1);
	i1[0] = 22;
	i2[1] = 22;
	print_arr<size_t>(i1, 4);
	print_arr<size_t>(i2, 3);
	print_arr<size_t>(i, 5);
}

void test_initlist(const std::initializer_list<size_t>& dim_l) {
	unsigned rank = dim_l.size();
	std::cout << "rank: " << rank << std::endl;
	size_t* jump = new size_t[rank];
	size_t* dim = new size_t[rank];
	size_t size = (bool)rank;
	unsigned i = 0;
	for (size_t d : dim_l) {
		dim[i] = d;
		size *= d;
		i++;
	}
	jump[rank - 1] = 1;
	for (unsigned i = 1; i < rank; i++) {
		jump[rank - i - 1] = dim[rank - i] * jump[rank - i];
	}
	print_arr<size_t>(dim, rank);
	print_arr<size_t>(jump, rank);
	delete[] jump;
	delete[] dim;
}

void test_opt() {
	Tensor<> a({ 3, 4 }, rnd); // Tensor<>::linspace(1, 5, 12).reshape({ 3,4 });
	Tensor<> b({ 3, 4 }, rnd);
	std::cout << a.to_string() << b.to_string() << std::endl;

	Tensor<> c = a + b ;
	std::cout << c.to_string() << std::endl;

}

void test_views() {
	Tensor<int> t({ 3, 3, 3 }, 1, 27);
	std::cout << "t: " << t.to_string() << '\n';
	Tensor<int> tt = t[1] = 33;
	std::cout << "t: " << t.to_string() << '\n';
	t[{0, 2}] = Tensor<int>::linspace(66,68,3);
	std::cout << "t: " << t.to_string() << '\n';
	t[{1, 1, 0}] = 99;
	std::cout << "t: " << t.to_string() << std::endl;
	std::cout << "tt: " << tt.to_string() << std::endl;
}

int main() {

	test_basics();
	test_initlist({ 3,2,3,2 });
	test_index();
	test_cnt();
	test_opt();
	test_views();

	std::cout << "Press Enter to end>";
	std::cin.get();
	return 0;
}