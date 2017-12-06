#include <iostream>
#include "svm.h"

int main(int argc, char **argv)
{
	uword n = 100000;
	uword m = 400;
	mat X = randu<mat>(n, m);
	mat y = randu<vec>(m);
	mat A(n, m);

	KRBF K(&X, 0.1);

	auto start = std::chrono::system_clock::now();
	dot(y, K);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	cout << "Elapse: " << elapsed.count() << " seconds" << endl;

	auto start2 = std::chrono::system_clock::now();
	dot(y, K);
	auto end2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - start2;
	cout << "Elapse: " << elapsed2.count() << " seconds" << endl;

//#pragma omp parallel for
	for (int i = 0; i < K.dim; i++) {
		for (int j = 0; j < K.dim; j++) {
			A(i,j) = K(i, j);
		}
	}

	auto start3 = std::chrono::system_clock::now();
	dot(y, A.t());
	auto end3 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed3 = end3 - start3;
	cout << "Elapse: " << elapsed3.count() << " seconds" << endl;


//	auto start = std::chrono::system_clock::now();
//#pragma omp parallel for //collapse(2)
//	for (int i = 0; i < K.dim; i++) {
//		for (int j = 0; j < K.dim; j++) {
//			A(i,j) = K(i, j);
//		}
//	}
//	auto end = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed = end - start;
//	cout << "Elapse: " << elapsed.count() << " seconds" << endl;
//
//	auto start2 = std::chrono::system_clock::now();
//#pragma omp parallel for //collapse(2)
//	for (int i = 0; i < K.dim; i++) {
//		for (int j = 0; j < K.dim; j++) {
//			auto tmp = A(i, j);
//		}
//	}
//	auto end2 = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed2 = end2 - start2;
//	cout << "Elapse: " << elapsed2.count() << " seconds" << endl;
//
//	auto start3 = std::chrono::system_clock::now();
//#pragma omp parallel for //collapse(2)
//	for (int i = 0; i < K.dim; i++) {
//		for (int j = 0; j < K.dim; j++) {
//			auto tmp = K(i, j);
//		}
//	}
//	auto end3 = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed3 = end3 - start3;
//	cout << "Elapse: " << elapsed3.count() << " seconds" << endl;

	system("pause");
	return 0;
}