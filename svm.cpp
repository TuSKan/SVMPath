#include <iostream>
#include "svm.h"

//int main(int argc, char **argv)
//{
//	mat X = randu<mat>(5, 5);
//	X.print("X: ");
//
//	KRBF K(&X, 0.1);
//
//	cout << endl;
//	cout << K(0, 0) << "\t" << K(0, 1) << "\t" << K(0, 2) << "\t" << K(0, 3) << "\t" << K(0, 4) << "\t" << endl;
//	cout << K(1, 0) << "\t" << K(1, 1) << "\t" << K(1, 2) << "\t" << K(1, 3) << "\t" << K(1, 4) << "\t" << endl;
//	cout << K(2, 0) << "\t" << K(2, 1) << "\t" << K(2, 2) << "\t" << K(2, 3) << "\t" << K(2, 4) << "\t" << endl;
//	cout << K(3, 0) << "\t" << K(3, 1) << "\t" << K(3, 2) << "\t" << K(3, 3) << "\t" << K(3, 4) << "\t" << endl;
//	cout << K(4, 0) << "\t" << K(4, 1) << "\t" << K(4, 2) << "\t" << K(4, 3) << "\t" << K(4, 4) << "\t" << endl;
//	cout << endl;
//
//	cout << "Cache size: " << K.Cache()->size() << endl;
//
//	system("pause");
//	return 0;
//}
//
int main(int argc, char **argv)
{
	mat X = randu<mat>(1000000, 500);
	mat A(X.n_cols, X.n_cols);

	KRBF K(&X, 0.1);

	auto start = std::chrono::system_clock::now();
#pragma omp parallel for //collapse(2)
	for (int i = 0; i < K.ncols; i++) {
		for (int j = 0; j < K.ncols; j++) {
			A(i,j) = K(i, j);
		}
	}
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	cout << "Elapse: " << elapsed.count() << " seconds" << endl;

	auto start2 = std::chrono::system_clock::now();
#pragma omp parallel for //collapse(2)
	for (int i = 0; i < K.ncols; i++) {
		for (int j = 0; j < K.ncols; j++) {
			auto tmp = A(i, j);
		}
	}
	auto end2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - start2;
	cout << "Elapse: " << elapsed2.count() << " seconds" << endl;

	auto start3 = std::chrono::system_clock::now();
#pragma omp parallel for //collapse(2)
	for (int i = 0; i < K.ncols; i++) {
		for (int j = 0; j < K.ncols; j++) {
			auto tmp = K(i, j);
		}
	}
	auto end3 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed3 = end3 - start3;
	cout << "Elapse: " << elapsed3.count() << " seconds" << endl;

	system("pause");
	return 0;
}