#include <iostream>
#include "svm.h"

int main(int argc, char **argv)
{
	uword n = 100000;
	uword m = 400;
	mat X = randu<mat>(n, m);
	vec y = randu<vec>(m);

	KRBF K(&X, 0.1);

	auto start = std::chrono::system_clock::now();
	auto r1 = K * y;
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	cout << "Elapse: " << elapsed.count() << " seconds" << endl;

	auto start2 = std::chrono::system_clock::now();
	auto r2 = K * y;
	auto end2 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed2 = end2 - start2;
	cout << "Elapse: " << elapsed2.count() << " seconds" << endl;

	mat A = K.copy();
	auto start3 = std::chrono::system_clock::now();
	auto r3 = A * y;
	auto end3 = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed3 = end3 - start3;
	cout << "Elapse: " << elapsed3.count() << " seconds" << endl;

	system("pause");
	return 0;
}