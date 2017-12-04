#include <armadillo>
#include "LRUCache.hpp"

using namespace std;
using namespace arma;

class KFun {
public:
	// Constructor
	KFun(const mat * X = NULL, const string name = NULL, const double param = NULL) : X(X), name(name), param(param), ncols(X->n_cols) {};
	~KFun() {};
	// Operator
	virtual double operator() (const unsigned int i, const unsigned int j) {};
	// Accessors
	const lru::Cache<unsigned int, double, mutex>* Cache() { return &KCache; }
	// Members
protected:
	string name;
	double param;
	lru::Cache<unsigned int, double, mutex> KCache;
	const mat * X;
public:
	const unsigned int ncols;
}; // KFun

   // Linear Kernel
class KLinear : public KFun {
public:
	KLinear(const mat * X) : KFun(X, "Linear") {}
	double operator() (const unsigned int i, const unsigned int j) {
		double res;
		unsigned int IJ = (i < j ? i : j) + (i > j ? i : j) * ncols;
		if (!KCache.tryGet(IJ, res)) {
			res = dot(X->col(i), X->col(j));
			KCache.insert(IJ, res);
		}
		return res;
	}
}; // KLinear

   // Gaussian Kernel
class KRBF : public KFun {
public:
	KRBF(const mat * X, double gamma) : KFun(X, "RBF", gamma) {}
	double operator() (const unsigned int i, const unsigned int j) {
		double res;
		unsigned int IJ = (i < j ? i : j) + (i > j ? i : j) * ncols;
		if (!KCache.tryGet(IJ, res)) {
			res = exp(-param*(dot(X->col(i), X->col(i)) + dot(X->col(j), X->col(j)) - 2 * dot(X->col(i), X->col(j))));
			KCache.insert(IJ, res);
		}
		return res;
	}
}; // KRBF

   // Polynomial Kernel
class KPolynomial : public KFun {
private:
	unsigned int degree;
	double coef0;
public:
	KPolynomial(const double gamma, const unsigned int degree, const double coef0 = NULL) : KFun(X, "Polynomial", gamma), degree(degree), coef0(coef0) {}
	double operator() (const unsigned int i, const unsigned int j) {
		double res;
		unsigned int IJ = (i < j ? i : j) + (i > j ? i : j) * ncols;
		if (!KCache.tryGet(IJ, res)) {
			res = pow(param * dot(X->col(i), X->col(j)) + coef0, degree);
			KCache.insert(IJ, res);
		}
		return res;
	}
};

// Sigmoid kernel
class KSigmoid : public KFun {
private:
	double coef0;
public:
	KSigmoid(double gamma, double coef0 = NULL) : coef0(coef0) { param = gamma; name = "Sigmoid"; }
	double operator() (const unsigned int i, const unsigned int j) {
		double res;
		unsigned int IJ = (i < j ? i : j) + (i > j ? i : j) * ncols;
		if (!KCache.tryGet(IJ, res)) {
			res = tanh(param * dot(X->col(i), X->col(j)) + coef0);
			KCache.insert(IJ, res);
		}
		return res;
	}
};

class Kernel {
public:
	Kernel() {};
	~Kernel() {};
	
public:

}; // Kernel


