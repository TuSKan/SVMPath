#include <armadillo>

using namespace std;
using namespace arma;

const uvec sequence(const uword& begin, const uword& end) {
	uvec r(end - begin + 1);
	for (uword i = begin; i < end + 1; ++i)
		r[i - begin] = i;
	return r;
}


class KMatrix {
public:
	// Constructor
	KMatrix(const mat * X = NULL, const string name = NULL, const double param = NULL) 
		: X(X), name(name), param(param), dim(X->n_cols), size((X->n_cols + X->n_cols*X->n_cols) / 2), count(0) {
		MTri.resize(size);
		MTri.fill(numeric_limits<double>::quiet_NaN());
	};
	~KMatrix() { };
	// Virtual
	virtual double KFun(const uword i, const uword j) const = 0;
	// Operator
	double operator() (const uword i, const uword j) {
		uword IJ;
		if (i < j) IJ = j + i * dim - accu(sequence(0, i));
		else IJ = i + j * dim - accu(sequence(0, j));
		if (isnan(MTri[IJ])) MTri[IJ] = KFun(i, j);
		return MTri[IJ];
	}

	const vec operator* (const vec& v) {
		arma_assert_mul_size(dim, dim, v.n_rows, v.n_cols, string("multiply KMatrix").c_str());
		if (count != size) fill(0, dim);
		vec res(dim);
#pragma omp parallel for
		for (uword j = 0; j < dim; ++j) {
			res(j) = dot(v, col(j));
		}
		return res;
	}
	// Accessors
	void fill(const uword begin, const uword end) {
#pragma omp parallel for
		for (uword j = begin; j < end; ++j) {
			for (uword i = j; i < dim; ++i) {
				const uword IJ = i + j * dim - accu(sequence(0, j));
				if (isnan(MTri[IJ])) MTri[IJ] = KFun(i, j);
			}
		}
	}
	const vec col(const uword j) const {
		uvec ID = zeros<uvec>(dim);
		for (uword i = 0; i < j; ++i) ID(i) = j + i * dim - accu(sequence(0, i));
		ID(sequence(j, dim - 1)) = sequence(j + j * dim - accu(sequence(0, j)), dim - 1 + j * dim - accu(sequence(0, j)));
		return MTri.elem(ID); 
	}
	const vec row(const uword i) const { return col(i); }
	const mat copy() {
		if (count != size) fill(0,dim);
		mat R(dim, dim);
#pragma omp parallel for
		for (uword j = 0; j < dim; ++j) {
			R.col(j) = col(j);
		}
		return R;
	}
	// Members
protected:
	string name;
	double param;
	vec MTri;
	const uword size;
	const mat * X;
	uword count;
public:
	const uword dim;
}; // KFun

   // Linear Kernel
class KLinear : public KMatrix {
public:
	KLinear(const mat * X) : KMatrix(X, "Linear") {}
	double KFun(const uword i, const uword j) const {
		return dot(X->col(i), X->col(j));
	}
}; // KLinear

   // RBF Kernel
class KRBF : public KMatrix {
public:
	KRBF(const mat * X, double gamma) : KMatrix(X, "RBF", gamma) {}
	double KFun(const uword i, const uword j) const {
		return exp(-param * (dot(X->col(i), X->col(i)) + dot(X->col(j), X->col(j)) - 2 * dot(X->col(i), X->col(j))));
	}
}; // KRBF

   // Polynomial Kernel
class KPolynomial : public KMatrix {
private:
	uword degree;
	double coef0;
public:
	KPolynomial(const double gamma, const uword degree, const double coef0 = NULL) : KMatrix(X, "Polynomial", gamma), degree(degree), coef0(coef0) {}
	double KFun (const uword i, const uword j) const {
		return pow(param * dot(X->col(i), X->col(j)) + coef0, degree);
	}
};

// Sigmoid kernel
class KSigmoid : public KMatrix {
private:
	double coef0;
public:
	KSigmoid(double gamma, double coef0 = NULL) : coef0(coef0) { param = gamma; name = "Sigmoid"; }
	double operator() (const uword i, const uword j) {
		return tanh(param * dot(X->col(i), X->col(j)) + coef0);
	}	
};

class Kernel {
public:
	Kernel() {};
	~Kernel() {};
	
public:

}; // Kernel


