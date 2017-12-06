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
	// Operator
	virtual double operator() (const uword i, const uword j) {};
	// Accessors
	const bool isFilled() const { return count == size; }
	const vec& getMtri() const { return MTri; }
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

const vec dot(const vec& v, KMatrix& K) {
	arma_assert_same_size(v.n_rows, v.n_cols, K.dim, 1, string("dot KMatrix").c_str());
	vec res(K.dim);
	if (K.isFilled()) {
		const vec& MTri = K.getMtri();
//#pragma omp parallel for
		for (uword j = 0; j < K.dim; ++j) {
			uvec ID = zeros<uvec>(K.dim);
			for (uword i = 0; i < j; ++i) ID(i) = j + i * K.dim - accu(sequence(0, i));
			ID(sequence(j, K.dim - 1)) = sequence(j + j * K.dim - accu(sequence(0, j)), K.dim - 1 + j * K.dim - accu(sequence(0, j)));
			res(j) = dot(v, MTri.elem(ID));
		}
	}
	else {		
//#pragma omp parallel for
		for (uword j = 0; j < K.dim; ++j) {
			double tmp(0);
			for (uword i = 0; i < K.dim; ++i) {
				tmp += v[i] * K(i, j);
			}
			res(j) = tmp;
		}
	}
	return res;
}

const vec dot(KMatrix& K, const vec& v) {  return dot(v, K); }

   // Linear Kernel
class KLinear : public KMatrix {
public:
	KLinear(const mat * X) : KMatrix(X, "Linear") {}
	double operator() (const uword i, const uword j) {
		uword IJ;
		if (i < j) IJ = j + i * dim - accu(sequence(0, i));
		else IJ = i + j * dim - accu(sequence(0, j));
		double res = MTri[IJ];
		if (isnan(res)) {
			res = dot(X->col(i), X->col(j));
			MTri[IJ] = res;
			++count;
		}
		return res;
	}
}; // KLinear

   // RBF Kernel
class KRBF : public KMatrix {
public:
	KRBF(const mat * X, double gamma) : KMatrix(X, "RBF", gamma) {}
	double operator() (const uword i, const uword j) {
		uword IJ;
		if (i < j) IJ = j + i * dim - accu(sequence(0,i));
		else IJ = i + j * dim - accu(sequence(0, j));
		double res = MTri[IJ];
		if (isnan(res)) {
			res = exp(-param*(dot(X->col(i), X->col(i)) + dot(X->col(j), X->col(j)) - 2 * dot(X->col(i), X->col(j))));
			MTri[IJ] = res;
			++count;
		}
		return res;
	}
}; // KRBF

   // Polynomial Kernel
class KPolynomial : public KMatrix {
private:
	uword degree;
	double coef0;
public:
	KPolynomial(const double gamma, const uword degree, const double coef0 = NULL) : KMatrix(X, "Polynomial", gamma), degree(degree), coef0(coef0) {}
	double operator() (const uword i, const uword j) {
		uword IJ;
		if (i < j) IJ = j + i * dim - accu(sequence(0, i));
		else IJ = i + j * dim - accu(sequence(0, i));
		double res = MTri[IJ];
		if (isnan(res)) {
			res = pow(param * dot(X->col(i), X->col(j)) + coef0, degree);
			MTri[IJ] = res;
			++count;
		}
		return res;
	}
};

// Sigmoid kernel
class KSigmoid : public KMatrix {
private:
	double coef0;
public:
	KSigmoid(double gamma, double coef0 = NULL) : coef0(coef0) { param = gamma; name = "Sigmoid"; }
	double operator() (const uword i, const uword j) {
		uword IJ;
		if (i < j) IJ = j + i * dim - accu(sequence(0, i));
		else IJ = i + j * dim - accu(sequence(0, i));
		double res = MTri[IJ];
		if (isnan(res)) {
			res = tanh(param * dot(X->col(i), X->col(j)) + coef0);
			MTri[IJ] = res;
			++count;
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


