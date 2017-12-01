#include <armadillo>
#include "LRUCache.hpp"

using namespace std;
using namespace arma;


class Kernel {
protected:
	class KFun {
	public:
		// Constructor
		KFun() {};
		~KFun() {};
		// Operator
		virtual double operator() (const vec &Xi, const vec &Xj) {};
		// Members
	protected:
		string name;
		double param;
	}; // KFun

	// Linear Kernel
	class KLinear : public KFun {
	public:
		KLinear() { name = "Linear"; }
		double operator() (const vec &Xi, const vec &Xj) { 
			return dot(Xi, Xj); 
		}
	}; // KLinear

	// Gaussian Kernel
	class KRBF : public KFun {
	public:
		KRBF(double gamma = NULL) { param = gamma; name = "RBF"; }
		double operator() (const vec &Xi, const vec &Xj) { 
			return exp(-param*(dot(Xi,Xi) + dot(Xj, Xj) - 2 * dot(Xi, Xj)));
		}
	}; // KRBF

	// Polynomial Kernel
	class KPolynomial : public KFun {
	private:
		unsigned int degree;
		double coef0;
	public:
		KPolynomial(double gamma = NULL, unsigned int degree = NULL, double coef0 = NULL) : degree(degree), coef0(coef0) { param = gamma; name = "Polynomial"; }
		double operator() (const vec &Xi, const vec &Xj) {
			return pow(param * dot(Xi, Xj) + coef0, degree);
		}
	};

	// Sigmoid kernel
	class KSigmoid : public KFun {
	private:
		double coef0;
	public:
		KSigmoid(double gamma = NULL, double coef0 = NULL) : coef0(coef0) { param = gamma; name = "Sigmoid"; }
		double operator() (const vec &Xi, const vec &Xj) {
			return tanh(param * dot(Xi, Xj) + coef0);
		}
	};

public: 
	Kernel(mat * const * X, vec const * Y, const KFun &K);
	~Kernel();
}; // Kernel


