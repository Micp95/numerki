#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "polynomial.h"

class Integral{
public:

	enum method {
		Rectangle,
		Trapezoidal,
		Simpson,
		Gauss,
		Newton
	};

	Integral(polynomial* pol, method myMethod);
	Integral(double(*fun)(double x), method myMethod);

	void setMethod(method newMethod);

	double calculate(double a, double b, int n);

private:
	polynomial* myPolynomial;
	double (*funPoint)(double x);

	bool isPolynomial;
	method myMethod;

	double PointValue(double x);

	double funRectangle(double a, double b, int n);
	double funTrapezoidal(double a, double b, int n);
	double funSimpson(double a, double b, int n);
	double funGauss(double a, double b, int n);
	double funNewton(double a, double b, int n);

	class GaussTab
	{
	public:
		GaussTab(int N);
		~GaussTab();

		int getSize();
		int getK(int k);
		double getXk(int k);
		double getAk(int k);


	private:
		int N;
		int*k;
		double* Xk;
		double* Ak;

	};


};



#endif