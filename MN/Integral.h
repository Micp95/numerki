#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "polynomial.h"

class Integral{
public:

	enum method {
		Rectangle,
		Trapezoidal,
		Simpson
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
};



#endif