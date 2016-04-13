#include "Integral.h"


Integral::Integral(polynomial * pol, method myMethod){
	this->myMethod = myMethod;
	this->myPolynomial = pol;
	isPolynomial = true;
}

Integral::Integral(double(*fun)(double x), method myMethod){
	this->myMethod = myMethod;
	this->funPoint = fun;
	isPolynomial = false;
}

void Integral::setMethod(method newMethod){
	this->myMethod = newMethod;
}

double Integral::PointValue(double x){
	if (isPolynomial)
		return myPolynomial->PointValue(x);
	else
		return funPoint(x);
}

double Integral::calculate(double a, double b, int n){
	switch (myMethod){
	case Integral::Rectangle:
		return funRectangle(a, b, n);
	case Integral::Trapezoidal:
		return funTrapezoidal(a, b, n);
	case Integral::Simpson:
		return funSimpson(a, b, n);
	default:
		return 0.0;
	}
}


double Integral::funRectangle(double a, double b, int n){
	double deltaX = (b - a) / (double)n;
	double res = 0,x = a;

	for (int k = 0; k < n; k++) {
		res += deltaX*PointValue(x);
		x += deltaX;
	}

	return res;
}

double Integral::funTrapezoidal(double a, double b, int n){
	double deltaX = (b - a) / (double)n;
	double res = 0, xa = a,xb =a;

	for (int k = 0; k < n; k++) {
		xb += deltaX;
		res += deltaX*((PointValue(xa)+PointValue(xb))/2);
		xa = xb;
	}
	return res;
}

double Integral::funSimpson(double a, double b, int n){
	double deltaX = (b - a) / (double)n;
	double res = 0, xa = a, xb=a;

	for (int k = 0; k < n; k++) {
		xb += deltaX;
		res += (deltaX/6)*(PointValue(xa) + 4*PointValue((xa+xb)/2) + PointValue(xb));
		xa = xb;
	}
	return res;
}


