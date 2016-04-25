#ifndef DIFFERENTIAL_H
#define DIFFERENTIAL_H

#include "Point.h"

class differential
{
public:
	enum Method {
		Euler, RungegKuttyIV, RungegKuttyII
	};

	differential(double(*function)(double x, double y),Method method) : function(function),myMethod(method) {};

	void setMethod(Method newMethod) { myMethod = newMethod; };
	double calculateDiffrential(double a, double b, Point start, double h, double x, int n = -1);

private:
	Method myMethod;

	double diffrentialEuler(double a, double b, Point start, double h,double x);
	double diffrentialRugegKuttIV(double a, double b, Point start, double h, double x);
	double diffrentialRugegKuttII(double a, double b, Point start, double h, double x);

	double(*function)(double x, double y);
	int compartment(double a, double b,double h);
	int n;

	template<class T>
	double diffrentialUniversalFunction(double a, double xx, Point start, double h, T fun);	
};

template<class T>
double differential::diffrentialUniversalFunction(double a,double xx, Point start, double h, T fun){
	double px = start.x, x = px;
	double py = start.y, y = py;

	int n;
	if (this->n == -1)
		n = compartment(a, xx, h);
	else
		n = this->n;

	for (int k = 0; k < n; k++) {
		x = px + h;
		y = py + h*fun(px, py);

		px = x;
		py = y;
	}

	return y;
}

#endif
