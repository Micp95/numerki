#include "differential.h"


double differential::calculateDiffrential(double a, double b, Point start, double h, double x, int n ){

	if (n != -1)
		h = (x - a) / (double)n;
	this->n = n;
	
	switch (myMethod){
	case differential::Euler:
		return diffrentialEuler(a, b, start, h, x);
	case differential::RungegKuttyIV:
		return diffrentialRugegKuttIV(a, b, start, h, x);
	case differential::RungegKuttyII:
		return diffrentialRugegKuttII(a, b, start, h, x);
	case differential::Heuna:
		return diffrentialHeun(a, b, start, h, x);
	case differential::modEulera:
		return diffrentialModEuler(a, b, start, h, x);
	default:
		break;
	}

	return 0.0;
}

double differential::diffrentialEuler(double a, double b, Point start, double h, double xx){

	if (xx > b)
		return -1;

	return diffrentialUniversalFunction(a, xx, start, h, function);
}

double differential::diffrentialRugegKuttIV(double a, double b, Point start, double h, double xx){

	double* k = new double[4];
	double(*fun)(double x, double y) = function;

	auto fi = [&k,&h,&fun](double x, double y) -> double {
		k[0] = fun(x, y);
		k[1] = fun(x + 0.5*h, y + 0.5*h*k[0]);
		k[2] = fun(x + 0.5*h, y + 0.5*h*k[1]);
		k[3] = fun(x + h, y + h*k[2]);
		return (1.0/6.0)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
	};

	double y = diffrentialUniversalFunction(a, xx, start, h, fi);
	delete[] k;

	return y;
}

double differential::diffrentialRugegKuttII(double a, double b, Point start, double h, double xx){

	double* k = new double[2];
	double(*fun)(double x, double y) = function;

	auto fi = [&k, &h, &fun](double x, double y) -> double {
		k[0] = fun(x, y);
		k[1] = fun(x + h, y + h*k[0]);
		return (1.0 / 2.0)*(k[0] +  k[1]);
	};

	double y = diffrentialUniversalFunction(a,xx,start,h,fi);
	delete[] k;

	return y;
}

double differential::diffrentialHeun(double a, double b, Point start, double h, double xx)
{
	double(*fun)(double x, double y) = function;

	auto fi = [&h, &fun](double x, double y) -> double {
		
		return (0.5)*(fun(x, y) + fun(x + h*y, y + h*fun(x, y)));
	};

	double y = diffrentialUniversalFunction(a, xx, start, h, fi);

	return y;
}

double differential::diffrentialModEuler(double a, double b, Point start, double h, double xx)
{
	double(*fun)(double x, double y) = function;

	auto fi = [&h, &fun](double x, double y) -> double {

		return fun(x + 0.5*h,y+0.5*h*fun(x,y));
	};

	double y = diffrentialUniversalFunction(a,	 xx, start, h, fi);

	return y;
}

int differential::compartment(double a, double b, double h) {
	double N = (b - a) / h;
	int n = (int)N;

	double err = N - (double)n;
	if (err > 0.5)
		n++;

	return n;
}
