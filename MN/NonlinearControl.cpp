#include "NonlinearControl.h"
#include "Point.h"
#include "derivative.h"
#include <cmath>

NonlinearControl::NonlinearControl(double a, double b, double(*funPointer)(double x)) {
	this->a = a;
	this->b = b;
	this->funPointer = funPointer;
}

NonlinearControl::NonlinearControl(double a, double b, double(*funPointer)(double x), double(*funDerPointer)(double x), double err) {
	this->a = a;
	this->b = b;
	this->funPointer = funPointer;
	this->funDerPointer = funDerPointer;
	this->epsilon = err;
}

NonlinearControl::NonlinearControl(double a, double b, double* pol, int n) {
	this->a = a;
	this->b = b;
	this->pol = pol;
	this->n = n;

}
double NonlinearControl::count(select sel) {
	if (sel == Newton)
		return methodNewton();
	else
		return methodSecant();
}
void NonlinearControl::changeInterval(double a, double b) {
	this->a = a;
	this->b = b;
}


double NonlinearControl::methodNewton() {

	double xk = a;
	double aa, bb;

	double err = -4E140;

	 while(true){
		xk = xk - funPointer( xk) / funDerPointer(xk);

		//if (modul(funPointer(xk) - funPointer(pre)) <= epsilon)
		if (modul(funPointer(xk)) <= epsilon)
			break;
	}

	return xk;
}
double NonlinearControl::methodSecant() {

	double xk1 = a;
	double xk2 = b;

	double xk = b;


	while(true){
		xk = xk2 - funPointer(xk2) * ((xk2 - xk1)/(funPointer(xk2) - funPointer(xk1)));
		xk1 = xk2;
		xk2 = xk;

		//if (modul(funPointer(xk) - funPointer(pre)) <= epsilon)
		if (modul(funPointer(xk)) <= epsilon)
			break;
	}

	return xk;
}