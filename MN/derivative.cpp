#include "derivative.h"

Derivative::Derivative(double* pol, int n):n(n){
	a = new double[n];
	b = new double[n];
	c = new double[n];
	for (int k = 0; k < n; k++)
		a[k] = pol[k];
	generatePQ();
	n--;
}
Derivative::~Derivative() {
	delete[] a;
	delete[] b;
	delete[] c;
}

double Derivative::dervativeIn(double x, int degree) {
	if (x == 0) 
		return T(n, degree, x);
	else
		return T(n, degree, x) / duplicate(x, degree%q);
}

double Derivative::T(int i, int j,double x) {

	if (x == 0)                     //by mozna bylo obliczyc pochodna w punkcie x=0 (?)
		return a[j];
	else if (j == -1) 
		return a[n - i - 1] * duplicate(x, s(i + 1));
	else if (j == i)
		return a[n] * duplicate(x, s(0));
	else 
		return T(i - 1, j - 1, x) + T(i - 1, j,x)*duplicate(x, r(i - j));
}

double Derivative::s(int j) {
	return (n - j) % q;
}
double Derivative::r(int j) {
	if (j%q == 0)
		return q;
	else
		return 0;
}


double Derivative::duplicate(double x, int k) {
	double acum = 1;
	for (int p = 0; p < k; p++)
		acum *= x;
	
	return acum;
}

void Derivative::generatePQ(){
	q = 0;

	do {
		if (q > n + 1)
			break;//error - throw exception!
		q++;
		p = (n + 1) / q;
	} while (p*q != n + 1);
}
