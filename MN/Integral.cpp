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
	case Integral::Gauss:
		return funGauss(a, b, --n);
	case Integral::Newton:
		return funNewton(a, b, n);
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

double Integral::funGauss(double a, double b, int n) {
	double* t = new double[n + 1];
	double* Ft = new double[n + 1];

	GaussTab myGaussTab(n);

	auto tFun = [&a, &b](double x) -> double {				//funckja anonimowa - lambda
		return (a+b)/2 + ((b-a)/2)*x;
	};

	for (int k = 0; k <= n; k++) {
		t[k] = tFun(myGaussTab.getXk(k));
		Ft[k] = PointValue(t[k]);
	}

	double Q = 0;
	for (int k = 0; k <= n; k++)
		Q += myGaussTab.getAk(k)*Ft[k];

	Q *= (b - a) / 2;

	return Q;
}


double Integral::funNewton(double a, double b, int n) {
	double res = 0;





	return res;
}


Integral::GaussTab::GaussTab(int N)
{
	N++;
	this->N = N;
	k = new int[N];
	Xk = new double[N];
	Ak = new double[N];

	for (int p = 0; p < N; p++)
		k[p] = p;

	if (N == 2) {
		Xk[0] = 0.57735;
		Xk[1] = -0.57735;

		Ak[0] = Ak[1] = 1;
	}
	else if (N == 3) {
		Xk[0] =  -0.774597;
		Xk[1] = 0;
		Xk[2] = 0.774597;

		Ak[0] = Ak[2] = 5.0/9.0;
		Ak[1] = 8.0 / 9.0;
	}
	else if (N == 4) {
		Xk[0] = -0.861136;
		Xk[1] = -0.339981;
		Xk[2] = 0.339981;
		Xk[3] = 0.861136;

		Ak[0] = Ak[3] = 0.347855;
		Ak[1] = Ak[2] = 0.652145;
	}
	else if (N == 5) {
		Xk[0] = -0.906180;
		Xk[1] = -0.538469;
		Xk[2] = 0;
		Xk[3] = 0.538469;
		Xk[4] = 0.906180;

		Ak[0] = Ak[4] = 0.236927;
		Ak[1] = Ak[3] = 0.478629;
		Ak[2] = 0.568889;
	}

}
Integral::GaussTab::~GaussTab()
{
	delete[] k;
	delete[] Xk;
	delete[] Ak;
}

int Integral::GaussTab::getSize() {
	return N;
}
int Integral::GaussTab::getK(int k) {
	if (k < N)
		return this->k[k];
	else
		return 0;
}
double Integral::GaussTab::getXk(int k) {
	if (k < N)
		return this->Xk[k];
	else
		return 0;
}
double Integral::GaussTab::getAk(int k) {
	if (k < N)
		return this->Ak[k];
	else
		return 0;
}