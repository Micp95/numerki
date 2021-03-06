#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "polynomial.h"

class interpolationLagrange: public polynomial {
	double* bn;
	double* xn;
	double* yn;
	int n;

	int degree;

	double* Interpolation(double* x, double* y, int n, int degree);
public:
	interpolationLagrange(double* x, double* y, int n, int degree) : xn(x), yn(y), degree(++degree), n(n) {

		bn = Interpolation(xn, yn, n, degree);
	}
	~interpolationLagrange();

	virtual double PointValue(double x);
	virtual double SSE();

	double* getFactors() { return bn; }
};


#endif


/*

Zadanie:

uk�ad r�wna� - a0 + a1x + a2x^2 = f(x)
potrzebne 3 punkty do rozwiazania



interplonacja Lagrangea
idea: szukanie wielomianu ktory w 1 punkcie przyjmie odpowiednia wartosc, w reszcie zero
rozwiazaniem jest suma tych wielomianow

wzor: na stronie - ogolno dosepny


Li - wielomian ktory przyjmuje kakretna wartosc w jednym punkcie - reszta 0

wzor interpolacyjny to wzor na:
Ln(x) = suma (f(x) * Li(x) )

liczenie wspolczynnikow - praktyczny wzro
wykozystanie postaci Newtona

Wn(x) = suma(bk*pk(x))
k : 0- > n
rozpisany:
p0(x) = 0
p1(x) = x-x0
p2(x) - (x-x0)(x-x1)
.
.
.

potrzebne bk - wz�r interpolacji z pdf

1. uzytkownik podaje ilosc punktow
2. wczytanie punktow i wartosci

3. obliczenie wspolczynnikow b
4. wartosci wielomianow w dowolnym punkcie
5. wartosc bledu interpolacji

wzor na blad:

analiza:
wyrysowanie punktow - odgadniecie wzoru
kolejno wybor stopnia wielomianu
suma kwadratow z odleglosci pomiedzy punktem podanym, a punktem wyliczonym (jesli ilosc punktow mniejsza niz


bn - pn

Wx = bo + b1(x - x0)+ b2(x-x0)(x - x1)




*/