#ifndef INTERPOLATION_H
#define INTERPOLATION_H


class polynomial {
	float* bn;
	float* xn;
	float* yn;
	int n;

	int degree;

	float* Interpolation(float* x, float* y, int n, int degree);
public:
	polynomial(float* x, float* y, int n, int degree) : xn(x), yn(y), degree(degree), n(n) {
		bn = Interpolation(xn, yn, n, degree);
	}
	float PointValue(float x);
	float SSE();
};




#endif


/*

Zadanie:

uk³ad równañ - a0 + a1x + a2x^2 = f(x)
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

potrzebne bk - wzór interpolacji z pdf

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