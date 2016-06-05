#ifndef  HERMITE_H
#define HERMITE_H

#include "global.h"
#include "polynomial.h"

//Struktura przechowujaca wejscie do algorytmu
struct HNode{
	HNode(int n,double x, double*var);
	HNode():var(NULL) {};
	HNode(const HNode& nh);

	HNode& operator= (const HNode& x);

	~HNode();

	int n;			//ilosc informacji
	double x;		//wartosci x
	double* var;		//informacje o x: 0-f(x), n-pochodne
};


class interpolationHermite: public polynomial{
public:
	interpolationHermite(int n, HNode* node);
	~interpolationHermite();

	virtual double PointValue(double x);
	virtual double SSE();

	double* getFactors() { return a; }
private:
	HNode* MyNodes;

	int n,nn;

	double* a;
	double* xi;
	int*  nodeIndex;

	void interpolation();
	void initialization();

	double thisSameFun(HNode x,int degree);
	double otherFun(double x1, double x2, double fx1, double fx2);

	double factorial(int x);
};


#endif