#ifndef  HERMITE_H
#define HERMITE_H

#include "global.h"

//Struktura przechowujaca wejscie do algorytmu
struct HNode{
	HNode(int n,float x,float*var);
	HNode():var(NULL) {};
	HNode(const HNode& nh);
	~HNode();

	HNode& operator= (const HNode& x);

	int n;			//ilosc informacji
	float x;		//wartosci x
	float* var;		//informacje o x: 0-f(x), n-pochodne
};


class Hermite{
public:
	Hermite(int n, HNode* node);
	~Hermite();

	float PointValue(float x);
private:
	HNode* MyNodes;

	int n,nn;

	float* a;
	float* xi;
	int*  nodeIndex;

	void interpolation();
	void initialization();

	float thisSameFun(HNode x,int degree);
	float otherFun(float x1, float x2, float fx1, float fx2);

	float factorial(int x);
};


#endif