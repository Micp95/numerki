#ifndef APPROXIMATION_H
#define APPROXIMATION_H

#include "polynomial.h"
#include "Point.h"

class Approximation: public polynomial{
public:

	enum sel {
		Cramer,Gauss
	};

	//Pomocnicze struktury
	class Matrix{
	public:
		Matrix(int size);
		Matrix(Matrix& other);
		~Matrix();
		void set(double var, int x, int y);
		double get(int x, int y);
		double getDeterminant();
	private:
		int size;
		double ** tab;
	};


	Approximation(Point* input, int size, int degree,sel mySel);
	~Approximation();

	virtual double PointValue(double x);
	double* getOutput() { return output; }
	virtual double SSE();

private:
	sel mySel;

	void approximation();
	double matrixSum(int k);
	double vectorSum(int k);

	void CramerFun();
	void GaussFun();

	int sizeIn;
	Point* input;

	int funDe;
	Matrix* tmpMatrix;
	double* vectr;

	double* output;

};

#endif