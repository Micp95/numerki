#ifndef EQUATIONS_H
#define EQUATIONS_H


class Equations {
public:
 enum sel {
		Jacob, Gauss
	};
private:

	class Matrix {
	public:
		Matrix();
		Matrix(int size);
		Matrix(Matrix& other);
		~Matrix();
		void set(double var, int x, int y);
		double get(int x, int y);
		Matrix& operator+(Matrix& other);
		Matrix& operator*(Matrix& other);
		Matrix& operator*(double  k);
		Matrix& operator=(Matrix& other);
		void invers();
	private:
		int size;
		double ** tab;
	};

public:
	Equations(double * matrix,double* b, int size);
	double* result(sel method);
	~Equations();

private :
	double*  GaussFun();
	double*  JacobFun();


	Matrix* data;
	double* vectr;

	int size;
};


#endif