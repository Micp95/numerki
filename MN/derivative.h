#ifndef DERIVATIVE_H
#define DERIVATIVE_H


class Derivative {

public:
	Derivative(double* pol, int n);
	~Derivative();

	double dervativeIn(double x,int degree);//Algorytm Show-Trauba

	double wartosc(double x, int degree) {
		double tmp;
		for (int i = 0; i <= n; i++) { b[i] = a[i]; }
		for (int j = 1; j <= degree; j++)
		{
			tmp = wspol(0,x);
			for (int i = 0; i <= n; i++) { b[i] = c[i]; }
		}
		return wart(0, x);
	}


private:
	double *a, *b, *c;
	int n, p, q;

	double s(int j);
	double r(int j);
	double T(int i, int j,double x);
	double duplicate(double x, int k);



	void generatePQ();



	float wart(int k,double x)            //wartosc wielmianu
	{
		if (k == n)
		{
			return b[n];
		}
		else
		{
			return wart(k + 1,x)*x + b[k];
		}
	}

	float wspol(int k, double x)           //wspolczynniki wielmianu po podzieleniu
	{
		if (k == n) { c[k - 1] = b[k]; return b[k]; }
		else { if (k>0) { c[k - 1] = wspol(k + 1,x)*x + b[k]; } return wspol(k + 1,x)*x + b[k]; }
	}

};


#endif