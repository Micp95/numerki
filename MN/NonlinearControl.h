#ifndef NONLINEARCONTROL_H
#define NONLINEARCONTROL_H



class NonlinearControl {

public:
	enum select {
		Newton, Secant
	};

	NonlinearControl(double a, double b, double(*funPointer)(double x));
	NonlinearControl(double a, double b, double(*funPointer)(double x), double(*funDerPointer)(double x),double err);
	NonlinearControl(double a, double b, double* pol, int n);
	double count(select sel);
	void changeInterval(double a, double b);

private:
	int a, b;
	double(*funPointer)(double x);
	double(*funDerPointer)(double x);


	double* pol;
	int n;
	double epsilon;

	double methodNewton();
	double methodSecant();

	double modul(double x) {
		if (x < 0)
			return -x;
		else
			return x;
	}

	double  w(double x) //algorytm Hornera - obliczanie wartosci wielomianu
	{

		double acummult = 1, acumsum = 0;

		for (int k = 0; k <n; k++) {

			acummult = 1;
			for (int p = 0; p < k; p++)
				acummult *= x;

			acumsum += pol[k] * acummult;
		}
		return acumsum;
	}

};




#endif