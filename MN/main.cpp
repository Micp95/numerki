//Defines - select your exercise

//#define INTER_LAGRANGE		//Lagrange interpolation
//#define INTER_HERMITE			//Hermite interpolation
//#define APPROX_LS				//Least-squares function approximation
//#define INTEGRAL				//Integral calculator (Methods: Rectangle, Trapezoidal, Simpson)
//#define DIFFERENTIAL			//Diffrential - Methods: Euler and Rugeg-Kutt (II, IV)
//#define EQUATIONS
#define NONLINEAR


//Standard includes
#include <iostream>
#include <string>
using namespace std;

//Includes to execises
#include "Point.h"
#include "polynomial.h" 

#include "interpolationLagrange.h"
#include "interpolationHermite.h"
#include "approximation.h"
#include "Integral.h"
#include "differential.h"
#include "equations.h"
#include "derivative.h"
#include "NonlinearControl.h"


//Math functions

//INTEGRAL
double funGauss(double x) {
	return 1 / x;
}

//DIFFERENTIAL
double funEuler(double x, double y) {
	return x*x + y;
}

//NONLINEAR 
double funx(double x) {//function
	return x*x -2 ;
}

double funP(double x){//derivative of the function
	return 2 * x;
}



//Main function
int main() {

#ifdef INTER_LAGRANGE
	double* x = new double[3]{ 1,2,4 };
	double* y = new double[3]{ 1,2,1 };
	int stopien = 2;

	interpolationLagrange wiel(x, y, 3, stopien);
	for (int k = -10; k < 10; k++)
		cout << k << "\t" << wiel.PointValue(k) << endl;


	cout << "Error\t" << wiel.SSE() << endl;


	double* wsp = wiel.getFactors();
	cout << "Wspolczynniki:";
	for (int k = 0; k <= stopien; k++)
		cout << "\t" << wsp[k] ;
	cout << endl;


	delete[] x;
	delete[] y;
#endif

#ifdef INTER_HERMITE
	double fx1[] = { 3,-1 };
	double fx2[] = { 1,-5,-8 };
	//float fx3[] = { 2,8,56 }; 

	int ilosc = 2;

	HNode x1(2, 0, fx1);
	HNode x2(3, 1, fx2);
	//HNode x3(3, 1, fx3);

	HNode* tmp = new HNode[ilosc];
	tmp[0] = x1;
	tmp[1] = x2;
	//tmp[2] = x3;

	interpolationHermite inter(ilosc, tmp);
	for (int k = -10; k < 10; k++)
		cout <<k <<"\t"<< inter.PointValue(k) << endl;
	
	cout <<"Error::\t"<<inter.SSE()<<endl;
	//cout << inter.PointValue(0) << endl;
	//cout << inter.PointValue(1) << endl;
	//cout << inter.PointValue(2) << endl;


	double* wsp = inter.getFactors();
	cout << "Wspolczynniki:";
	for (int k = 0; k <= ilosc; k++)
		cout << "\t" << wsp[k];
	cout << endl;


	delete[] tmp;
#endif

#ifdef APPROX_LS

	Point tab[] = { Point(1,-1), Point(3,101), Point(5,739),
		Point(6,1499), Point(7,2729)};

	int ilosc = 5;
	int stopien = 4 ;

	Approximation aproks(tab, ilosc, stopien, Approximation::sel::Gauss);

	for (int k = -10; k < 10; k++)
		cout << k << "\t" << aproks.PointValue(k) << endl;


	cout << "Error::\t" << aproks.SSE() << endl;


	double* wsp = aproks.getOutput();
	cout << "Wspolczynniki:";
	for (int k = 0; k <= stopien; k++)
		cout << "\t" << wsp[k];
	cout << endl;


	delete[] wsp;

#endif

#ifdef INTEGRAL

	Integral integral(funGauss, Integral::method::Rectangle);
	double a = 1, b = 2;

	cout.precision(30);
	cout << "Pole liczone metoda prostokatow wynosi:\t" << integral.calculate(a, b, 1000)<<endl;
	integral.setMethod(Integral::method::Trapezoidal);
	cout << "Pole liczone metoda trapezow wynosi:\t" << integral.calculate(a, b, 1000) << endl;
	integral.setMethod(Integral::method::Simpson);
	cout << "Pole liczone metoda Simpsona wynosi:\t" << integral.calculate(a, b	, 1000) << endl;

	integral.setMethod(Integral::method::Gauss);


	cout << "\n\n";
	for (int k = 2; k <=5;k++)
		cout << "Pole liczone metoda Gaussa (wezly "<<k<<") wynosi:\t" << integral.calculate(a, b, k) << endl;


#endif
	
#ifdef DIFFERENTIAL

	differential rozniczka(funEuler,differential::Euler);

	cout.precision(30);
	cout << "Rozniczka policzona metoda Eulera:\n"
		<<"\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 1) << endl << endl;

	cout << "Rozniczka policzona metoda Eulera:\n"
		<<"\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\tw punkcie 0.3\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 0.3) << endl << endl;


	rozniczka.setMethod(differential::Heuna);
	cout << "Rozniczka policzona metoda Heuna:\n"
		<< "\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\tw punkcie 0.3\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 0.3) << endl <<endl;

	rozniczka.setMethod(differential::modEulera);
	cout << "Rozniczka policzona metoda modyfikowana Eulera:\n"
		<< "\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\tw punkcie 0.3\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 0.3) << endl << endl;


	rozniczka.setMethod(differential::RungegKuttyII);
	cout << "Rozniczka policzona metoda RK II:\n"
		<< "\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\tw punkcie 0.3\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 0.3) << endl << endl;

	rozniczka.setMethod(differential::RungegKuttyIV);
	cout << "Rozniczka policzona metoda RK IV:\n"
		<< "\tPrzedzial 0-1\tkrok: 0.1\tstart point (0,0.1)\tw punkcie 0.3\n\t\t";
	cout << rozniczka.calculateDiffrential(0, 1, Point(0, 0.1), 0.1, 0.3) << endl << endl;


#endif

#ifdef EQUATIONS
	int wymiar = 3;
	double matrix[3][3] = { 3,1,1,
							0,5,1,
							1,1,6};

	double vector[3] = {5,6,8};

	Equations equations(matrix[0], vector, wymiar);


	double* res = equations.result(Equations::Jacob);

	cout << "Obliczone metoda Jacobiego:" << endl;
	for (int k = 0; k < wymiar; k++)
		cout << "\t" << res[k];
	cout << endl;

	delete[] res;


	cout << "Obliczone metoda Gaussa:" << endl;
	res = equations.result(Equations::Gauss);
	for (int k = 0; k < wymiar; k++)
		cout << "\t" << res[k];
	cout << endl;

	delete[] res;

#endif

#ifdef NONLINEAR
	double x;
	double err = 1E-6;

	double A = 1;
	double B = 4;

	cout << "Error:" << err << endl<<endl;

	NonlinearControl ag(A, B, funx, funP, err);


	cout << "Obliczone metoda Newton:" << endl;
	x = ag.count(NonlinearControl::select::Newton);
	cout <<"f( " << x << " )=\t" << funx(x) << endl;
	cout << endl;


	cout << "Obliczone metoda siecznych:" << endl;
	x = ag.count(NonlinearControl::select::Secant);
	cout << "f( " << x << " )=\t" << funx(x) << endl;
#endif

	system("pause");
	return 0;
}



