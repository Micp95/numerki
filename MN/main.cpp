//#define INTER_HERMITE			//Hermite interpolation
//#define INTER_LAGRANGE		//Lagrange interpolation
//#define APPROX_LS				//Least-squares function approximation
#define INTEGRAL				//Integral calculator (Methods: Rectangle, Trapezoidal, Simpson)
//#define ALL_LIB


#include <iostream>
#include <string>
using namespace std;


//Includes to execises
#if defined(INTER_LAGRANGE) || defined(ALL_LIB)
#include "interpolation.h"
#endif
#if defined(INTER_HERMITE) || defined(ALL_LIB)
#include "hermite.h"
#endif
#if defined(APPROX_LS) || defined(ALL_LIB)
#include "approximation.h"
#endif
#if defined(INTEGRAL) || defined(ALL_LIB)
#include "Integral.h"
#include "polynomial.h"
#endif



//Others definitions
#if defined(INTEGRAL) || defined(ALL_LIB)
double fun(double x) {
	return x*x+2*x+3;
}
double funGauss(double x) {
	return 1 / x;
}
#endif

//Main function
int main() {

#ifdef INTER_LAGRANGE
	float* x = new float[3]{ 1,2,4 };
	float* y = new float[3]{ 1,2,1 };


	polynomial wiel(x, y, 3, 3);
	for (int k = 0; k < 5; k++)
		cout << k << "\t" << wiel.PointValue(k) << endl;

	cout << "Error\t" << wiel.SSE() << endl;

	delete[] x;
	delete[] y;
#endif

#ifdef INTER_HERMITE
	float fx1[] = { 3,-1 };
	float fx2[] = { 1,-5,-8 };
	//float fx3[] = { 2,8,56 };

	HNode x1(2, 0, fx1);
	HNode x2(3, 1, fx2);
	//HNode x3(3, 1, fx3);

	HNode* tmp = new HNode[2];
	tmp[0] = x1;
	tmp[1] = x2;
	//tmp[2] = x3;

	Hermite inter(2, tmp);
	cout << inter.PointValue(0) << endl;
	cout << inter.PointValue(1) << endl;
	//cout << inter.PointValue(2) << endl;

	delete[] tmp;
#endif

#ifdef APPROX_LS

	Approximation::Point tab[] = { Approximation::Point(1,-1), Approximation::Point(3,101), Approximation::Point(5,739),
		Approximation::Point(6,1499), Approximation::Point(7,2729)};
	/*
	Approximation::Point tab[] = {
		Approximation::Point(1,62),Approximation::Point(2,232),
		Approximation::Point(3,1330),Approximation::Point(4,5984),
		Approximation::Point(5, 20590), Approximation::Point(6, 57952),
		Approximation::Point(7,140642),Approximation::Point(8,305080),
		Approximation::Point(9,606334),Approximation::Point(10,1123640),
		Approximation::Point(11,1966642),Approximation::Point(12,3282352),
		Approximation::Point(13,5262830),Approximation::Point(14,8153584),
	};
	

	Approximation aproks(tab, 14, 6, Approximation::sel::Gauss);
	
	Approximation::Point tab[] = { Approximation::Point(1,1), Approximation::Point(3,2), Approximation::Point(4,4),
		Approximation::Point(6,4), Approximation::Point(8,5), Approximation::Point(9,7),
		Approximation::Point(11,8), Approximation::Point(14,9) };
	*/
	Approximation aproks(tab, 5, 5, Approximation::sel::Gauss);

	double* out = aproks.getOutput();
	
	cout << endl << endl;
	for (int k = 0; k < 6; k++)
		cout << out[k] << " "; 
	
	cout << endl<<endl;

	for (int k = 1; k < 15; k++)
		cout << k << ".\t" << aproks.PointValue(k) << endl;
	
	cout <<endl <<"Error:\t"<< aproks.SSE() << endl;

	delete[] out;

#endif

#ifdef INTEGRAL

	Integral integral(funGauss, Integral::method::Rectangle);
	double a = 1, b = 2;

	cout << "Pole liczone metoda prostokatow wynosi:\t" << integral.calculate(a, b, 1000)<<endl;
	integral.setMethod(Integral::method::Trapezoidal);
	cout << "Pole liczone metoda trapezow wynosi:\t" << integral.calculate(a, b, 1000) << endl;
	integral.setMethod(Integral::method::Simpson);
	cout << "Pole liczone metoda Simpsona wynosi:\t" << integral.calculate(a, b	, 1000) << endl;

	integral.setMethod(Integral::method::Gauss);
	cout << "Pole liczone metoda Gaussa (stopien 1) wynosi:\t" << integral.calculate(a, b, 1) << endl;
	cout << "Pole liczone metoda Gaussa (stopien 3) wynosi:\t" << integral.calculate(a, b, 3) << endl;


#endif
	



#ifdef ALL_LIB




#endif


	system("pause");
	return 0;
}


