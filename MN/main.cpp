//#define INTER_HERMITE		//Hermite interpolation
//#define INTER_LAGRANGE		//Lagrange interpolation
#define APPROX_LS			//Least-squares function approximation


#include <iostream>
#include <string>


#ifdef INTER_LAGRANGE
#include "interpolation.h"
#endif
#ifdef INTER_HERMITE
#include "hermite.h"
#endif
#ifdef APPROX_LS
#include "approximation.h"
#endif


using namespace std;

int main() {

#ifdef INTER_LAGRANGE
	float* x = new float[3]{ 1,2,4 };
	float* y = new float[3]{ 1,2,1 };


	polynomial wiel(x, y, 3, 3);
	for (int k = 0; k < 5; k++)
		cout << k << "\t" << wiel.PointValue(k) << endl;

	cout << "Error\t" << wiel.SSE() << endl;
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
#endif

#ifdef APPROX_LS
	/*
	Approximation::Point tab[] = { Approximation::Point(1,1), Approximation::Point(3,2), Approximation::Point(4,4),
		Approximation::Point(6,4), Approximation::Point(8,5), Approximation::Point(9,7),
		Approximation::Point(11,8), Approximation::Point(14,9) };

	Approximation::Point tab[] = { Approximation::Point(1,-1), Approximation::Point(3,101), Approximation::Point(5,739),
		Approximation::Point(6,1499), Approximation::Point(7,2739)};
	*/
	
	Approximation::Point tab[] = {
		Approximation::Point(1,62),Approximation::Point(2,232),
		Approximation::Point(3,1330),Approximation::Point(4,5984),
		Approximation::Point(5, 20590), Approximation::Point(6, 57952),
		Approximation::Point(7,140642),Approximation::Point(8,305080),
		Approximation::Point(9,606334),Approximation::Point(10,1123640),
		Approximation::Point(11,1966642),Approximation::Point(12,3282352),
		Approximation::Point(13,5262830),Approximation::Point(14,8153584),
	};

	Approximation aproks(tab, 14, 7, Approximation::sel::Gausse);
	
	double* out = aproks.getOutput();
	
	cout << endl << endl;
	for (int k = 0; k < 7; k++)
		cout << out[k] << " "; 
	
	cout << endl<<endl;

	for (int k = 1; k < 15; k++)
		cout << k << ".\t" << aproks.PointValue(k) << endl;
	
	cout <<endl <<"Error:\t"<< aproks.SSE() << endl;

#endif

	system("pause");
	return 0;
}


