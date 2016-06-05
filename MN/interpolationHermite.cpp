#include "interpolationHermite.h"
#include "global.h"


HNode::HNode(int n, double x, double * var):n(n),x(x){
	this->var = new double[n];

	for (int k = 0; k < n; k++)
		this->var[k] = var[k];

}

HNode::HNode(const HNode& nx){
	x = nx.x;
	n = nx.n;
	var = new double[n];
	for (int k = 0; k < n; k++)
		this->var[k] = nx.var[k];
}

HNode::~HNode(){
	if (var != NULL)
		delete[] var;
}

HNode & HNode::operator=(const HNode & nx){
	if (var != NULL)
		delete[] var;
	x = nx.x;
	n = nx.n;
	var = new double[n];
	for (int k = 0; k < n; k++)
		this->var[k] = nx.var[k];

	return *this;
}


interpolationHermite::interpolationHermite(int n, HNode * node):nn(n){
	MyNodes = new HNode[n];
	for (int k = 0; k < n; k++)
		MyNodes[k] = node[k];

	initialization();
	interpolation();
}

interpolationHermite::~interpolationHermite(){
	if (a != NULL)
		delete[] a;

	if (xi != NULL)
		delete[] xi;

	if (MyNodes != NULL)
		delete[] MyNodes;

	if ( nodeIndex != NULL)
		delete[] nodeIndex;
}

void interpolationHermite::initialization() {
	int acum = 0;
	for (int k = 0; k < nn; k++) {
		acum += MyNodes[k].n;
	}
	n = acum;

	a = new double[n];
	xi = new double[n];
	nodeIndex = new int[n];

	//uzupelnienie tablicy xi - powielona informacjia o x w zaleznosci od jego opisu
	for (int k = 0, p =0, s; k < nn; k++)
		for (s = 0; s < MyNodes[k].n; s++) {
			nodeIndex[p] = k;
			xi[p++] = MyNodes[k].x;
		}

}

double interpolationHermite::thisSameFun(HNode x, int degree) {
	return x.var[degree] / factorial(degree);
}

double interpolationHermite::otherFun(double x1, double x2, double fx1, double fx2) {
	return (fx2 - fx1) / (x2 - x1);
}

double interpolationHermite::factorial(int x) {
	if (x == 0)
		return 1;
	return x * factorial(x - 1);
}


void interpolationHermite::interpolation(){
	double* zi = new double[n];
	int act = n, count = 0;


	//zerowa iteracja - wypelnienie tablicy wartosiacmi f(xi)
	for (int k = 0, p = 0, s; k < nn; k++)
		for (s = 0; s < MyNodes[k].n; s++)
			zi[p++] = MyNodes[k].var[0];

	do {
		//przechowanie wyniku (elementy przekatnej)
		a[count++] = zi[0];

		
		//pomocnicze wypisanie tablicy - kontrola
		//for (int k = 0; k < act; k++)
		//	cout << zi[k] << " ";
		//cout << endl;
		

		//wykonanie algorytmu dla tablicy
		for (int k = 1; k < act; k++) {

			if (zi[k] == zi[k - 1]) 
				zi[k - 1] = thisSameFun(MyNodes[nodeIndex[k-1+(count -1)]], count);
			else
				zi[k - 1] = otherFun(xi[k + (count - 1)], xi[k - 1], zi[k], zi[k - 1]);
		}

		act--;//zmniejszenie tablicy
	} while (act > 0);//dopoki tablica nie jest "pusta"


	delete[] zi;
}


double interpolationHermite::PointValue(double x) {
	double acum = 0,tmp;
	for (int k = 0,p; k < n; k++) {
		tmp = a[k];
		for (p = 1; p <= k; p++)
			tmp *= (x - xi[p - 1]);
		acum += tmp;
	}
	return acum;
}


double interpolationHermite::SSE() {
	double acum = 0, y2;
	for (int k = 0; k < nn; k++) {
		y2 = PointValue(MyNodes[k].x);
		acum += (MyNodes[k].var[0] - y2)*(MyNodes[k].var[0] - y2);
	}
	return acum / n;;
}

