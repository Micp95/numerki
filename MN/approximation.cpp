#include "approximation.h"


//Matrix
Approximation::Matrix::Matrix(int size) :size(size) {
	tab = new double*[size];
	for (int k = 0; k < size; k++)
		tab[k] = new double[size];
}

Approximation::Matrix::Matrix(Matrix & other){
	size = other.size;
	tab = new double*[size];
	for (int k = 0; k < size; k++)
		tab[k] = new double[size];

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			tab[x][y] = other.get(x, y);
}

Approximation::Matrix::~Matrix() {
	for (int k = 0; k < size; k++)
		delete[] tab[k];
	delete[] tab;
}

void Approximation::Matrix::set(double var, int x, int y){
	if ( x < size && x >= 0 && y < size && y >= 0)
		tab[x][y] = var;
}

double Approximation::Matrix::get(int x, int y){
	if (x < size && x >= 0 && y < size && y >= 0)
		return tab[x][y];
	return 0.0f;
}

double Approximation::Matrix::getDeterminant(){
	if (size == 1)
		return get(0, 0);
	else if (size == 2) {
		return get(0, 1)*get(1, 0) - get(0, 0)*get(1, 1);
	}
	else
	{
		double acumpl = 0, acummin = 0;
		int actjump;

		double tmp;
		//plusy
		for (int k = 0; k < size; k++) {
			actjump = 0;
			tmp = 1;
			for (int y = k, x = 0; actjump < size; actjump++,x++, y--) {
				if (y < 0)
					y += size;
				tmp *= get(x, y);
			}
			acumpl += tmp;
		}

		//minusy
		for (int k = 0; k < size; k++) {
			actjump = 0;
			tmp = 1;
			for (int y = k, x = 0; actjump < size; actjump++, x++, y++) {
				if (y >= size)
					y = 0;
				tmp *= get(x, y);
			}
			acummin += tmp;
		}

		return acumpl - acummin;
	}
}


//Approximation
Approximation::Approximation(Point * input,int size, int degree, sel mySel):sizeIn(size),funDe(degree+1), mySel(mySel){
	tmpMatrix = new Matrix(funDe);
	this->input = new Point[sizeIn];

	for (int k = 0; k < size; k++)
		this->input[k] = input[k];

	output = new double[funDe];
	vectr = new double[funDe];

	approximation();
}

void Approximation::approximation() {

	//wyliczenie wspolczynnikow macierzy
	for (int x = 0; x < funDe; x++)
		for (int y = 0; y < funDe; y++)
			tmpMatrix->set(matrixSum(x + y), x, y);

	//wyliczenie wektora
	for (int k = 0; k < funDe; k++)
		vectr[k] = vectorSum(k);

	if (mySel == sel::Cramer)
		CramerFun();
	else if (mySel == sel::Gausse)
		GausseFun();

}


void Approximation::CramerFun() {
	//rozwiazanie ukladu - Cramer

	//wyliczenie wyznacznika glownego + stworzenie pomocniczych zmiennych
	double W = tmpMatrix->getDeterminant(), Wx;
	Matrix* newMatrix;

	//petla wyliczajaca wszystkie niewiadome
	for (int k = 0; k < funDe; k++) {

		//Tworzenie macierzy z podstawiona kolumna
		newMatrix = new Matrix(*tmpMatrix);
		for (int p = 0; p < funDe; p++)
			newMatrix->set(vectr[p], k, p);

		//liczenie wyznacznika macierzy pomocniczej
		Wx = newMatrix->getDeterminant();

		//Zapisanie wyniku do tablicy wyjsciowej
		output[k] = Wx / W;

		//Usuniecie zbednej macierzy
		delete newMatrix;
	}
}
#include <iostream>
using namespace std;

void Approximation::GausseFun() {
	double* line = new double[funDe];
	double tmp, tmp2;
	int cont = 0;


	for (int k = 0; k < funDe; k++) {
		for (int p = k+1; p < funDe; p++) {
			tmp = tmpMatrix->get(k, p) / tmpMatrix->get(k, k);
			for (int r = k; r < funDe; r++) {
				tmp2 = tmp * tmpMatrix->get(r, k);
				tmpMatrix->set(tmpMatrix->get(r,p)- tmp2, r, p);
			}
			tmp2 = tmp * vectr[k];
			vectr[p] = vectr[p] - tmp2;
		}

	}

	for (int k = funDe - 1; k >= 0; k--) {
 		tmp = vectr[k];
		for (int p = funDe - 1; p != k; p--)
			tmp -= (tmpMatrix->get(p,k) *output[p]);
		tmp /= tmpMatrix->get(k, k);
		output[k] = tmp;
	}

	delete line;
}

double Approximation::matrixSum(int k) {
	double sum = 0;
	double acum;

	//Suma operacji na calym wejsciu
	for (int m= 0; m < sizeIn; m++) {
		//potegowanie wejscia
		acum = 1;
		for (int p = 0; p < k; p++)
			acum *= input[m].x;
		//Sumowanie
		sum += acum;
	}
	return sum;
}
double Approximation::vectorSum(int k) {
	double sum = 0;
	double acum;

	//Suma operacji na calym wejsciu
	for (int m = 0; m < sizeIn; m++) {
		//potegowanie wejscia
		acum = 1;
		for (int p = 0; p < k; p++) 
			acum *= input[m].x;
		//przemnozenie x^n przez y
		acum *= input[m].y;
		//Sumowanie
		sum += acum;
	}
	return sum;
}

Approximation::~Approximation(){
	delete[] input;
	delete tmpMatrix;
	delete[] vectr;
	delete[] output;
}

double Approximation::PointValue(double x){

	double acummult = 1, acumsum = 0;

	for (int k = 0; k <funDe; k++) {

		acummult = 1;
		for (int p = 0; p < k; p++)
			acummult *= x;

		acumsum += output[k] * acummult;
	}
	return acumsum;
}

double Approximation::SSE() {
	double acum = 0, y2;
	for (int k = 0; k < sizeIn; k++) {
		y2 = PointValue(input[k].x);
		acum += (input[k].y - y2)*(input[k].y - y2);
	}
	return acum/sizeIn;
}
