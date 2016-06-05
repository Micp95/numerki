#include "equations.h"


Equations::Matrix::Matrix() {
	size = 0;
	tab = 0;
}
Equations::Matrix::Matrix(int size) :size(size) {
	tab = new double*[size];
	for (int k = 0; k < size; k++)
		tab[k] = new double[size];

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			set(0.0, x, y);
}
Equations::Matrix::Matrix(Matrix& other) {
	size = other.size;
	tab = new double*[size];
	for (int k = 0; k < size; k++)
		tab[k] = new double[size];

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			tab[x][y] = other.get(x, y);
}
Equations::Matrix::~Matrix() {
	for (int k = 0; k < size; k++)
		delete[] tab[k];
	delete[] tab;
}
void Equations::Matrix::set(double var, int x, int y) {
	if (x < size && x >= 0 && y < size && y >= 0)
		tab[x][y] = var;
}
double Equations::Matrix::get(int x, int y) {
	if (x < size && x >= 0 && y < size && y >= 0)
		return tab[x][y];
	return 0.0f;
}

Equations::Matrix & Equations::Matrix::operator+(Matrix & other){
	Matrix* res = new Matrix(other.size);

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			res->set(get(x, y) + other.get(x, y), x, y);

	return *res;
}

Equations::Matrix & Equations::Matrix::operator*(Matrix & other){
	Matrix* res = new Matrix(other.size);


	for (int i = 0; i<other.size; i++){
		for (int j = 0; j<other.size; j++){
			res->set(0,i,j);
			for (int k = 0; k<other.size; k++){
				res->set(res->get(i,j)+(get(i,k)*other.get(k,j)), i, j);
			}
		}
	}
	return *res;
}

Equations::Matrix & Equations::Matrix::operator*(double k){
	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++) {
			if (get(x, y)*k != -0.0)
				set(get(x, y)*k, x, y);
		}
	return *this;
}

Equations::Matrix & Equations::Matrix::operator=(Matrix & other){
	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			set(other.get(x, y), x, y);

	return *this;
}

void Equations::Matrix::invers(){
	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++) {
			if ( get(x,y) != 0.0)
				set(1.0 / get(x, y), x, y);
		}
}


Equations::~Equations() {
	if (data != 0)
		delete data;
	if (vectr != 0)
		delete vectr;
}
Equations::Equations(double * matrix, double* b, int size):size(size){
	data = new Matrix(size);

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			data->set(matrix[x+size*y], y, x);

	vectr = new double[size];
	for (int k = 0; k < size; k++)
		vectr[k] = b[k];

}

double * Equations::result(sel method)
{
	if (method == sel::Gauss)
		return GaussFun();
	else if (method == sel::Jacob)
		return JacobFun();

	return nullptr;
}


double * Equations::GaussFun(){

	Matrix* tmpMatrix = new Matrix(size);

	for (int x = 0; x < size; x++)
		for (int y = 0; y < size; y++)
			tmpMatrix->set(data->get(x, y), y,x);

	double* vec = new double[size];
	for (int x = 0; x < size; x++)
		vec[x] = vectr[x];

	double tmp, tmp2;
	double* output = new double[size];


	for (int k = 0; k < size; k++) {
		for (int p = k + 1; p < size; p++) {

			tmp = tmpMatrix->get(k, p) / tmpMatrix->get(k, k);

			for (int r = k; r < size; r++) {
				tmp2 = tmp * tmpMatrix->get(r, k);
				tmpMatrix->set(tmpMatrix->get(r, p) - tmp2, r, p);
			}


			tmp2 = tmp * vec[k];
			vec[p] = vec[p] - tmp2;

		}
	}

	
	for (int k = size - 1; k >= 0; k--) {
		tmp = vec[k];
		for (int p = size - 1; p != k; p--)
			tmp -= (tmpMatrix->get(p, k) *output[p]);
		tmp /= tmpMatrix->get(k, k);
		output[k] = tmp;
	}
	delete tmpMatrix;
	delete[] vec;


	return output;
}



double * Equations::JacobFun(){
	Matrix* L = new Matrix(size);
	Matrix* D = new Matrix(size);
	Matrix* U = new Matrix(size);

	for (int x = 0; x < size; x++) {
		for (int y = 0; y < size; y++) {
			if (x < y) 
				L->set(data->get(x,y), x, y);
			else if (x > y)
				U->set(data->get(x, y), x, y);
			else 
				D->set(data->get(x, y), x, y);
		}
	}

	Matrix* LU = new Matrix(*L + *U);
	Matrix* Dinv = new Matrix(*D);
	Dinv->invers();

	Matrix* vectorMatrix = new Matrix(size);
	for (int y = 0; y < size; y++)
		vectorMatrix->set(vectr[y], y, 0);

	Matrix* xMatrix = new Matrix(size);

	for (int k = 1; k < size+10; k++) 
		*xMatrix = (*Dinv)*(*LU)*(*xMatrix)*(-1.0) + (*Dinv)*(*vectorMatrix);
	
	delete L;
	delete D;
	delete U;
	delete LU;
	delete Dinv;
	delete vectorMatrix;

	double* outbput = new double[size];
	for (int x = 0; x < size; x++)
		outbput[x] = xMatrix->get(x, 0);

	delete xMatrix;

	return outbput;
}
