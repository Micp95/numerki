#include "interpolation.h"
#include "global.h"



float* polynomial::Interpolation(float* x, float* y, int n, int degree) {
	float* b = new float[degree];

	if (degree > n)
		return NULL;
	float acummult, acumsum;

	for (int k = 0, i, j; k < degree; k++) {//kolejne wyrazy

		acumsum = 0;
		for (i = 0; i <= k; i++) {//Suma

			acummult = 1;
			for (j = 0; j <= k; j++)//Mnozenie
				if (i != j)
					acummult *= (x[i] - x[j]);

			acumsum += y[i] / acummult;
		}
		b[k] = acumsum;
	}

	return b;
}


float polynomial::PointValue(float x) {

	float acummult = 1, acumsum = 0;

	for (int k = 0; k < degree; k++) {

		acummult = 1;
		for (int p = 0; p < k; p++)
			acummult *= (x - xn[p]);

		acumsum += bn[k] * acummult;
	}
	return acumsum;
}

float polynomial::SSE() {
	float acum = 0, y2;
	for (int k = 0; k < n; k++) {
		y2 = PointValue(xn[k]);
		acum += (yn[k] - y2)*(yn[k] - y2);
	}
	return acum;
}
