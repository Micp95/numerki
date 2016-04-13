#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H


class polynomial
{
public:
	virtual double PointValue(double x) = 0;
	virtual double SSE() = 0;
private:
};

#endif