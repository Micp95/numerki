#ifndef POINT_H
#define POINT_H

struct Point{
	Point() {};
	Point(double x, double y) : x(x), y(y) {};

	double x, y;
};


#endif