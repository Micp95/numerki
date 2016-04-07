#ifndef APPROXIMATION_H
#define APPROXIMATION_H




class Approximation{
public:

	enum sel {
		Cramer,Gausse
	};
	//Pomocnicze struktury
	struct Point{
		double x, y;
		Point(double x, double y) :x(x), y(y) {}
		Point(Point& other) {
			x = other.x;
			y = other.y;
		}
		Point() {}
		Point& operator=(Point& other) {
			x = other.x;
			y = other.y;
			return *this;
		}

	};

	class Matrix{
	public:
		Matrix(int size);
		Matrix(Matrix& other);
		~Matrix();
		void set(double var, int x, int y);
		double get(int x, int y);
		double getDeterminant();
	private:
		int size;
		double ** tab;
	};


	Approximation(Point* input, int size, int degree,sel mySel);
	~Approximation();

	double PointValue(double x);
	double* getOutput() { return output; }
	double SSE();

private:
	sel mySel;

	void approximation();
	double matrixSum(int k);
	double vectorSum(int k);

	void CramerFun();
	void GausseFun();

	int sizeIn;
	Point* input;

	int funDe;
	Matrix* tmpMatrix;
	double* vectr;

	double* output;


};

#endif