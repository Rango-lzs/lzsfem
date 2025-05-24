
#pragma once

#include "datastructure/Vector2d.h"

class Matrix2d
{
public:
	// constructors
	Matrix2d(){}
	Matrix2d(double a00, double a01, double a10, double a11)
	{
		d[0][0] = a00; d[0][1] = a01;
		d[1][0] = a10; d[1][1] = a11;
	}

	// access operators
	double& operator () (int i, int j) { return d[i][j]; }
	double operator () (int i, int j) const { return d[i][j]; }
	double* operator [] (int i) { return d[i]; }

public: // arithmetic operations
	Matrix2d operator + (const Matrix2d& m) { return Matrix2d(d[0][0]+m.d[0][0], d[0][1]+m.d[0][1], d[1][0]+m.d[1][0], d[1][1]+m.d[1][1]); }
	Matrix2d operator - (const Matrix2d& m) { return Matrix2d(d[0][0]-m.d[0][0], d[0][1]-m.d[0][1], d[1][0]-m.d[1][0], d[1][1]-m.d[1][1]); }
	Matrix2d operator * (double g) { return Matrix2d(d[0][0]*g, d[0][1]*g, d[1][0]*g, d[1][1]*g); }
	Matrix2d operator / (double g) { return Matrix2d(d[0][0]/g, d[0][1]/g, d[1][0]/g, d[1][1]/g); }

	Matrix2d& operator += (const Matrix2d& m) { d[0][0] += m.d[0][0]; d[0][1] += m.d[0][1]; d[1][0] += m.d[1][0]; d[1][1] += m.d[1][1]; return *this; }
	Matrix2d& operator -= (const Matrix2d& m) { d[0][0] -= m.d[0][0]; d[0][1] -= m.d[0][1]; d[1][0] -= m.d[1][0]; d[1][1] -= m.d[1][1]; return *this; }
	Matrix2d& operator *= (double g) { d[0][0] *= g; d[0][1] *= g; d[1][0] *= g; d[1][1] *= g; return *this; }
	Matrix2d& operator /= (double g) { d[0][0] /= g; d[0][1] /= g; d[1][0] /= g; d[1][1] /= g; return *this; }

	Matrix2d operator * (const Matrix2d& m) { 
		return Matrix2d(
			d[0][0]*m.d[0][0]+d[0][1]*m.d[1][0],
			d[0][0]*m.d[0][1]+d[0][1]*m.d[1][1],
			d[1][0]*m.d[0][0]+d[1][1]*m.d[1][0],
			d[1][0]*m.d[0][1]+d[1][1]*m.d[1][1]); 
	}

public:	// Matrix operations
	Matrix2d inverse() const
	{
		double Di = 1/(d[0][0]*d[1][1] - d[0][1]*d[1][0]);
		return Matrix2d(d[1][1]*Di, -d[0][1]*Di, -d[1][0]*Di, d[0][0]*Di);
	}

	Matrix2d transpose() const
	{
		return Matrix2d(d[0][0], d[1][0], d[0][1], d[1][1]);
	}

	void zero()
	{
		d[0][0] = d[0][1] = d[1][0] = d[1][1] = 0.0;
	}

	void identity()
	{
		d[0][0] = d[1][1] = 1.0;
		d[0][1] = d[1][0] = 0.0;
	}
	
protected:
	double	d[2][2];
};

// Matrix-vector operations
inline Vector2d operator * (Matrix2d& m, Vector2d& a) { return Vector2d(m[0][0]*a[0]+m[0][1]*a[1], m[1][0]*a[0]+m[1][1]*a[1]); }

// dyadic product
inline Matrix2d dyad(Vector2d& a, Vector2d& b) { return Matrix2d(a.r[0]*b.r[0], a.r[0]*b.r[1], a.r[1]*b.r[0], a.r[1]*b.r[1]); }
