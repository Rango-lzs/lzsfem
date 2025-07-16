#pragma once
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix2d.h"
#include <assert.h>
#include "Eigen/Dense"

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class Matrix3d;	// general 3D Matrixrix of doubles
class Matrix3ds;	// symmetric 3D Matrixrix of doubles
class Matrix3da;	// anti-symmetric 3D Matrixrix of doubles
class Matrix3dd;	// diagonal Matrixrix of doubles
class Matrix2d;


using StressTensor = Matrix3ds;
using StrainTensor = Matrix3ds;

//-----------------------------------------------------------------------------
//! This class describes a diagonal Matrixrix of doubles in 3D

class Matrix3dd
{
public:
	// default constructor
	Matrix3dd(){}

	// constructors
	explicit Matrix3dd(double a);
	Matrix3dd(double a0, double a1, double a2);

	// assignment operators
	Matrix3dd& operator = (const Matrix3dd& m);
	Matrix3dd& operator = (double a);

	// access operators
	double operator () (int i, int j) const;
	double& diag(int i);
	const double& diag(int i) const;

	// arithmetic operators
	Matrix3dd operator + (const Matrix3dd& m) const;
	Matrix3dd operator - (const Matrix3dd& m) const;
	Matrix3dd operator * (const Matrix3dd& m) const;
	Matrix3dd operator * (double a) const;
	Matrix3dd operator / (double a) const;

	Matrix3dd operator - () const;

	// arithmetic operators for Matrix3ds
	Matrix3ds operator + (const Matrix3ds& m) const;
	Matrix3ds operator - (const Matrix3ds& m) const;
	Matrix3ds operator * (const Matrix3ds& m) const;

	// arithmetic operators for Matrix3d
	Matrix3d operator + (const Matrix3d& m) const;
	Matrix3d operator - (const Matrix3d& m) const;
	Matrix3d operator * (const Matrix3d& m) const;

	// arithmetic operators for Matrix3da const;
	Matrix3d operator + (const Matrix3da& m) const;
	Matrix3d operator - (const Matrix3da& m) const;
	Matrix3d operator * (const Matrix3da& m) const;

	// arithmetic assignment operators
	Matrix3dd& operator += (const Matrix3dd& m);
	Matrix3dd& operator -= (const Matrix3dd& m);
	Matrix3dd& operator *= (const Matrix3dd& m);
	Matrix3dd& operator *= (double a);
	Matrix3dd& operator /= (double a);

	// Matrixrix-vector multiplication
	Vector3d operator * (const Vector3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

	double xx() const { return d[0]; }
	double yy() const { return d[1]; }
	double zz() const { return d[2]; }

	// TODO: Make these constexpr
	double xy() const { return 0.0; }
	double yz() const { return 0.0; }
	double xz() const { return 0.0; }

protected:
	double	d[3];	// the diagonal elements

	friend class Matrix3d;
	friend class Matrix3ds;
	friend class Matrix3da;
};

inline Matrix3dd operator * (double a, const Matrix3dd& d) { return d*a; }

//-----------------------------------------------------------------------------
//! This class describes a symmetric 3D Matrixrix of doubles

class Matrix3ds
{
protected:
	// This enumeration can be used to remember the order
	// in which the components are stored.
	enum {
		XX = 0,
		XY = 1,
		YY = 2,
		XZ = 3,
		YZ = 4,
		ZZ = 5 };
public:
	// default constructor
	Matrix3ds(){}

	// constructors
	explicit Matrix3ds(double a);
	Matrix3ds(double xx, double yy, double zz, double xy, double yz, double xz);
	Matrix3ds(const Matrix3dd& d);
	Matrix3ds(const Matrix3ds& d);

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;

	double& xx() { return m[XX]; }
	double& yy() { return m[YY]; }
	double& zz() { return m[ZZ]; }
	double& xy() { return m[XY]; }
	double& yz() { return m[YZ]; }
	double& xz() { return m[XZ]; }

	const double& xx() const { return m[XX]; }
	const double& yy() const { return m[YY]; }
	const double& zz() const { return m[ZZ]; }
	const double& xy() const { return m[XY]; }
	const double& yz() const { return m[YZ]; }
	const double& xz() const { return m[XZ]; }

	// arithmetic operators for Matrix3dd objects
	Matrix3ds operator + (const Matrix3dd& d) const;
	Matrix3ds operator - (const Matrix3dd& d) const;
	Matrix3ds operator * (const Matrix3dd& d) const;

	// arithmetic operators
	Matrix3ds operator + (const Matrix3ds& t) const;
	Matrix3ds operator - (const Matrix3ds& t) const;
	Matrix3d  operator * (const Matrix3ds& t) const;
	Matrix3ds operator * (double g) const;
	Matrix3ds operator / (double g) const;

	// arithmetic operators for Matrix3d objects
	Matrix3d operator + (const Matrix3d& t) const;
	Matrix3d operator - (const Matrix3d& t) const;
	Matrix3d operator * (const Matrix3d& t) const;
	
    // arithmetic operators for Matrix3da objects
    Matrix3d operator + (const Matrix3da& t) const;
    Matrix3d operator - (const Matrix3da& t) const;
    Matrix3d operator * (const Matrix3da& t) const;
    
	// unary operators
	Matrix3ds operator - () const;
	
	// arithmetic assignment operators
	Matrix3ds& operator += (const Matrix3ds& t);
	Matrix3ds& operator -= (const Matrix3ds& t);
	Matrix3ds& operator *= (const Matrix3ds& t);
	Matrix3ds& operator *= (double g);
	Matrix3ds& operator /= (double g);

	// arithmetic assignment operators for Matrix3dd
	Matrix3ds& operator += (const Matrix3dd& d);
	Matrix3ds& operator -= (const Matrix3dd& d);

	// Matrixrix-vector multiplication
	Vector3d operator * (const Vector3d& r) const;

	// trace
	double tr() const;

	// determinant
	double det() const;

	// intialize to zero
	void zero();

	// initialize to unit tensor
	void unit();

	// deviator
	Matrix3ds dev() const;

	// isotropic part
	Matrix3ds iso() const;

	// return the square 
	Matrix3ds sqr() const;

	// calculates the inverse
	Matrix3ds inverse() const;
    double invert(Matrix3ds& Ai);
	
	// determine eigen values and vectors
	void eigen(double d[3], Vector3d r[3] = 0) const;
	void exact_eigen(double l[3]) const;
	void eigen2(double d[3], Vector3d r[3] = 0) const;

	// L2-norm 
	double norm() const;

	// double contraction
	double dotdot(const Matrix3ds& S) const;

	// "effective" or von-Mises norm
	double effective_norm() const;

	// the "max shear" value
	double max_shear() const;

protected:
	double m[6];	// stores data in the order xx, xy, yy, xz, yz, zz

	friend class Matrix3dd;
	friend class Matrix3d;
};

inline Matrix3ds operator * (double a, const Matrix3ds& m) { return m*a; }

//-----------------------------------------------------------------------------
//! This class describes an anti-symmetric 3D Matrixrix of doubles
//! The Matrixrix is defined such that for a vector b the following is true:
//! A.b = a x b where A = Matrix3da(a).
//!
//     | 0 -z  y |   |   0  d0  d2 |
// A = | z  0 -x | = | -d0   0  d1 |
//     |-y  x  0 |   | -d2 -d1   0 |
//

class Matrix3da
{
public:
	// default constructor
	Matrix3da(){}

	// constructors
	Matrix3da(double xy, double yz, double xz);

	// calculates the antisymmetric Matrixrix from a vector
	// A.b = a x b where A = Matrix3da(a).
	explicit Matrix3da(const Vector3d& a);

	// access operator
	double operator () (int i, int j) const;

	double& xy() { return d[0]; }
	double& yz() { return d[1]; }
	double& xz() { return d[2]; }

	const double& xy() const { return d[0]; }
	const double& yz() const { return d[1]; }
	const double& xz() const { return d[2]; }

	Matrix3da operator + (const Matrix3da& a);
	Matrix3da operator - (const Matrix3da& a);

	Matrix3da operator - () const;

	Matrix3da operator * (double g) const;

	Matrix3da transpose() const;

	// Matrixrix algebra
	Matrix3d operator * (const Matrix3d& a);

	// return the equivalent vector
	Vector3d vec() const { return Vector3d(-d[1], d[2], -d[0]); }

	Vector3d operator * (const Vector3d& a);

    // arithmetic operators for Matrix3ds objects
    Matrix3d operator + (const Matrix3ds& t) const;
    Matrix3d operator - (const Matrix3ds& t) const;
    
protected:
	double	d[3];	// stores xy, yz, xz

	friend class Matrix3dd;
	friend class Matrix3ds;
	friend class Matrix3d;
};


//-----------------------------------------------------------------------------
//! This class describes a general 3D Matrixrix of doubles
class Matrix3d
{
public:
	// default constructor
	Matrix3d() {}

	explicit Matrix3d(double a);

	// constructors
	Matrix3d(double a00, double a01, double a02,
		  double a10, double a11, double a12,
		  double a20, double a21, double a22);

	Matrix3d(double m[3][3]);
	Matrix3d(double a[9]);

	Matrix3d(const Matrix3dd& m);
	Matrix3d(const Matrix3ds& m);
	Matrix3d(const Matrix3da& m);

	Matrix3d(const Matrix2d& m);

	Matrix3d(const Vector3d& e1, const Vector3d& e2, const Vector3d& e3);

	// assignment operators
	Matrix3d& operator = (const Matrix3dd& m);
	Matrix3d& operator = (const Matrix3ds& m);
	Matrix3d& operator = (const Matrix3d& m);
	Matrix3d& operator = (const double m[3][3]);

	// Matrix3d
	Matrix3d operator - () 
	{
		return Matrix3d(-d[0][0], -d[0][1], -d[0][2], \
					 -d[1][0], -d[1][1], -d[1][2], \
					 -d[2][0], -d[2][1], -d[2][2]);
	}

	// access operators
	double& operator () (int i, int j);
	const double& operator () (int i, int j) const;
	double* operator [] (int i);
	const double* operator [] (int i) const;

	// arithmetic operators
	Matrix3d operator + (const Matrix3d& m) const;
	Matrix3d operator - (const Matrix3d& m) const;
	Matrix3d operator * (const Matrix3d& m) const;
	Matrix3d operator * (double a) const;
	Matrix3d operator / (double a) const;

	// arithmetic operators for Matrix3dd
	Matrix3d operator + (const Matrix3dd& m) const;
	Matrix3d operator - (const Matrix3dd& m) const;
	Matrix3d operator * (const Matrix3dd& m) const;

	// arithmetic operators for Matrix3ds
	Matrix3d operator + (const Matrix3ds& m) const;
	Matrix3d operator - (const Matrix3ds& m) const;
	Matrix3d operator * (const Matrix3ds& m) const;

	// arithmetic assignment operators
	Matrix3d& operator += (const Matrix3d& m);
	Matrix3d& operator -= (const Matrix3d& m);
	Matrix3d& operator *= (const Matrix3d& m);
	Matrix3d& operator *= (double a);
	Matrix3d& operator /= (double a);

	// arithmetic assignment operators for Matrix3dd
	Matrix3d& operator += (const Matrix3dd& m);
	Matrix3d& operator -= (const Matrix3dd& m);
	Matrix3d& operator *= (const Matrix3dd& m);

	// arithmetic assignment operators for Matrix3ds
	Matrix3d& operator += (const Matrix3ds& m);
	Matrix3d& operator -= (const Matrix3ds& m);
	Matrix3d& operator *= (const Matrix3ds& m);

	// Matrixrix-vector muliplication
	Vector3d operator * (const Vector3d& r) const;

	// determinant
	double det() const;

	// trace
	double trace() const;

	// zero the Matrixrix
	void zero();

	// make unit Matrixrix
	void unit();

	// return a column vector from the Matrixrix
	Vector3d col(int j) const;

	// return a row vector from the Matrixrix
	Vector3d row(int j) const;

	// set the column of the Matrixrix
	void setCol(int i, const Vector3d& a);

	// set the row of the Matrixrix
	void setRow(int i, const Vector3d& a);

	// return the symmetric Matrixrix 0.5*(A+A^T)
	Matrix3ds sym() const;

	// return the antisymmetric Matrixrix 0.5*(A-A^T)
	Matrix3da skew() const;

	// calculates the inverse
	Matrix3d inverse() const;
    
	// inverts the Matrixrix.
	bool invert();

	// calculates the transpose
	Matrix3d transpose() const;

	// calculates the transposed inverse
	Matrix3d transinv() const;

	// calculate the skew-symmetric Matrixrix from a vector
	void skew(const Vector3d& v);

	// calculate the one-norm
	double norm() const;

	// double contraction
	double dotdot(const Matrix3d& T) const;

	// polar decomposition
	void right_polar(Matrix3d& R, Matrix3ds& U) const;
	void left_polar(Matrix3ds& V, Matrix3d& R) const;

	// return identity Matrixrix
	static Matrix3d identity() { return Matrix3d(1,0,0, 0,1,0, 0,0,1); }

protected:
	//double d[3][3];	// Matrixrix data
    Eigen::Matrix3d m_matrix;

	friend class Matrix3dd;
	friend class Matrix3ds;
	friend class Matrix3da;
};

// outer product for vectors
inline Matrix3d operator & (const Vector3d& a, const Vector3d& b)
{
	return Matrix3d(a.x*b.x, a.x*b.y, a.x*b.z,
				 a.y*b.x, a.y*b.y, a.y*b.z,
				 a.z*b.x, a.z*b.y, a.z*b.z);
}

inline Matrix3ds dyad(const Vector3d& a)
{
	return Matrix3ds(a.x*a.x, a.y*a.y, a.z*a.z, a.x*a.y, a.y*a.z, a.x*a.z);
}

// c_ij = a_i*b_j + a_j*b_i
inline Matrix3ds dyads(const Vector3d& a, const Vector3d& b)
{
	return Matrix3ds(2.0*a.x*b.x, 2.0*a.y*b.y, 2.0*a.z*b.z, a.x*b.y + a.y*b.x, a.y*b.z + a.z*b.y, a.x*b.z + a.z*b.x);
}

// skew-symmetric Matrixrix of dual vector
inline Matrix3d skew(const Vector3d& a)
{
    return Matrix3d(   0, -a.z,  a.y,
                  a.z,    0, -a.x,
                 -a.y,  a.x,    0);
}

//-----------------------------------------------------------------------------
// This class stores a 2nd order diagonal tensor
class Matrix3fd
{
public:
	Matrix3fd() { x = y = z = 0.f; }
	Matrix3fd(float X, float Y, float Z) { x = X; y = Y; z = Z; }

public:
	float x, y, z;
};

//-----------------------------------------------------------------------------
// Matrix3fs stores a 2nd order symmetric tensor
//
class Matrix3fs
{
public:
	// constructors
	Matrix3fs() { x = y = z = xy = yz = xz = 0; }
	Matrix3fs(float fx, float fy, float fz, float fxy, float fyz, float fxz)
	{
		x = fx; y = fy; z = fz;
		xy = fxy; yz = fyz; xz = fxz;
	}

	// operators
	Matrix3fs& operator += (const Matrix3fs& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		xy += v.xy;
		yz += v.yz;
		xz += v.xz;

		return (*this);
	}

	// operators
	Matrix3fs& operator -= (const Matrix3fs& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		xy -= v.xy;
		yz -= v.yz;
		xz -= v.xz;

		return (*this);
	}

	Matrix3fs& operator *= (float g)
	{
		x *= g;
		y *= g;
		z *= g;
		xy *= g;
		yz *= g;
		xz *= g;

		return (*this);
	}

	Matrix3fs& operator /= (float g)
	{
		x /= g;
		y /= g;
		z /= g;
		xy /= g;
		yz /= g;
		xz /= g;

		return (*this);
	}

	Matrix3fs operator + (const Matrix3fs& a) { return Matrix3fs(x + a.x, y + a.y, z + a.z, xy + a.xy, yz + a.yz, xz + a.xz); }
	Matrix3fs operator - (const Matrix3fs& a) { return Matrix3fs(x - a.x, y - a.y, z - a.z, xy - a.xy, yz - a.yz, xz - a.xz); }

	Matrix3fs operator * (float a)
	{
		return Matrix3fs(x * a, y * a, z * a, xy * a, yz * a, xz * a);
	}

	Matrix3fs operator / (float g)
	{
		return Matrix3fs(x / g, y / g, z / g, xy / g, yz / g, xz / g);
	}

	Vector3f operator*(Vector3f& r)
	{
        return Vector3f(
			x * r.x + xy * r.y + xz * r.z,
			xy * r.x + y * r.y + yz * r.z,
			xz * r.x + yz * r.y + z * r.z);
	}

	// Effective or von-mises value
	float von_mises() const
	{
		float vm;
		vm = x * x + y * y + z * z;
		vm -= x * y + y * z + x * z;
		vm += 3 * (xy * xy + yz * yz + xz * xz);
		vm = (float)sqrt(vm >= 0.0 ? vm : 0.0);
		return vm;
	}

	// principle values
	void Principals(float e[3]) const;

	// principle directions
	Vector3f PrincDirection(int l);

	// deviatroric principle values
	void DeviatoricPrincipals(float e[3]) const;

	// max-shear value
	float MaxShear() const;

	// eigen-vectors and values
    void eigen(Vector3f e[3], float l[3]) const;

	// trace
	float tr() const { return x + y + z; }

	// determinant
	float det() const { return (x * y * z + xy * yz * xz + xz * xy * yz - y * xz * xz - x * yz * yz - z * xy * xy); }

	// L2 norm
	float norm() const {
		double d = x * x + y * y + z * z + 2 * (xy * xy + yz * yz + xz * xz);
		return (float)sqrt(d);
	}

public:
	float x, y, z;
	float xy, yz, xz;
};

double fractional_anisotropy(const Matrix3fs& m);

///////////////////////////////////////////////////////////////////
// Matrix3f

class Matrix3f
{
public:
	Matrix3f() { zero(); }

	Matrix3f(float a00, float a01, float a02, float a10, float a11, float a12, float a20, float a21, float a22)
	{
		d[0][0] = a00; d[0][1] = a01; d[0][2] = a02;
		d[1][0] = a10; d[1][1] = a11; d[1][2] = a12;
		d[2][0] = a20; d[2][1] = a21; d[2][2] = a22;
	}

	Matrix3f(const Matrix3fs& a)
	{
		d[0][0] = a.x; d[0][1] = a.xy; d[0][2] = a.xz;
		d[1][0] = a.xy; d[1][1] = a.y; d[1][2] = a.yz;
		d[2][0] = a.xz; d[2][1] = a.yz; d[2][2] = a.z;
	}

	float* operator [] (int i) { return d[i]; }
	float& operator () (int i, int j) { return d[i][j]; }
	float operator () (int i, int j) const { return d[i][j]; }

	Matrix3f operator * (Matrix3f& m)
	{
		Matrix3f a;

		int k;
		for (k = 0; k < 3; k++)
		{
			a[0][0] += d[0][k] * m[k][0]; a[0][1] += d[0][k] * m[k][1]; a[0][2] += d[0][k] * m[k][2];
			a[1][0] += d[1][k] * m[k][0]; a[1][1] += d[1][k] * m[k][1]; a[1][2] += d[1][k] * m[k][2];
			a[2][0] += d[2][k] * m[k][0]; a[2][1] += d[2][k] * m[k][1]; a[2][2] += d[2][k] * m[k][2];
		}

		return a;
	}

	Vector3f operator*(const Vector3f& a) const
	{
        return Vector3f(
			d[0][0] * a.x + d[0][1] * a.y + d[0][2] * a.z,
			d[1][0] * a.x + d[1][1] * a.y + d[1][2] * a.z,
			d[2][0] * a.x + d[2][1] * a.y + d[2][2] * a.z
			);
	}

	Matrix3f& operator *= (float g)
	{
		d[0][0] *= g;	d[0][1] *= g; d[0][2] *= g;
		d[1][0] *= g;	d[1][1] *= g; d[1][2] *= g;
		d[2][0] *= g;	d[2][1] *= g; d[2][2] *= g;
		return (*this);
	}

	Matrix3f& operator /= (float g)
	{
		d[0][0] /= g;	d[0][1] /= g; d[0][2] /= g;
		d[1][0] /= g;	d[1][1] /= g; d[1][2] /= g;
		d[2][0] /= g;	d[2][1] /= g; d[2][2] /= g;
		return (*this);
	}

	Matrix3f operator += (const Matrix3f& a)
	{
		d[0][0] += a.d[0][0]; d[0][1] += a.d[0][1]; d[0][2] += a.d[0][2];
		d[1][0] += a.d[1][0]; d[1][1] += a.d[1][1]; d[1][2] += a.d[1][2];
		d[2][0] += a.d[2][0]; d[2][1] += a.d[2][1]; d[2][2] += a.d[2][2];
		return (*this);
	}

	Matrix3f operator -= (const Matrix3f& a)
	{
		d[0][0] -= a.d[0][0]; d[0][1] -= a.d[0][1]; d[0][2] -= a.d[0][2];
		d[1][0] -= a.d[1][0]; d[1][1] -= a.d[1][1]; d[1][2] -= a.d[1][2];
		d[2][0] -= a.d[2][0]; d[2][1] -= a.d[2][1]; d[2][2] -= a.d[2][2];
		return (*this);
	}

	Matrix3fs sym() const
	{
		return Matrix3fs(d[0][0], d[1][1], d[2][2], 0.5f * (d[0][1] + d[1][0]), 0.5f * (d[1][2] + d[2][1]), 0.5f * (d[0][2] + d[2][0]));
	}

	void zero()
	{
		d[0][0] = d[0][1] = d[0][2] = 0.f;
		d[1][0] = d[1][1] = d[1][2] = 0.f;
		d[2][0] = d[2][1] = d[2][2] = 0.f;
	}

	Vector3f col(int i) const
	{
        Vector3f r;
		switch (i)
		{
		case 0: r.x = d[0][0]; r.y = d[1][0]; r.z = d[2][0]; break;
		case 1: r.x = d[0][1]; r.y = d[1][1]; r.z = d[2][1]; break;
		case 2: r.x = d[0][2]; r.y = d[1][2]; r.z = d[2][2]; break;
		}
		return r;
	}

	Vector3f row(int i) const
	{
        Vector3f r;
		switch (i)
		{
		case 0: r.x = d[0][0]; r.y = d[0][1]; r.z = d[0][2]; break;
		case 1: r.x = d[1][0]; r.y = d[1][1]; r.z = d[1][2]; break;
		case 2: r.x = d[2][0]; r.y = d[2][1]; r.z = d[2][2]; break;
		}
		return r;
	}

	Matrix3f transpose() const
	{
		return Matrix3f(
			d[0][0], d[1][0], d[2][0],
			d[0][1], d[1][1], d[2][1],
			d[0][2], d[1][2], d[2][2]
		);
	}

	// inverts the Matrixrix.
	bool invert();

public:
	float d[3][3];
};

// The following file contains the actual definition of the class functions
#include "Matrix3d.hpp"
