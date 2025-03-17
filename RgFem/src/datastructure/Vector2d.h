#pragma once
#include <math.h>
#include "femcore/fem_export.h"

class Vector2d
{
public:
	// constructor
	Vector2d() { r[0] = r[1] = 0; }
	explicit Vector2d(double v) { r[0] = r[1] = v; }
	Vector2d(double x, double y) { r[0] = x; r[1] = y; }

	bool operator == (const Vector2d& r) const { return (r[0] == r.r[0]) && (r[1] == r.r[1]); }

	// access operators
	double operator [] (int i) const { return r[i]; }
	double& operator [] (int i) { return r[i]; }

	double& x() { return r[0]; }
	double& y() { return r[1]; }

	double x() const { return r[0]; }
	double y() const { return r[1]; }

	double norm() const { return sqrt(r[0] * r[0] + r[1] * r[1]); }
	double norm2() const { return r[0] * r[0] + r[1] * r[1]; }

	double unit() { double R = norm(); if (R != 0) { r[0] /= R; r[1] /= R; }; return R; }

public: // arithmetic operators

	Vector2d operator + (const Vector2d& v) const { return Vector2d(r[0]+v.r[0], r[1]+v.r[1]); }
	Vector2d operator - (const Vector2d& v) const { return Vector2d(r[0]-v.r[0], r[1]-v.r[1]); }
	Vector2d operator * (double g) const { return Vector2d(r[0]*g, r[1]*g); }
	Vector2d operator / (double g) const { return Vector2d(r[0]/g, r[1]/g); }

	Vector2d& operator += (const Vector2d& v) { r[0] += v.r[0]; r[1] += v.r[1]; return *this; }
	Vector2d& operator -= (const Vector2d& v) { r[0] -= v.r[0]; r[1] -= v.r[1]; return *this; }
	Vector2d& operator *= (double g) { r[0] *= g; r[1] *= g; return *this; }
	Vector2d& operator /= (double g) { r[0] /= g; r[1] /= g; return *this; }

    Vector2d operator - () const { return Vector2d(-r[0], -r[1]); }
    
	// dot product
	double operator * (const Vector2d& v) const { return r[0]*v[0] + r[1]*v[1]; }

public:
	double r[2];
};

//-----------------------------------------------------------------------------
class vec2i
{
public:
	vec2i() { x = y = 0; }
	vec2i(int X, int Y) { x = X; y = Y; }

    bool operator == (const vec2i& r) const { return (x == r.x) && (y == r.y); }

public:
	int		x, y;
};

//-----------------------------------------------------------------------------
class vec2f
{
public:
	vec2f() { x = y = 0.f; }
	vec2f(float rx, float ry) { x = rx; y = ry; }

public:
	float	x, y;
};
