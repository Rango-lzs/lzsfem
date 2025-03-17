/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <math.h>
#include "datastructure/Vector2d"

class Vector3d
{
public:
	// constructors
	Vector3d() : x(0), y(0), z(0) {}
	explicit Vector3d(double a) : x(a), y(a), z(a) {}
	Vector3d(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
	Vector3d(const Vector2d& v) { x = v.r[0]; y = v.r[1]; z = 0.0; }

	// operators
	Vector3d operator + (const Vector3d& r) const { return Vector3d(x+r.x, y+r.y, z+r.z); }
	Vector3d operator - (const Vector3d& r) const { return Vector3d(x-r.x, y-r.y, z-r.z); }

	Vector3d operator * (double a) const { return Vector3d(x*a, y*a, z*a); }
	Vector3d operator / (double a) const { return Vector3d(x/a, y/a, z/a); }

	Vector3d& operator += (const Vector3d& r) { x += r.x; y += r.y; z += r.z; return (*this); }
	Vector3d& operator -= (const Vector3d& r) { x -= r.x; y -= r.y; z -= r.z; return (*this); }

	Vector3d& operator *= (double a) { x*=a; y*=a; z*=a; return (*this); }
	Vector3d& operator /= (double a) { x/=a; y/=a; z/=a; return (*this); }

	Vector3d operator - () const { return Vector3d(-x, -y, -z); }

    double& operator() (int i)
    {
        switch(i)
        {
            case 0: {return x; break;}
            case 1: {return y; break;}
            case 2: {return z; break;}
            default: {return x; break;}
        }
    }
    
    double operator() (int i) const
    {
        switch(i)
        {
            case 0: {return x; break;}
            case 1: {return y; break;}
            case 2: {return z; break;}
            default: {return x; break;}
        }
    }
    
	// dot product
	double operator * (const Vector3d& r) const { return (x*r.x + y*r.y + z*r.z); }

	// cross product
	Vector3d operator ^ (const Vector3d& r) const { return Vector3d(y*r.z-z*r.y,z*r.x-x*r.z,x*r.y-y*r.x); }

	// normalize the vector
	double unit()
	{
		double d = sqrt(x*x+y*y+z*z);
		if (d != 0) { x/=d; y/=d; z/=d; }
		return d;
	}

	// return a normalized version of this vector
	Vector3d normalized() const { 
		double d = sqrt(x*x + y*y + z*z); 
		d = (d == 0.0? d = 1.0 : d = 1.0/d); 
		return Vector3d(x*d, y*d, z*d);
	}

	// length of vector
	double norm() const { return sqrt(x*x+y*y+z*z); }

	// length square of vector
	double norm2() const { return (x*x + y*y + z*z); }

public:
	// NOTE: Added to simplify integration with FEBio Studio
	Vector3d Normalize() { unit(); return *this; }
	Vector3d Normalized() const { Vector3d v(x, y, z); v.unit(); return v; }
	double Length() const { return norm(); }
	double SqrLength() const { return norm2(); }
	bool operator == (const Vector3d& a) const { return ((a.x == x) && (a.y == y) && (a.z == z)); }

 public:
	double x, y, z;
};


//-----------------------------------------------------------------------------
class Vector3f
{
public:
	Vector3f() { x = y = z = 0; }
	Vector3f(float rx, float ry, float rz) { x = rx; y = ry; z = rz; }

	Vector3f operator + (const Vector3f& v) const { return Vector3f(x + v.x, y + v.y, z + v.z); }
	Vector3f operator - (const Vector3f& v) const { return Vector3f(x - v.x, y - v.y, z - v.z); }
	Vector3f operator ^ (const Vector3f& v) const
	{
		return Vector3f(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
	}

	float operator * (const Vector3f& v) const { return (x * v.x + y * v.y + z * v.z); }

	Vector3f operator * (const float g) const { return Vector3f(x * g, y * g, z * g); }
	Vector3f operator / (const float g) const { return Vector3f(x / g, y / g, z / g); }

	const Vector3f& operator += (const Vector3f& v) { x += v.x; y += v.y; z += v.z; return (*this); }
	const Vector3f& operator -= (const Vector3f& v) { x -= v.x; y -= v.y; z -= v.z; return (*this); }
	const Vector3f& operator /= (const float& f) { x /= f; y /= f; z /= f; return (*this); }
	const Vector3f& operator /= (const int& n) { x /= n; y /= n; z /= n; return (*this); }
	const Vector3f& operator *= (const float& f) { x *= f; y *= f; z *= f; return (*this); }

	Vector3f operator - () { return Vector3f(-x, -y, -z); }

	float Length() const { return (float)sqrt(x * x + y * y + z * z); }

	float SqrLength() const { return (float)(x * x + y * y + z * z); }

	Vector3f& Normalize()
	{
		float L = Length();
		if (L != 0) { x /= L; y /= L; z /= L; }

		return (*this);
	}

public:
	float x, y, z;
};

inline Vector3d to_Vector3d(const Vector3f& r) { return Vector3d((double)r.x, (double)r.y, (double)r.z); }
inline Vector3f to_Vector3f(const Vector3d& r) { return Vector3f((float)r.x, (float)r.y, (float)r.z); }
