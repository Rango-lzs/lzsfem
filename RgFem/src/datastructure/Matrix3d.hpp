

#pragma once
// NOTE: This file is automatically included from Matrix3d.h
// Users should not include this file manually!

//-----------------------------------------------------------------------------
// class Matrix3dd : class describing diagonal 3D matrices of doubles
//-----------------------------------------------------------------------------

// constructor
inline Matrix3dd::Matrix3dd(double a) { d[0] = d[1] = d[2] = a; }
inline Matrix3dd::Matrix3dd(double a0, double a1, double a2) { d[0] = a0; d[1] = a1; d[2] = a2; }

// assignment operators
inline Matrix3dd& Matrix3dd::operator = (const Matrix3dd& m) { d[0] = m.d[0]; d[1] = m.d[1]; d[2] = m.d[2]; return (*this); }
inline Matrix3dd& Matrix3dd::operator = (double a) { d[0] = d[1] = d[2] = a; return (*this); }

// access operators
inline double Matrix3dd::operator () (int i, int j) const { return (i==j? d[i] : 0); }
inline double& Matrix3dd::diag(int i) { return d[i]; }
inline const double& Matrix3dd::diag(int i) const { return d[i]; }

// arithmetic operators
inline Matrix3dd Matrix3dd::operator + (const Matrix3dd& m) const { return Matrix3dd(d[0]+m.d[0],d[1]+m.d[1],d[2]+m.d[2]);}
inline Matrix3dd Matrix3dd::operator - (const Matrix3dd& m) const { return Matrix3dd(d[0]-m.d[0],d[1]-m.d[1],d[2]-m.d[2]);}
inline Matrix3dd Matrix3dd::operator * (const Matrix3dd& m) const { return Matrix3dd(d[0]*m.d[0],d[1]*m.d[1],d[2]*m.d[2]);}
inline Matrix3dd Matrix3dd::operator * (double a) const { return Matrix3dd(d[0]*a, d[1]*a, d[2]*a); }
inline Matrix3dd Matrix3dd::operator / (double a) const { a = 1./a; return Matrix3dd(d[0]*a, d[1]*a, d[2]*a); }
inline Matrix3dd Matrix3dd::operator - () const { return Matrix3dd(-d[0], -d[1], -d[2]); }

// arithmetic operators with Matrix3ds
inline Matrix3ds Matrix3dd::operator + (const Matrix3ds& m) const
{
	return Matrix3ds(
		d[0]+m.m[Matrix3ds::XX],
		d[1]+m.m[Matrix3ds::YY],
		d[2]+m.m[Matrix3ds::ZZ],
		m.m[Matrix3ds::XY],
		m.m[Matrix3ds::YZ],
		m.m[Matrix3ds::XZ]
		);
}

inline Matrix3ds Matrix3dd::operator - (const Matrix3ds& m) const
{
	return Matrix3ds(
		d[0]-m.m[Matrix3ds::XX],
		d[1]-m.m[Matrix3ds::YY],
		d[2]-m.m[Matrix3ds::ZZ],
		-m.m[Matrix3ds::XY],
		-m.m[Matrix3ds::YZ],
		-m.m[Matrix3ds::XZ]
		);
}


inline Matrix3ds Matrix3dd::operator * (const Matrix3ds& m) const
{
	return Matrix3ds(
		d[0]*m.m[Matrix3ds::XX],
		d[1]*m.m[Matrix3ds::YY],
		d[2]*m.m[Matrix3ds::ZZ],
		d[0]*m.m[Matrix3ds::XY],
		d[1]*m.m[Matrix3ds::YZ],
		d[0]*m.m[Matrix3ds::XZ]
	);
}

// arithmetic operators for Matrix3d
inline Matrix3d Matrix3dd::operator + (const Matrix3d& m) const
{
	return Matrix3d(d[0]+m[0][0], m[0][1], m[0][2],
				 m[1][0], d[1]+m[1][1], m[1][2],
				 m[2][0], m[2][1], d[2]+m[2][2]);
}

inline Matrix3d Matrix3dd::operator - (const Matrix3d& m) const
{
	return Matrix3d(d[0]-m[0][0], -m[0][1], -m[0][2],
				 -m[1][0], d[1]-m[1][1], -m[1][2],
				 -m[2][0], -m[2][1], d[2]-m[2][2]);
}

inline Matrix3d Matrix3dd::operator * (const Matrix3d& m) const
{
	return Matrix3d(d[0]*m[0][0], d[0]*m[0][1], d[0]*m[0][2],
				 d[1]*m[1][0], d[1]*m[1][1], d[1]*m[1][2],
				 d[2]*m[2][0], d[2]*m[2][1], d[2]*m[2][2]);
}

// arithmetic operators for Matrix3da
inline Matrix3d Matrix3dd::operator + (const Matrix3da& m) const
{
	return Matrix3d(   d[0],  m.d[0], m.d[2],
				 -m.d[0],    d[1], m.d[1],
				 -m.d[2], -m.d[1],   d[2]);
}

inline Matrix3d Matrix3dd::operator - (const Matrix3da& m) const
{
	return Matrix3d(  d[0],-m.d[0], -m.d[2],
				 m.d[0],   d[1], -m.d[1],
				 m.d[2], m.d[1],   d[2]);
}

inline Matrix3d Matrix3dd::operator *(const Matrix3da& m) const
{
	return Matrix3d(           0,  d[0]*m.d[0], d[0]*m.d[2],
				 -d[1]*m.d[0],            0, d[1]*m.d[1],
				 -d[2]*m.d[2], -d[2]*m.d[1],           0);
}

// arithmetic assignment operators
inline Matrix3dd& Matrix3dd::operator += (const Matrix3dd& m) { d[0] += m.d[0]; d[1] += m.d[1]; d[2] += m.d[2]; return (*this); }
inline Matrix3dd& Matrix3dd::operator -= (const Matrix3dd& m) { d[0] -= m.d[0]; d[1] -= m.d[1]; d[2] -= m.d[2]; return (*this); }
inline Matrix3dd& Matrix3dd::operator *= (const Matrix3dd& m) { d[0] *= m.d[0]; d[1] *= m.d[1]; d[2] *= m.d[2]; return (*this); }
inline Matrix3dd& Matrix3dd::operator *= (double a) { d[0] *= a; d[1] *= a; d[2] *= a; return (*this); }
inline Matrix3dd& Matrix3dd::operator /= (double a) { a = 1./a; d[0] *= a; d[1] *= a; d[2] *= a; return (*this); }

// Matrix-vector multiplication
inline Vector3d Matrix3dd::operator * (const Vector3d& r) const { return Vector3d(r.x*d[0], r.y*d[1], r.z*d[2]); }

// trace
inline double Matrix3dd::tr() const { return d[0]+d[1]+d[2]; }

// determinant
inline double Matrix3dd::det() const { return d[0]*d[1]*d[2]; }

//-----------------------------------------------------------------------------
// class Matrix3ds : this class describes a symmetric 3D Matrix of doubles
//-----------------------------------------------------------------------------

// constructor
inline Matrix3ds::Matrix3ds(double a)
{
	m[XX] = a;
	m[XY] = a;
	m[YY] = a;
	m[XZ] = a;
	m[YZ] = a;
	m[ZZ] = a;
}

inline Matrix3ds::Matrix3ds(double xx, double yy, double zz, double xy, double yz, double xz)
{
	m[XX] = xx;
	m[XY] = xy;
	m[YY] = yy;
	m[XZ] = xz;
	m[YZ] = yz;
	m[ZZ] = zz;
}

inline Matrix3ds::Matrix3ds(const Matrix3dd& d)
{
	m[XX] = d.d[0];
	m[YY] = d.d[1];
	m[ZZ] = d.d[2];
	m[XY] = m[YZ] = m[XZ] = 0.;
}

inline Matrix3ds::Matrix3ds(const Matrix3ds& d)
{
	m[0] = d.m[0];
	m[1] = d.m[1];
	m[2] = d.m[2];
	m[3] = d.m[3];
	m[4] = d.m[4];
	m[5] = d.m[5];
}

// access operator
inline double& Matrix3ds::operator ()(int i, int j)
{
	const int n[] = {0, 1, 3};
	return (i<=j? m[n[j]+i] : m[n[i]+j]);
}

// access operator for const objects
inline const double& Matrix3ds::operator ()(int i, int j) const
{
	const int n[] = {0, 1, 3};
	return (i<=j? m[n[j]+i] : m[n[i]+j]);
}

// operator + for Matrix3dd objects
inline Matrix3ds Matrix3ds::operator + (const Matrix3dd& d) const
{
    return Matrix3ds(m[XX] + d.d[0], m[YY] + d.d[1], m[ZZ] + d.d[2], m[XY], m[YZ], m[XZ]);
}

// operator - for Matrix3dd objects
inline Matrix3ds Matrix3ds::operator - (const Matrix3dd& d) const
{
    return Matrix3ds(m[XX] - d.d[0], m[YY] - d.d[1], m[ZZ] - d.d[2], m[XY], m[YZ], m[XZ]);
}

// operator * for Matrix3dd objects
inline Matrix3ds Matrix3ds::operator * (const Matrix3dd& d) const
{
    return Matrix3ds(m[XX] * d.d[0], m[YY] * d.d[1], m[ZZ] * d.d[2], m[XY] * d.d[1], m[YZ] * d.d[2], m[XZ] * d.d[2]);
}


// operator +
inline Matrix3ds Matrix3ds::operator +(const Matrix3ds& t) const
{
	return Matrix3ds(m[XX]+t.m[XX], m[YY]+t.m[YY], m[ZZ]+t.m[ZZ], m[XY]+t.m[XY], m[YZ]+t.m[YZ],m[XZ]+t.m[XZ]);
}

// operator -
inline Matrix3ds Matrix3ds::operator -(const Matrix3ds& t) const
{
	return Matrix3ds(m[XX]-t.m[XX], m[YY]-t.m[YY], m[ZZ]-t.m[ZZ], m[XY]-t.m[XY], m[YZ]-t.m[YZ],m[XZ]-t.m[XZ]);
}

// operator *
inline Matrix3d Matrix3ds::operator *(const Matrix3ds& t) const
{
	return Matrix3d(
		m[XX] * t.m[XX] + m[XY] * t.m[XY] + m[XZ] * t.m[XZ],
		m[XX] * t.m[XY] + m[XY] * t.m[YY] + m[XZ] * t.m[YZ],
		m[XX] * t.m[XZ] + m[XY] * t.m[YZ] + m[XZ] * t.m[ZZ],

		m[XY] * t.m[XX] + m[YY] * t.m[XY] + m[YZ] * t.m[XZ],
		m[XY] * t.m[XY] + m[YY] * t.m[YY] + m[YZ] * t.m[YZ],
		m[XY] * t.m[XZ] + m[YY] * t.m[YZ] + m[YZ] * t.m[ZZ],

		m[XZ] * t.m[XX] + m[YZ] * t.m[XY] + m[ZZ] * t.m[XZ],
		m[XZ] * t.m[XY] + m[YZ] * t.m[YY] + m[ZZ] * t.m[YZ],
		m[XZ] * t.m[XZ] + m[YZ] * t.m[YZ] + m[ZZ] * t.m[ZZ]
	);
}

// operator *
inline Matrix3ds Matrix3ds::operator * (double g) const
{
	return Matrix3ds(m[XX]*g, m[YY]*g, m[ZZ]*g, m[XY]*g, m[YZ]*g, m[XZ]*g);
}

// operator /
inline Matrix3ds Matrix3ds::operator / (double g) const
{
	g = 1.0/g;
	return Matrix3ds(m[XX]*g, m[YY]*g, m[ZZ]*g, m[XY]*g, m[YZ]*g, m[XZ]*g);
}

// operator + for Matrix3d objects
inline Matrix3d Matrix3ds::operator + (const Matrix3d& d) const
{
	return Matrix3d(m[XX]+d[0][0], m[XY]+d[0][1], m[XZ]+d[0][2],
				 m[XY]+d[1][0], m[YY]+d[1][1], m[YZ]+d[1][2],
				 m[XZ]+d[2][0], m[YZ]+d[2][1], m[ZZ]+d[2][2]);
}

// operator - for Matrix3d objects
inline Matrix3d Matrix3ds::operator - (const Matrix3d& mat3d) const
{
	return Matrix3d(m[XX]-mat3d[0][0], m[XY]-mat3d[0][1], m[XZ]-mat3d[0][2],
				 m[XY]-mat3d[1][0], m[YY]-mat3d[1][1], m[YZ]-mat3d[1][2],
				 m[XZ]-mat3d[2][0], m[YZ]-mat3d[2][1], m[ZZ]-mat3d[2][2]);
}

// operator * for Matrix3d objects
inline Matrix3d Matrix3ds::operator * (const Matrix3d& d) const
{
	return Matrix3d(d[0][0]*m[XX] + d[1][0]*m[XY] + d[2][0]*m[XZ], 
				 d[0][1]*m[XX] + d[1][1]*m[XY] + d[2][1]*m[XZ], 
				 d[0][2]*m[XX] + d[1][2]*m[XY] + d[2][2]*m[XZ],
				 d[0][0]*m[XY] + d[1][0]*m[YY] + d[2][0]*m[YZ], 
				 d[0][1]*m[XY] + d[1][1]*m[YY] + d[2][1]*m[YZ], 
				 d[0][2]*m[XY] + d[1][2]*m[YY] + d[2][2]*m[YZ],
				 d[0][0]*m[XZ] + d[1][0]*m[YZ] + d[2][0]*m[ZZ], 
				 d[0][1]*m[XZ] + d[1][1]*m[YZ] + d[2][1]*m[ZZ], 
				 d[0][2]*m[XZ] + d[1][2]*m[YZ] + d[2][2]*m[ZZ]);
}

// operator + for Matrix3da objects
inline Matrix3d Matrix3ds::operator + (const Matrix3da& d) const
{
    return Matrix3d(m[XX]       , m[XY]+d.xy(), m[XZ]+d.xz(),
                 m[XY]-d.xy(), m[YY]       , m[YZ]+d.yz(),
                 m[XZ]-d.xz(), m[YZ]-d.yz(), m[ZZ]       );
}

// operator - for Matrix3da objects
inline Matrix3d Matrix3ds::operator - (const Matrix3da& d) const
{
    return Matrix3d(m[XX]       , m[XY]-d.xy(), m[XZ]-d.xz(),
                 m[XY]+d.xy(), m[YY]       , m[YZ]-d.yz(),
                 m[XZ]+d.xz(), m[YZ]+d.yz(), m[ZZ]       );
}

// operator * for Matrix3d objects
inline Matrix3d Matrix3ds::operator * (const Matrix3da& d) const
{
    return Matrix3d(
                 -d.xy()*m[XY]-d.xz()*m[XZ],d.xy()*m[XX]-d.yz()*m[XZ],d.xz()*m[XX]+d.yz()*m[XY],
                 -d.xy()*m[YY]-d.xz()*m[YZ],d.xy()*m[XY]-d.yz()*m[YZ],d.xz()*m[XY]+d.yz()*m[YY],
                 -d.xy()*m[YZ]-d.xz()*m[ZZ],d.xy()*m[XZ]-d.yz()*m[ZZ],d.xz()*m[XZ]+d.yz()*m[YZ]
                 );
}


// unary operator -
inline Matrix3ds Matrix3ds::operator - () const
{
	return Matrix3ds(-m[XX], -m[YY], -m[ZZ], -m[XY], -m[YZ], -m[XZ]);
}

// assignment operator +=
inline Matrix3ds& Matrix3ds::operator += (const Matrix3ds& t)
{
	m[XX] += t.m[XX]; m[YY] += t.m[YY]; m[ZZ] += t.m[ZZ];
	m[XY] += t.m[XY]; m[YZ] += t.m[YZ]; m[XZ] += t.m[XZ];
	return (*this);
}

// assignment operator -=
inline Matrix3ds& Matrix3ds::operator -= (const Matrix3ds& t)
{
	m[XX] -= t.m[XX]; m[YY] -= t.m[YY]; m[ZZ] -= t.m[ZZ];
	m[XY] -= t.m[XY]; m[YZ] -= t.m[YZ]; m[XZ] -= t.m[XZ];
	return (*this);
}

// assignment operator *=
inline Matrix3ds& Matrix3ds::operator *= (const Matrix3ds& t)
{
	double xx = m[XX]*t.m[XX]+m[XY]*t.m[XY]+m[XZ]*t.m[XZ];
	double yy = m[XY]*t.m[XY]+m[YY]*t.m[YY]+m[YZ]*t.m[YZ];
	double zz = m[XZ]*t.m[XZ]+m[YZ]*t.m[YZ]+m[ZZ]*t.m[ZZ];
	double xy = m[XX]*t.m[XY]+m[XY]*t.m[YY]+m[XZ]*t.m[YZ];
	double yz = m[XY]*t.m[XZ]+m[YY]*t.m[YZ]+m[YZ]*t.m[ZZ];
	double xz = m[XX]*t.m[XZ]+m[XY]*t.m[YZ]+m[XZ]*t.m[ZZ];

	m[XX] = xx; m[YY] = yy; m[ZZ] = zz;
	m[XY] = xy; m[YZ] = yz; m[XZ] = xz;

	return (*this);
}

// assignment operator *=
inline Matrix3ds& Matrix3ds::operator *= (double g)
{
	m[XX] *= g; m[YY] *= g; m[ZZ] *= g;
	m[XY] *= g; m[YZ] *= g; m[XZ] *= g;
	return (*this);
}

// assignment operator /=
inline Matrix3ds& Matrix3ds::operator /= (double g)
{
	g = 1./g;
	m[XX] *= g; m[YY] *= g; m[ZZ] *= g;
	m[XY] *= g; m[YZ] *= g; m[XZ] *= g;
	return (*this);
}

// arithmetic assignment operators for Matrix3dd
inline Matrix3ds& Matrix3ds::operator += (const Matrix3dd& d)
{
	m[XX] += d.d[0];
	m[YY] += d.d[1];
	m[ZZ] += d.d[2];
	return (*this);
}

inline Matrix3ds& Matrix3ds::operator -= (const Matrix3dd& d)
{
    m[XX] -= d.d[0];
    m[YY] -= d.d[1];
    m[ZZ] -= d.d[2];
	return (*this);
}

// Matrix-vector multiplication
inline Vector3d Matrix3ds::operator* (const Vector3d& r) const
{
	return Vector3d(
		r.x*m[XX]+r.y*m[XY]+r.z*m[XZ],
		r.x*m[XY]+r.y*m[YY]+r.z*m[YZ],
		r.x*m[XZ]+r.y*m[YZ]+r.z*m[ZZ]
	);
}

// trace
inline double Matrix3ds::tr() const
{
	return m[XX]+m[YY]+m[ZZ];
}

// determinant
inline double Matrix3ds::det() const
{
	return (m[XX]*(m[YY]*m[ZZ] - m[YZ]*m[YZ])
		  + m[XY]*(m[YZ]*m[XZ] - m[ZZ]*m[XY])
		  + m[XZ]*(m[XY]*m[YZ] - m[YY]*m[XZ]));
}

// zero
inline void Matrix3ds::zero()
{
	m[0] = m[1] = m[2] = m[3] = m[4] = m[5] = 0;
}

// unit tensor
inline void Matrix3ds::unit()
{
	m[XX] = m[YY] = m[ZZ] = 1.0;
	m[XY] = m[YZ] = m[XZ] = 0.0;
}

// deviator
inline Matrix3ds Matrix3ds::dev() const
{
	double t = (m[XX]+m[YY]+m[ZZ])/3.0;
	return Matrix3ds(m[XX]-t, m[YY]-t, m[ZZ]-t, m[XY], m[YZ], m[XZ]);
}

// isotropic part
inline Matrix3ds Matrix3ds::iso() const
{
	double t = (m[XX]+m[YY]+m[ZZ])/3.0;
	return Matrix3ds(t, t, t, 0, 0, 0);
}

// return the square 
inline Matrix3ds Matrix3ds::sqr() const
{
	return Matrix3ds(
		m[XX] * m[XX] + m[XY] * m[XY] + m[XZ] * m[XZ],
		m[XY] * m[XY] + m[YY] * m[YY] + m[YZ] * m[YZ],
		m[XZ] * m[XZ] + m[YZ] * m[YZ] + m[ZZ] * m[ZZ],
		m[XX] * m[XY] + m[XY] * m[YY] + m[XZ] * m[YZ], 
		m[XY] * m[XZ] + m[YY] * m[YZ] + m[YZ] * m[ZZ],
		m[XX] * m[XZ] + m[XY] * m[YZ] + m[XZ] * m[ZZ]
	);
}

// inverse
inline Matrix3ds Matrix3ds::inverse() const
{
	double D = det();
	assert(D != 0);
	D = 1/D;
	
	return Matrix3ds(D*(m[YY]*m[ZZ]-m[YZ]*m[YZ]), 
				  D*(m[XX]*m[ZZ]-m[XZ]*m[XZ]), 
				  D*(m[XX]*m[YY]-m[XY]*m[XY]),
				  D*(m[XZ]*m[YZ]-m[XY]*m[ZZ]), 
				  D*(m[XY]*m[XZ]-m[XX]*m[YZ]), 
				  D*(m[XY]*m[YZ]-m[YY]*m[XZ]));
}

// invert
inline double Matrix3ds::invert(Matrix3ds& Ai)
{
    double D = det();
    if (D != 0) {
        double Di = 1/D;
        
        Ai = Matrix3ds(Di*(m[YY]*m[ZZ]-m[YZ]*m[YZ]),
                    Di*(m[XX]*m[ZZ]-m[XZ]*m[XZ]),
                    Di*(m[XX]*m[YY]-m[XY]*m[XY]),
                    Di*(m[XZ]*m[YZ]-m[XY]*m[ZZ]),
                    Di*(m[XY]*m[XZ]-m[XX]*m[YZ]),
                    Di*(m[XY]*m[YZ]-m[YY]*m[XZ]));
    }
    return D;
}

// L2-norm
inline double Matrix3ds::norm() const
{ 
	double D = m[XX]*m[XX] + m[YY]*m[YY] + m[ZZ]*m[ZZ] + 2*(m[XY]*m[XY] + m[YZ]*m[YZ] + m[XZ]*m[XZ]);
	return sqrt(D); 
}

// double contraction
inline double Matrix3ds::dotdot(const Matrix3ds& B) const
{
	const double* n = B.m;
	return m[XX]*n[XX] + m[YY]*n[YY] + m[ZZ]*n[ZZ] + 2.0*(m[XY]*n[XY] + m[YZ]*n[YZ] + m[XZ]*n[XZ]);
}

// Effective or von-mises value
inline double Matrix3ds::effective_norm() const
{
	double vm;
	vm = m[XX] * m[XX] + m[YY] * m[YY] + m[ZZ] * m[ZZ];
	vm -= m[XX] * m[YY] + m[YY] * m[ZZ] + m[XX] * m[ZZ];
	vm += 3 * (m[XY] * m[XY] + m[YZ] * m[YZ] + m[XZ] * m[XZ]);
	vm = sqrt(vm >= 0.0 ? vm : 0.0);
	return vm;
}

//-----------------------------------------------------------------------------
// class Matrix3da : anti-symmetric 3D Matrix of doubles
//-----------------------------------------------------------------------------

// constructors
inline Matrix3da::Matrix3da(double xy, double yz, double xz)
{
	d[0] = xy; d[1] = yz; d[2] = xz;
}

// calculates the antisymmetric Matrix from a vector such that for any b,
// A.b = a x b where A = Matrix3da(a).
inline Matrix3da::Matrix3da(const Vector3d& a)
{
	d[0] = -a.z; d[1] = -a.x; d[2] = a.y;
}

// access operator
inline double Matrix3da::operator ()(int i, int j) const
{
	return (i==j? 0 : (i<j? d[((j-1)<<1)-i] : -d[((i-1)<<1)-j]));
}

// scalar multiplication
inline Matrix3da Matrix3da::operator * (double g) const
{
	return Matrix3da(d[0]*g, d[1]*g, d[2]*g);
}

inline Matrix3da Matrix3da::operator + (const Matrix3da& a)
{
	return Matrix3da(d[0]+a.d[0], d[1]+a.d[1], d[2]+a.d[2]);
}

inline Matrix3da Matrix3da::operator - (const Matrix3da& a)
{
	return Matrix3da(d[0]-a.d[0], d[1]-a.d[1], d[2]-a.d[2]);
}

inline Matrix3da Matrix3da::operator - () const
{
	return Matrix3da(-d[0], -d[1], -d[2]);
}

inline Matrix3da Matrix3da::transpose() const
{
	return Matrix3da(-d[0], -d[1], -d[2]);
}

// Matrix multiplication
inline Matrix3d Matrix3da::operator * (const Matrix3d& mat3d)
{
	return Matrix3d(
		 d[0]*mat3d[1][0] + d[2]*mat3d[2][0],  d[0]*mat3d[1][1] + d[2]*mat3d[2][1],  d[0]*mat3d[1][2] + d[2]*mat3d[2][2],
		-d[0]*mat3d[0][0] + d[1]*mat3d[2][0], -d[0]*mat3d[0][1] + d[1]*mat3d[2][1], -d[0]*mat3d[0][2] + d[1]*mat3d[2][2],
		-d[2]*mat3d[0][0] - d[1]*mat3d[1][0], -d[2]*mat3d[0][1] - d[1]*mat3d[1][1], -d[2]*mat3d[0][2] - d[1]*mat3d[1][2]
	);
}

inline Matrix3d Matrix3da::operator + (const Matrix3ds& a) const
{
    return Matrix3d(
                 a.xx(),a.xy()+xy(),a.xz()+xz(),
                 a.xy()-xy(),a.yy(),a.yz()+yz(),
                 a.xz()-xz(),a.yz()-yz(),a.zz()
                 );
}

inline Matrix3d Matrix3da::operator - (const Matrix3ds& a) const
{
    return Matrix3d(
                  -a.xx(),-a.xy()+xy(),-a.xz()+xz(),
                  -a.xy()-xy(),-a.yy(),-a.yz()+yz(),
                  -a.xz()-xz(),-a.yz()-yz(),-a.zz()
                  );
}

inline Vector3d Matrix3da::operator * (const Vector3d& a)
{
	return Vector3d(
		 d[0] * a.y + d[2] * a.z, \
		-d[0] * a.x + d[1] * a.z, \
		-d[2] * a.x - d[1] * a.y);
}
