#pragma once

#include "datastructure/Matrix3d.h"
#include "tensor_base.h"

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class tens3ds;	// symmetric 3o tensor
class tens3drs;	// right-conjugate symmetric 3o tensor
class tens3dls;	// left-conjugate symmetric 3o tensor
class tens3d;	// general 3o tensor (no symmetry)

//-----------------------------------------------------------------------------
// traits for these classes defining the number of components
template <> class tensor_traits<tens3ds > {public: enum { NNZ = 10}; };
template <> class tensor_traits<tens3drs> {public: enum { NNZ = 18}; };
template <> class tensor_traits<tens3dls> {public: enum { NNZ = 18}; };
template <> class tensor_traits<tens3d  > {public: enum { NNZ = 27}; };

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with full symmetry Tijk = Tjik = Tkji = Tikj = Tkij = Tjki (only 10 out of 27 components are unique)

// We store this tensor as a 1x10 array.
// [T] = [T111 T112 T113 T122 T123 T133 T222 T223 T233 T333]
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9

class tens3ds : public tensor_base<tens3ds>
{
public:
	// constructors
	tens3ds(){}

	// access operator
	double operator () (int i, int j, int k) const;

	Vector3d contractdyad1(const Vector3d& v);
	double tripledot(const tens3ds& H);
};

tens3ds dyad3s(const Vector3d& l, const Vector3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with right-conjugate symmetry Gijk = Gikj (only 18 out of 27 components are unique)

// We store this tensor as a 1x18 array.
// [G] = [G111 G112 G113 G122 G123 G133 G211 G212 G213 G222 G223 G233 G311 G312 G313 G322 G323 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3drs : public tensor_base<tens3drs>
{
public:
	// constructors
	explicit tens3drs(double a);
	tens3drs(){}

	// access operator
	double  operator () (int i, int j, int k) const;
	double& operator () (int i, int j, int k);

	Vector3d contractdyad1(const Vector3d& v) const;
	Vector3d contract2s(const Matrix3ds& s) const;
	double tripledot(const tens3drs& H) const;
	Vector3d contractdyad2(const Vector3d& v, const Vector3d& w);
	tens3dls transpose();
	void contractleg2(const Matrix3d& F, int leg);
};

tens3drs operator * (const Matrix3d& F, const tens3drs& t);
tens3drs dyad3rs(const Vector3d& l, const Vector3d& r);
tens3drs dyad3rs(const Matrix3d& L, const Vector3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with left-conjugate symmetry Gijk = Gjik (only 18 out of 27 components are unique)

// We store this tensor as a 1x18 array.
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3dls : public tensor_base<tens3dls>
{
public:
	// constructors
	tens3dls(){}

	tens3dls operator * (const Matrix3d& F) const;
    tens3dls operator * (const double& f) const;

	// transpose
	tens3drs transpose();
    Vector3d trace();
    tens3d generalize();
};

tens3dls dyad3ls(const Matrix3ds& L, const Vector3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with no symmetry (27 components)

// Due to symmetry we can store this tensor as a 1x27 array.
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26

class tens3d : public tensor_base<tens3d>
{
public:
	// constructors
	tens3d(){}
    explicit tens3d(double a);

	// access operators
	double operator () (int i, int j, int k) const;
	double& operator () (int i, int j, int k);

	// return symmetric tens3ds
	tens3ds symm();
    
    // right transpose
    tens3d transposer();
    
    //Contract by 2nd order tensor
    Vector3d contract2(const Matrix3d& s) const;

    //Contract on right by vector
    Matrix3d contract1(const Vector3d& v) const;
};

tens3d operator + (const tens3dls& l, const tens3drs& r);
inline tens3d operator + (const tens3drs& r, const tens3dls& l) { return l+r; }

// The following file contains the actual definition of the class functions
#include "tens3ds.hpp"
#include "tens3drs.hpp"
#include "tens3dls.hpp"
#include "tens3d.hpp"
