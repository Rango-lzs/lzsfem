#pragma once
#include "datastructure/Vector3d.h"
#include "datastructure/quatd.h"
#include "datastructure/MathUtils.h"

//-----------------------------------------------------------------------------
// Class that defines an affine transformation (scale, rotate, translate).
// This currently applies the transformation as follows: 
// 1. scale : the scale is applied in the local coordinate system
// 2. rotate: rotation from local to global coordinates
// 3. translate: translate to a global position
// 
class Transform
{
public:
	Transform();

	// Reset the transform
	void Reset();

	// set the scale factors
	void SetScale(double sx, double sy, double sz);

	// set the scale of the object
	void SetScale(const Vector3d& s) { m_scl = s; }

	//! get scale of the object
	const Vector3d& GetScale() const { return m_scl; }

	// set the position (or translation)
	void SetPosition(const Vector3d& t);

	// get position of object
	const Vector3d& GetPosition() const { return m_pos; }

	// set the rotation quaternion
	void SetRotation(const quatd& q);

	// set the rotation vector (uses degrees)
	void SetRotation(const Vector3d& r);

	// set rotation via Euler angles Tait-Bryan (Z,Y,X) convention (in degrees)
	void SetRotation(double X, double Y, double Z);

	//! get orientation
	const quatd& GetRotation() const { return m_rot; }

	// get inverse of rotation
	quatd GetRotationInverse() const { return m_roti; }

	// apply transformation
	Vector3d Apply(const Vector3d& r) const;

	// translate the transform
	void Translate(const Vector3d& dr);

	// scale an object
	void Scale(double s, Vector3d r, Vector3d rc);

	// rotate around the center rc
	void Rotate(quatd q, Vector3d rc);

	// Rotate angle w around an axis defined by the position vectors a, b.
	void Rotate(const Vector3d& a, const Vector3d& b, double w);

	// comparison
	bool operator == (const Transform& T) const;

public:
	// convert from local to global coordinates
	Vector3d LocalToGlobal(const Vector3d& r) const;

	// convert from global to local coordinates
	Vector3d GlobalToLocal(const Vector3d& r) const;

	//! get a normal-like vector from global to local
	Vector3d LocalToGlobalNormal(const Vector3d& n) const;

	//! get a normal-like vector from global to local
	Vector3d GlobalToLocalNormal(const Vector3d& n) const;

private:
	Vector3d	m_scl;		// scale factors
	Vector3d	m_pos;		// translation (global space)
	quatd	m_rot;		// rotation
	quatd	m_roti;		// inverse rotation
};

inline bool Transform::operator == (const Transform& T) const
{
	return ((m_pos == T.m_pos) && (m_scl == T.m_scl) && (m_rot == T.m_rot));
}

inline void Transform::Translate(const Vector3d& dr) { m_pos += dr; }

// convert from local to global coordinates
inline Vector3d Transform::LocalToGlobal(const Vector3d& r) const
{
	return m_pos + m_rot * Vector3d(r.x * m_scl.x, r.y * m_scl.y, r.z * m_scl.z);
}

// convert from global to local coordinates
inline Vector3d Transform::GlobalToLocal(const Vector3d& r) const
{
	Vector3d p = m_roti * (r - m_pos);
	return Vector3d(p.x / m_scl.x, p.y / m_scl.y, p.z / m_scl.z);
}

//! get a normal-like vector from global to local
inline Vector3d Transform::LocalToGlobalNormal(const Vector3d& n) const
{
	// NOTE: scaling is turned off because this is used in the generation of material axes.
	//       If I use scaling the axes may no longer be orthogonal. Maybe I should create another
	//       function for this since this is now inconsistent with the reverse operation.
//		return m_rot*Vector3d(n.x / m_scl.x, n.y / m_scl.y, n.z / m_scl.z);
	return m_rot * Vector3d(n.x, n.y, n.z);
}

//! get a normal-like vector from global to local
inline Vector3d Transform::GlobalToLocalNormal(const Vector3d& n) const
{
	Vector3d m = m_roti * n;
	m.x /= m_scl.x; m.y /= m_scl.y; m.z /= m_scl.z;
	m.Normalize();
	return m;
}

// scale 
inline void Transform::Scale(double s, Vector3d r, Vector3d rc)
{
	Vector3d r0 = GlobalToLocal(rc);

	double a = s - 1;
	m_roti.RotateVector(r);
	r.Normalize();

	r.x = 1 + a * fabs(r.x);
	r.y = 1 + a * fabs(r.y);
	r.z = 1 + a * fabs(r.z);

	m_scl.x *= r.x;
	m_scl.y *= r.y;
	m_scl.z *= r.z;

	m_pos -= LocalToGlobal(r0) - rc;
}

// rotate around the center rc
inline void Transform::Rotate(quatd q, Vector3d rc)
{
	m_rot = q * m_rot;
	m_roti = m_rot.Inverse();

	m_rot.MakeUnit();
	m_roti.MakeUnit();

	m_pos = rc + q * (m_pos - rc);
}

// Rotate angle w around an axis defined by the position vectors a, b.
inline void Transform::Rotate(const Vector3d& a, const Vector3d& b, double w)
{
	double wr = PI * w / 180.0;
	Vector3d n = (b - a); n.Normalize();
	quatd q(wr, n);
	Rotate(q, a);
}
