#include "femcore/FETransform.h"

Transform::Transform()
{
	Reset();
}

void Transform::Reset()
{
	m_scl = Vector3d(1, 1, 1);
	m_pos = Vector3d(0, 0, 0);
	m_rot = quatd(0, 0, 0, 1);
	m_roti = quatd(0, 0, 0, 1);
}

void Transform::SetPosition(const Vector3d& t)
{
	m_pos = t;
}

void Transform::SetScale(double sx, double sy, double sz)
{
	m_scl.x = sx;
	m_scl.y = sy;
	m_scl.z = sz;
}

void Transform::SetRotation(const quatd& q)
{
	m_rot = q;
	m_roti = m_rot.Inverse();
}

void Transform::SetRotation(const Vector3d& r)
{
	m_rot = quatd(r*DEG2RAD);
	m_roti = m_rot.Inverse();
}

void Transform::SetRotation(double X, double Y, double Z)
{
	X *= DEG2RAD;
	Y *= DEG2RAD;
	Z *= DEG2RAD;
	m_rot.SetEuler(X, Y, Z);
	m_roti = m_rot.Inverse();
}

Vector3d Transform::Apply(const Vector3d& r) const
{
	Vector3d p(m_scl.x * r.x, m_scl.y * r.y, m_scl.z * r.z);
	m_rot.RotateVector(p);
	p += m_pos;
	return p;
}
