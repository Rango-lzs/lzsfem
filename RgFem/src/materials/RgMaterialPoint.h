#pragma once

#include "femcore/fem_export.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/quatd.h"
#include "femcore/FETimeInfo.h"
#include <vector>


class RgMaterialPointData;
class RgElement;

//-----------------------------------------------------------------------------
class FEM_EXPORT RgMaterialPoint
{
public:
	RgMaterialPoint(RgMaterialPointData* data = nullptr);
	virtual ~RgMaterialPoint();

	//! The init function is used to intialize data
	virtual void Init();

	virtual RgMaterialPoint* Copy();

	//! The Update function is used to update material point data
	//! Note that this gets called at the start of the time step during PreSolveUpdate
	virtual void Update(const FETimeInfo& timeInfo);

	virtual void Serialize(DumpStream& ar);

	void Append(RgMaterialPointData* pt);

public:
	//! Extract data (\todo Is it safe for a plugin to use this function?)
	template <class T> T* ExtractData();
	template <class T> const T* ExtractData() const;

public:
	Vector3d		m_r0;		//!< material point position
	Vector3d		m_rt;		//!< current point position
	double		m_J0;		//!< reference Jacobian
	double		m_Jt;		//!< current Jacobian
	RgElement* m_elem;		//!< Element where this material point is
	int			m_index;	//!< local integration point index 

	// pointer to element's shape function values
	double* m_shape;

protected:
	RgMaterialPointData* m_data;
};


//-----------------------------------------------------------------------------
template <class T> 
inline T* RgMaterialPoint::ExtractData()
{
	return (m_data ? m_data->ExtractData<T>() : nullptr);
}

//-----------------------------------------------------------------------------
template <class T> 
inline const T* RgMaterialPoint::ExtractData() const
{
	return (m_data ? m_data->ExtractData<T>() : nullptr);
}
