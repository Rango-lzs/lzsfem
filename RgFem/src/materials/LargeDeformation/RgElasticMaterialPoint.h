#pragma once
#include "materials/FEMaterialPoint.h"
#include "datastructure/tens4ds.hpp"

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic materials.
class FEM_EXPORT FEElasticMaterialPoint : public FEMaterialPointData
{
public:
	//! constructor
	FEElasticMaterialPoint(FEMaterialPointData* mp = nullptr);

	//! Initialize material point data
	void Init() override;

	//! create a shallow copy
	FEMaterialPointData* Copy() override;

	//! serialize material point data
	void Serialize(DumpStream& ar) override;

public:
	Matrix3ds Strain() const;
	Matrix3ds SmallStrain() const;

	Matrix3ds RightCauchyGreen() const;
	Matrix3ds LeftCauchyGreen () const;

	Matrix3ds DevRightCauchyGreen() const;
	Matrix3ds DevLeftCauchyGreen () const;
    
    Matrix3ds RightStretch() const;
    Matrix3ds LeftStretch () const;
    
    Matrix3ds RightStretchInverse() const;
    Matrix3ds LeftStretchInverse () const;
    
    Matrix3ds RightHencky() const;
    Matrix3ds LeftHencky () const;
    
    Matrix3d Rotation() const;
    
    Matrix3ds RateOfDeformation() const { return m_L.sym(); }

	Matrix3ds pull_back(const Matrix3ds& A) const;
	Matrix3ds push_forward(const Matrix3ds& A) const;

	tens4ds pull_back(const tens4ds& C) const;
	tens4ds push_forward(const tens4ds& C) const;

public:
    bool    m_buncoupled;   //!< set to true if this material point was created by an uncoupled material
    
	// deformation data at intermediate time
	Matrix3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
    Vector3d   m_gradJ;  //!< gradient of J
    Vector3d   m_v;    //!< velocity
    Vector3d   m_a;    //!< acceleration
    Matrix3d   m_L;    //!< spatial velocity gradient

	// solid material data
	Matrix3ds		m_s;		//!< Cauchy stress
    
    // uncoupled pressure
    double      m_p;        //!< only for uncoupled materials
    
    // current time data
    double      m_Wt;       //!< strain energy density at current time
    
    // previous time data
    double      m_Wp;       //!< strain energy density
};
