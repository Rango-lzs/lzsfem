#pragma once
#include "material/FEMaterialPoint.h"
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
	mat3ds Strain() const;
	mat3ds SmallStrain() const;

	mat3ds RightCauchyGreen() const;
	mat3ds LeftCauchyGreen () const;

	mat3ds DevRightCauchyGreen() const;
	mat3ds DevLeftCauchyGreen () const;
    
    mat3ds RightStretch() const;
    mat3ds LeftStretch () const;
    
    mat3ds RightStretchInverse() const;
    mat3ds LeftStretchInverse () const;
    
    mat3ds RightHencky() const;
    mat3ds LeftHencky () const;
    
    mat3d Rotation() const;
    
    mat3ds RateOfDeformation() const { return m_L.sym(); }

	mat3ds pull_back(const mat3ds& A) const;
	mat3ds push_forward(const mat3ds& A) const;

	tens4ds pull_back(const tens4ds& C) const;
	tens4ds push_forward(const tens4ds& C) const;

public:
    bool    m_buncoupled;   //!< set to true if this material point was created by an uncoupled material
    
	// deformation data at intermediate time
	mat3d	m_F;	//!< deformation gradient
	double	m_J;	//!< determinant of F
    vec3d   m_gradJ;  //!< gradient of J
    vec3d   m_v;    //!< velocity
    vec3d   m_a;    //!< acceleration
    mat3d   m_L;    //!< spatial velocity gradient

	// solid material data
	mat3ds		m_s;		//!< Cauchy stress
    
    // uncoupled pressure
    double      m_p;        //!< only for uncoupled materials
    
    // current time data
    double      m_Wt;       //!< strain energy density at current time
    
    // previous time data
    double      m_Wp;       //!< strain energy density
};
