#pragma once
#include "materials/FEElasticMaterialPoint.h"
#include "datastructure/tens4ds.hpp"

//-----------------------------------------------------------------------------
//! This class defines material point data for elastic-plastic materials.
class FEM_EXPORT FEElasticPlasticMaterialPoint : public FEElasticMaterialPoint
{
public:
	//! constructor
	FEElasticPlasticMaterialPoint(FEMaterialPointData* mp = nullptr);

	//! Initialize material point data
	void Init() override;

	//! create a shallow copy
	FEMaterialPointData* Copy() override;

	//! serialize material point data
	void Serialize(DumpStream& ar) override;

public:
    //! Get plastic strain
    Matrix3ds PlasticStrain() const;
    
    //! Get elastic strain
    Matrix3ds ElasticStrain() const;
    
    //! Set plastic strain
    void SetPlasticStrain(const Matrix3ds& ep);
    
    //! Set elastic strain
    void SetElasticStrain(const Matrix3ds& ee);
    
    //! Get accumulated plastic strain
    double PlasticStrainMagnitude() const;
    
    //! Set accumulated plastic strain
    void SetPlasticStrainMagnitude(double ep_mag);
    
    //! Check if the material point is in plastic state
    bool IsPlastic() const;

public:
    // Plasticity data
    Matrix3ds   m_ep;       //!< Plastic strain
    Matrix3ds   m_ee;       //!< Elastic strain
    double      m_ep_mag;   //!< Accumulated plastic strain magnitude
    bool        m_bPlastic; //!< Flag indicating if in plastic state
    double      m_yield_stress; //!< Current yield stress (considering hardening)
};