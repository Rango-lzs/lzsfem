#pragma once
#include "materials/SmallDeformation/RgSmallDefMaterial.h"
#include "datastructure/Matrix.h"
#include "materials/RgMaterialPointData.h"


namespace SmallDef {

/// This class defines material point data for elastic-plastic materials.
class RgElastoPlasticMaterialPoint : public SmallDefRgMaterialPointData {
public:
    //! constructor
    RgElastoPlasticMaterialPoint(RgMaterialPointData* mp = nullptr);

    //! Initialize material point data
    void init() override;

    //! create a shallow copy
    RgMaterialPointData* copy() override;

    //! serialize material point data
    void serialize(DumpStream& ar) override;

public:
    //! Get total strain
    const Matrix& getTotalStrain() const;
    
    //! Set total strain
    void setTotalStrain(const Matrix& e);
    
    //! Get plastic strain
    const Matrix& getPlasticStrain() const;
    
    //! Set plastic strain
    void setPlasticStrain(const Matrix& ep);
    
    //! Get elastic strain
    const Matrix& getElasticStrain() const;
    
    //! Calculate elastic strain from total and plastic strains
    void calculateElasticStrain();
    
    //! Get accumulated plastic strain
    double getAccumulatedPlasticStrain() const;
    
    //! Set accumulated plastic strain
    void setAccumulatedPlasticStrain(double ep_bar);
    
    //! Get stress
    const Matrix& getStress() const;
    
    //! Set stress
    void setStress(const Matrix& stress);
    
    //! Check if the material point is in plastic state
    bool isPlastic() const;
    
    //! Set plastic state flag
    void setPlastic(bool plastic);

public:
    // Plasticity data
    Matrix m_e;         //!< Total strain (small strain tensor in Voigt notation)
    Matrix m_ep;        //!< Plastic strain (small strain tensor in Voigt notation)
    Matrix m_ee;        //!< Elastic strain (small strain tensor in Voigt notation)
    Matrix m_stress;    //!< Stress tensor (Cauchy stress in Voigt notation)
    double m_ep_bar;    //!< Accumulated plastic strain magnitude
    bool   m_bPlastic;  //!< Flag indicating if in plastic state
    double m_yield_stress; //!< Current yield stress (considering hardening)
};

} // namespace SmallDef
