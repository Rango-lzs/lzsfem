#pragma once
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

#include <FECore/FEDofList.h>
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
//! domain described by Lagrange-type 3D volumetric elements
//!
class FEM_EXPORT FEElasticSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
    //! constructor
    FEElasticSolidDomain(FEModel* pfem);

    //! assignment operator
    FEElasticSolidDomain& operator=(FEElasticSolidDomain& d);

    //! activate
    void Activate() override;

    //! initialize elements
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;

    //! Unpack solid element data
    void UnpackLM(FEElement& el, vector<int>& lm) override;

    //! Set flag for update for dynamic quantities
    void SetDynamicUpdateFlag(bool b);

    //! serialization
    void Serialize(DumpStream& ar) override;

    // get the total dof list
    const FEDofList& GetDOFList() const override;

public:  // overrides from FEDomain
    //! get the material
    FEMaterial* GetMaterial() override
    {
        return m_pMat;
    }

    //! set the material
    void SetMaterial(FEMaterial* pm) override;

public:  // overrides from FEElasticDomain
    // update stresses
    void Update(const FETimeInfo& tp) override;

    // update the element stress
    virtual void UpdateElementStress(int iel, const FETimeInfo& tp);

    //! intertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, vector<double>& F) override;

    //! internal stress forces
    void InternalForces(FEGlobalVector& R) override;

    //! body forces
    void BodyForce(FEGlobalVector& R, FEBodyForce& BF) override;

    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FELinearSystem& LS) override;

    //! calculates inertial stiffness
    void MassMatrix(FELinearSystem& LS, double scale) override;

    //! body force stiffness
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;

public:
    // --- S T I F F N E S S ---

    //! calculates the solid element stiffness matrix
    virtual void ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke);

    //! geometrical stiffness (i.e. initial stress)
    virtual void ElementGeometricalStiffness(FESolidElement& el, matrix& ke);

    //! material stiffness component
    virtual void ElementMaterialStiffness(FESolidElement& el, matrix& ke);

    // --- R E S I D U A L ---

    //! Calculates the internal stress vector for solid elements
    void ElementInternalForce(FESolidElement& el, vector<double>& fe);

    //! Calculates the inertial force vector for solid elements
    void ElementInertialForce(FESolidElement& el, vector<double>& fe);

protected:
    double m_alphaf;
    double m_alpham;
    double m_beta;
    bool m_update_dynamic;  //!< flag for updating quantities only used in dynamic analysis

    bool m_secant_stress;   //!< use secant approximation to stress
    bool m_secant_tangent;  //!< flag for using secant tangent

protected:
    FEDofList m_dofU;   // displacement dofs
    FEDofList m_dofR;   // rigid rotation rofs
    FEDofList m_dofSU;  // shell displacement dofs
    FEDofList m_dofV;   // velocity dofs
    FEDofList m_dofSV;  // shell velocity dofs
    FEDofList m_dofSA;  // shell acceleration dofs
    FEDofList m_dof;    // total dof list

    FESolidMaterial* m_pMat;
};

class FEStandardElasticSolidDomain : public FEElasticSolidDomain
{
public:
    FEStandardElasticSolidDomain(FEModel* fem);

private:
    std::string m_elemType;

    DECLARE_FECORE_CLASS();
};
