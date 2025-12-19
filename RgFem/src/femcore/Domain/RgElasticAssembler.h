#pragma once
#include "RgAssembler.h"
#include "RgDomain.h"
#include <vector>

// Forward declarations
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;
class FESolidElement;

//-----------------------------------------------------------------------------
//! Interface class for elastic domain assemblers.
class FEM_EXPORT RgElasticAssembler : public RgAssembler
{
public:
    RgElasticAssembler(FEModel* pfem);
    virtual ~RgElasticAssembler() {}

    // Implementations of elastic domain interface
    void InternalForces(FEGlobalVector& R) override;
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;
    void StiffnessMatrix(FELinearSystem& LS) override;
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    void MassMatrix(FELinearSystem& LS, double scale) override;

protected:
    // Virtual functions that derived classes must implement
    virtual void ElementInternalForce(FESolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementBodyForce(FEBodyForce& BF, FESolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementStiffness(int iel, Matrix& ke) = 0;
    virtual void ElementMassMatrix(int iel, Matrix& ke, double scale) = 0;
    
    // Helper functions that derived classes can override
    virtual void UnpackLM(FESolidElement& el, std::vector<int>& lm);
    
protected:
    // Data members
    RgDomain* m_domain;  // Pointer to the domain that owns this assembler
};