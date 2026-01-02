#pragma once
#include "RgAssembler.h"
#include "RgDomain.h"
#include <vector>

// Forward declarations
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;
class RgSolidElement;
class FEElementMatrix;

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
    virtual void ElementInternalForce(RgSolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementBodyForce(FEBodyForce& BF, RgSolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementStiffness(int iel, FEElementMatrix& ke) = 0;
    virtual void ElementMassMatrix(int iel, FEElementMatrix& ke, double scale) = 0;
    
    // Helper functions that derived classes can override
    virtual void UnpackLM(RgSolidElement& el, std::vector<int>& lm);
    
protected:
    // Data members
    RgDomain* m_domain;  // Pointer to the domain that owns this assembler
};