#pragma once
#include "femcore/fem_export.h"
#include "FEElasticDomain.h"
#include "FESolidDomain.h"
#include <vector>

// Forward declarations
class FEModel;
class FEGlobalVector;
class FEBodyForce;
class FELinearSystem;
class FESolidElement;

//-----------------------------------------------------------------------------
//! Abstract interface class for elastic solid assemblers.
class FEM_EXPORT FEElasticSolidAssembler : public FEElasticAssembler
{
public:
    FEElasticSolidAssembler(FEModel* pfem);
    virtual ~FEElasticSolidAssembler() {}

    // Implementations of FEElasticDomain interface
    void InternalForces(FEGlobalVector& R) override;
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    void InertialForces(FEGlobalVector& R, std::vector<double>& F) override;
    void StiffnessMatrix(FELinearSystem& LS) override;
    void BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf) override;
    void MassMatrix(FELinearSystem& LS, double scale) override;

protected:
    // Pure virtual functions that derived classes must implement
    virtual void ElementInternalForce(FESolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementBodyForce(FEBodyForce& BF, FESolidElement& el, std::vector<double>& fe) = 0;
    virtual void ElementStiffness(int iel, matrix& ke) = 0;
    virtual void ElementMassMatrix(int iel, matrix& ke, double scale) = 0;

protected:
    FEModel* m_pfem;
};