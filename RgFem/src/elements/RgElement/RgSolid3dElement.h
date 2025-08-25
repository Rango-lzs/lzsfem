#pragma once
// RgSolid3dElement.h
// Derived 3D solid element for RgFem
//
// Minimal interface derived from FESolidElement.
// Add/adjust methods to match your project's FE core.

#include "FESolidElement.h"
#include <vector>

class FEElementMatrix;
class FEElementVector;
class FEMaterialPoint;
class DumpStream;
struct vec3d;

class RgSolid3dElement : public FESolidElement
{
public:
    // constructors/destructor
    RgSolid3dElement();
    explicit RgSolid3dElement(int ntype);
    virtual ~RgSolid3dElement();

    // initialization / lifecycle
    virtual bool Init() override;
    virtual void Reset() override;
    virtual void Cleanup();

    // persistence
    virtual void Serialize(DumpStream& ar) override;

    // element calculations (high-level interfaces)
    // assemble element stiffness matrix
    virtual void AssembleStiffness(FEElementMatrix& ke);

    // assemble element mass matrix (consistent or lumped inside)
    virtual void AssembleMass(FEElementMatrix& me);

    // assemble internal (material) forces for the element
    virtual void AssembleInternalForces(FEElementVector& fe);

    // assemble external body forces (gravity, etc.)
    virtual void AssembleBodyForces(FEElementVector& fe, const vec3d& bodyForce);

    // update element state (called each step / load increment)
    virtual void Update();

    // query element result at integration point
    // stress/strain stored per integration point
    virtual void GetGaussPointStress(int ip, std::vector<double>& stress) const;
    virtual void GetGaussPointStrain(int ip, std::vector<double>& strain) const;

    // convenience utilities for subclasses or external callers
    // compute B-matrix (strain-displacement) at integration point
    virtual void CalculateBMatrix(int ip, std::vector<double>& B) const;

    // compute shape functions and derivatives at integration point
    virtual void ShapeFunctions(int ip, std::vector<double>& H) const;
    virtual void ShapeDerivatives(int ip, std::vector<std::vector<double>>& dH) const;

protected:
    // number of integration points used by this element
    int m_nint; 

    // cached shape functions and derivatives per integration point
    // layout: m_SH[ip][a] -> shape function a at ip
    std::vector<std::vector<double>> m_SH;
    // layout: m_dSH[ip][a][i] -> derivative of shape a wrt xi_i at ip
    std::vector<std::vector<std::vector<double>>> m_dSH;

    // per-integration point storage for strains/stresses (flattened)
    std::vector<std::vector<double>> m_strain;
    std::vector<std::vector<double>> m_stress;

    // helper: ensure caches sized to current integration rule
    virtual void ResizeIntegrationData(int nint);

    // low-level helpers used by Assemble* methods
    virtual void ElementStiffness(int ip, FEElementMatrix& ke_local);
    virtual void ElementMass(int ip, FEElementMatrix& me_local);
    virtual void ElementInternalForce(int ip, FEElementVector& fe_local);
};