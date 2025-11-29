#ifndef RGHEX8GEOMNLELEMENT_H
#define RGHEX8GEOMNLELEMENT_H

#include "RgHex8Element.h"
#include <array>
#include <vector>

namespace RgFem {

// forward declarations
class RgMaterial;
class Matrix3d;

/*
    RgHex8GeomNLElement
    - 8-node linear hexahedral element with geometric nonlinearity
    - Derives from RgHex8Element (inherits linear shape functions)
    - Accounts for large deformations and rotations
    - Uses updated Lagrangian formulation
    - Features:
      * Deformation gradient F computation
      * Green-Lagrange strain (material description)
      * Second Piola-Kirchhoff stress (material description)
      * Geometric (initial stress) stiffness matrix
      * Supports hyperelastic and elastoplastic materials
*/
class RgHex8GeomNLElement : public RgNLSolid3dElement
{
public:
    static constexpr int kNodeCount = 8;

    // Constructors and Destructors
    RgHex8GeomNLElement();
    explicit RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds);
    RgHex8GeomNLElement(const RgHex8GeomNLElement& other);
    RgHex8GeomNLElement& operator=(const RgHex8GeomNLElement& other);
    virtual ~RgHex8GeomNLElement();

    // Element Type Identification - Override to indicate geometric nonlinearity
    virtual ElementType elementType() const override;
    
    // Override stiffness calculations for geometric nonlinearity
    virtual void calculateStiffnessMatrix(Matrix& K) const override;
    virtual void calculateTangentStiffnessMatrix(Matrix& Kt) const override;
    virtual void calculateInternalForceVector(Vector& F) const override;

    // Override strain/stress calculations for large deformations
    virtual void calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress) override;
    virtual void calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain) override;

    // Nonlinear mechanics specific methods
    // Compute deformation gradient F = ∂x/∂X = I + ∂u/∂X
    void computeDeformationGradient(const Vector3d& naturalCoord,
                                     const std::vector<double>& nodalDispX,
                                     const std::vector<double>& nodalDispY,
                                     const std::vector<double>& nodalDispZ,
                                     Matrix3d& F) const;

    // Compute Green-Lagrange strain E = 0.5(C - I) where C = F^T * F
    void computeGreenLagrangeStrain(const Matrix3d& F, Matrix3ds& E) const;

    // Compute right Cauchy-Green deformation tensor C = F^T * F
    void computeRightCauchyGreen(const Matrix3d& F, Matrix3ds& C) const;

    // Compute left Cauchy-Green deformation tensor (Finger tensor) B = F * F^T
    void computeLeftCauchyGreen(const Matrix3d& F, Matrix3ds& B) const;

    // Compute Euler-Almansi strain e = 0.5(I - B^{-1})
    void computeEulerAlmansiStrain(const Matrix3d& F, Matrix3ds& e) const;

    // Material stress-strain relationship for hyperelastic materials
    // Returns second Piola-Kirchhoff stress S from Green-Lagrange strain E
    void computeSecondPiolaKirchhoffStress(const Matrix3ds& E, Matrix3ds& S) const;

    // Compute Cauchy (true) stress σ = (1/det(F)) * F * S * F^T
    void computeCauchyStress(const Matrix3d& F, const Matrix3ds& S, Matrix3ds& sigma) const;

    // Modified B-matrix for geometric nonlinearity
    // For nonlinear analysis, includes nonlinear strain terms
    void computeNonlinearBMatrix(const Vector3d& naturalCoord, const Matrix3d& F, Matrix& B_nl) const;

    // Geometric (initial stress) stiffness matrix
    // Kg = integral of B_geo^T * sigma * B_geo dV
    void computeGeometricStiffness(const std::vector<Matrix3ds>& stressAtGauss, Matrix& Kg) const;

    // Update displacement field for current iteration
    void updateCurrentDisplacement(const std::vector<double>& displacement);
    const std::vector<double>& getCurrentDisplacement() const { return m_currentDisplacement; }

    // Displacement at previous converged state
    void updatePreviousDisplacement(const std::vector<double>& displacement);
    const std::vector<double>& getPreviousDisplacement() const { return m_previousDisplacement; }

    // Serialization
    virtual void Serialize(DumpStream& ar) override;

protected:
    // Helper: Extract nodal displacement components
    void getNodalDisplacements(const std::vector<double>& u,
                               std::vector<double>& ux,
                               std::vector<double>& uy,
                               std::vector<double>& uz) const;

    // Helper: Compute displacement derivatives at Gauss point
    void computeDisplacementGradient(const Vector3d& naturalCoord,
                                      const std::vector<double>& nodalDispX,
                                      const std::vector<double>& nodalDispY,
                                      const std::vector<double>& nodalDispZ,
                                      Matrix3d& gradu) const;

private:
    // Current nodal displacement (for updated Lagrangian formulation)
    std::vector<double> m_currentDisplacement;
    
    // Previous converged displacement (for incremental analysis)
    std::vector<double> m_previousDisplacement;

    // Cached stress and strain at Gauss points
    mutable std::vector<Matrix3ds> m_stressAtGauss;
    mutable std::vector<Matrix3ds> m_strainAtGauss;
    mutable std::vector<Matrix3d> m_deformationGradient;
};

} // namespace RgFem

#endif // RGHEX8GEOMNLELEMENT_H
