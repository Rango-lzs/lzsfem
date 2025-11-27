#include "RgHex8GeomNLElement_new.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>

namespace RgFem {

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8GeomNLElement::RgHex8GeomNLElement()
    : RgHex8Element()
{
    m_currentDisplacement.resize(kNodeCount * 3, 0.0);
    m_previousDisplacement.resize(kNodeCount * 3, 0.0);
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const std::array<int, kNodeCount>& nodeIds)
    : RgHex8Element(nodeIds)
{
    m_currentDisplacement.resize(kNodeCount * 3, 0.0);
    m_previousDisplacement.resize(kNodeCount * 3, 0.0);
}

RgHex8GeomNLElement::RgHex8GeomNLElement(const RgHex8GeomNLElement& other)
    : RgHex8Element(other),
      m_currentDisplacement(other.m_currentDisplacement),
      m_previousDisplacement(other.m_previousDisplacement),
      m_stressAtGauss(other.m_stressAtGauss),
      m_strainAtGauss(other.m_strainAtGauss),
      m_deformationGradient(other.m_deformationGradient)
{
}

RgHex8GeomNLElement::~RgHex8GeomNLElement()
{
}

RgHex8GeomNLElement& RgHex8GeomNLElement::operator=(const RgHex8GeomNLElement& other)
{
    if (this != &other) {
        RgHex8Element::operator=(other);
        m_currentDisplacement = other.m_currentDisplacement;
        m_previousDisplacement = other.m_previousDisplacement;
        m_stressAtGauss = other.m_stressAtGauss;
        m_strainAtGauss = other.m_strainAtGauss;
        m_deformationGradient = other.m_deformationGradient;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgHex8GeomNLElement::elementType() const
{
    // Indicate geometric nonlinear hex8 element
    return ElementType::FE_HEX8G8;  // Could be FE_HEX8NL if available
}

// ============================================================================
// FEM Matrix Calculations (Geometric Nonlinear)
// ============================================================================

void RgHex8GeomNLElement::calculateStiffnessMatrix(Matrix& K) const
{
    // For geometric nonlinearity: K = Km + Kg
    // Km: material (elastic) stiffness
    // Kg: geometric (stress-dependent) stiffness
    
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Call parent method for material stiffness
    RgHex8Element::calculateStiffnessMatrix(K);
    
    // Add geometric stiffness
    if (!m_stressAtGauss.empty()) {
        Matrix Kg;
        computeGeometricStiffness(m_stressAtGauss, Kg);
        K += Kg;
    }
}

void RgHex8GeomNLElement::calculateTangentStiffnessMatrix(Matrix& Kt) const
{
    // Tangent stiffness for nonlinear Newton-Raphson iteration
    // Kt = dF/du = Km + Kg(σ)
    calculateStiffnessMatrix(Kt);
}

void RgHex8GeomNLElement::calculateInternalForceVector(Vector& F) const
{
    // Internal force: F_int = integral of B^T * σ dV
    int ndofs = kNodeCount * 3;
    F.resize(ndofs);
    F.zero();
    
    if (m_stressAtGauss.empty()) {
        return;  // No stress computed yet
    }
    
    int npts = getNumberOfGaussPoints();
    
    for (int gp = 0; gp < npts; ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[gp] * jdet;
        
        // Compute B matrix at this Gauss point (nonlinear B)
        Matrix3d F = m_deformationGradient[gp];
        Matrix B_nl;
        computeNonlinearBMatrix(natCoord, F, B_nl);
        
        // Convert stress from matrix form to Voigt vector
        double sigma_voigt[6];
        const Matrix3ds& sigma = m_stressAtGauss[gp];
        sigma_voigt[0] = sigma.xx();
        sigma_voigt[1] = sigma.yy();
        sigma_voigt[2] = sigma.zz();
        sigma_voigt[3] = sigma.xy();
        sigma_voigt[4] = sigma.yz();
        sigma_voigt[5] = sigma.xz();
        
        // F_int += B^T * σ * weight
        for (int i = 0; i < ndofs; ++i) {
            for (int j = 0; j < 6; ++j) {
                F(i) += weight * B_nl(j, i) * sigma_voigt[j];
            }
        }
    }
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgHex8GeomNLElement::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    // For geometric nonlinearity, use Cauchy stress (true stress)
    // σ = (1/J) * F * S * F^T where S is second Piola-Kirchhoff stress
    
    // This would be computed from current strain at material point
    // Placeholder for actual material model integration
}

void RgHex8GeomNLElement::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    // For geometric nonlinearity, use Green-Lagrange strain
    // E = 0.5(C - I) where C = F^T * F
    
    // This would be computed from current displacement at material point
    // Placeholder for actual computation
}

// ============================================================================
// Kinematics - Deformation Gradient and Strain
// ============================================================================

void RgHex8GeomNLElement::computeDeformationGradient(const Vector3d& naturalCoord,
                                                     const std::vector<double>& nodalDispX,
                                                     const std::vector<double>& nodalDispY,
                                                     const std::vector<double>& nodalDispZ,
                                                     Matrix3d& F) const
{
    // F = I + ∂u/∂X = I + (∂u/∂x)(∂x/∂X) = I + (∂u/∂x) * J^{-1}
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    // Compute displacement gradient in physical coordinates
    Matrix3d gradu(0, 0, 0, 0, 0, 0, 0, 0, 0);
    
    for (int i = 0; i < kNodeCount; ++i) {
        double dN_dx = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
        double dN_dy = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
        double dN_dz = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
        
        gradu.m[0][0] += nodalDispX[i] * dN_dx;
        gradu.m[0][1] += nodalDispX[i] * dN_dy;
        gradu.m[0][2] += nodalDispX[i] * dN_dz;
        
        gradu.m[1][0] += nodalDispY[i] * dN_dx;
        gradu.m[1][1] += nodalDispY[i] * dN_dy;
        gradu.m[1][2] += nodalDispY[i] * dN_dz;
        
        gradu.m[2][0] += nodalDispZ[i] * dN_dx;
        gradu.m[2][1] += nodalDispZ[i] * dN_dy;
        gradu.m[2][2] += nodalDispZ[i] * dN_dz;
    }
    
    // F = I + ∇u
    F.m[0][0] = 1.0 + gradu.m[0][0];
    F.m[0][1] = gradu.m[0][1];
    F.m[0][2] = gradu.m[0][2];
    
    F.m[1][0] = gradu.m[1][0];
    F.m[1][1] = 1.0 + gradu.m[1][1];
    F.m[1][2] = gradu.m[1][2];
    
    F.m[2][0] = gradu.m[2][0];
    F.m[2][1] = gradu.m[2][1];
    F.m[2][2] = 1.0 + gradu.m[2][2];
}

void RgHex8GeomNLElement::computeGreenLagrangeStrain(const Matrix3d& F, Matrix3ds& E) const
{
    // Compute C = F^T * F
    Matrix3ds C = F.transpose() * F;
    
    // Compute E = 0.5(C - I)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            double Cij = C(i, j);
            if (i == j) {
                E(i, j) = 0.5 * (Cij - 1.0);
            } else {
                E(i, j) = 0.5 * Cij;
            }
        }
    }
}

void RgHex8GeomNLElement::computeRightCauchyGreen(const Matrix3d& F, Matrix3ds& C) const
{
    // C = F^T * F
    C = F.transpose() * F;
}

void RgHex8GeomNLElement::computeLeftCauchyGreen(const Matrix3d& F, Matrix3ds& B) const
{
    // B = F * F^T (Finger tensor)
    B = F * F.transpose();
}

void RgHex8GeomNLElement::computeEulerAlmansiStrain(const Matrix3d& F, Matrix3ds& e) const
{
    // e = 0.5(I - B^{-1}) where B = F * F^T
    Matrix3ds B = F * F.transpose();
    Matrix3ds Binv = B.inverse();
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                e(i, j) = 0.5 * (1.0 - Binv(i, j));
            } else {
                e(i, j) = -0.5 * Binv(i, j);
            }
        }
    }
}

// ============================================================================
// Material Constitutive Relations
// ============================================================================

void RgHex8GeomNLElement::computeSecondPiolaKirchhoffStress(const Matrix3ds& E, Matrix3ds& S) const
{
    // St. Venant-Kirchhoff material: S = λ*tr(E)*I + 2*μ*E
    // λ and μ are Lamé parameters (from material)
    
    // For now, use placeholder elastic material
    double lambda = 100.0;  // Would come from material model
    double mu = 80.0;       // Would come from material model
    
    double traceE = E.tr();
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            double Sij = 2.0 * mu * E(i, j);
            if (i == j) {
                Sij += lambda * traceE;
            }
            S(i, j) = Sij;
        }
    }
}

void RgHex8GeomNLElement::computeCauchyStress(const Matrix3d& F, const Matrix3ds& S, Matrix3ds& sigma) const
{
    // σ = (1/det(F)) * F * S * F^T
    double detF = F.det();
    if (std::abs(detF) < 1e-10) {
        sigma.zero();
        return;
    }
    
    // Compute F * S (product with symmetric matrix)
    Matrix3d FS;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            FS.m[i][j] = 0.0;
            for (int k = 0; k < 3; ++k) {
                FS.m[i][j] += F.m[i][k] * S(k, j);
            }
        }
    }
    
    // Compute (F * S) * F^T
    Matrix3ds temp;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j <= i; ++j) {
            double val = 0.0;
            for (int k = 0; k < 3; ++k) {
                val += FS.m[i][k] * F.m[j][k];
            }
            temp(i, j) = val / detF;
        }
    }
    sigma = temp;
}

// ============================================================================
// B-Matrix for Nonlinear Analysis
// ============================================================================

void RgHex8GeomNLElement::computeNonlinearBMatrix(const Vector3d& naturalCoord,
                                                   const Matrix3d& F,
                                                   Matrix& B_nl) const
{
    // For geometric nonlinearity, B includes nonlinear strain-displacement terms
    // This is typically used in the updated Lagrangian formulation
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    // Compute physical shape function derivatives
    std::vector<double> dN_dx(kNodeCount), dN_dy(kNodeCount), dN_dz(kNodeCount);
    for (int i = 0; i < kNodeCount; ++i) {
        dN_dx[i] = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
        dN_dy[i] = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
        dN_dz[i] = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
    }
    
    // B_nl matrix for nonlinear strain terms
    // For the updated Lagrangian formulation with Cauchy stress
    // B_nl relates incremental displacement to incremental strain
    
    int ndofs = kNodeCount * 3;
    B_nl.resize(6, ndofs);
    B_nl.zero();
    
    // Linear strain-displacement part
    for (int i = 0; i < kNodeCount; ++i) {
        B_nl(0, 3*i + 0) = dN_dx[i];  // exx
        B_nl(1, 3*i + 1) = dN_dy[i];  // eyy
        B_nl(2, 3*i + 2) = dN_dz[i];  // ezz
        B_nl(3, 3*i + 0) = dN_dy[i];  // exy
        B_nl(3, 3*i + 1) = dN_dx[i];
        B_nl(4, 3*i + 1) = dN_dz[i];  // eyz
        B_nl(4, 3*i + 2) = dN_dy[i];
        B_nl(5, 3*i + 0) = dN_dz[i];  // exz
        B_nl(5, 3*i + 2) = dN_dx[i];
    }
    
    // Nonlinear terms could be added here for specific formulations
    // e.g., for truly nonlinear strain-displacement relations
}

// ============================================================================
// Geometric Stiffness
// ============================================================================

void RgHex8GeomNLElement::computeGeometricStiffness(const std::vector<Matrix3ds>& stressAtGauss,
                                                    Matrix& Kg) const
{
    // Geometric stiffness: Kg = integral of B_geo^T * σ * B_geo dV
    // This accounts for the effect of initial stress on the stiffness
    
    int ndofs = kNodeCount * 3;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    int npts = getNumberOfGaussPoints();
    
    for (int gp = 0; gp < npts && gp < (int)stressAtGauss.size(); ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[gp] * jdet;
        
        std::vector<double> dN_dr, dN_ds, dN_dt;
        evaluateShapeDerivatives(natCoord.x, natCoord.y, natCoord.z,
                                dN_dr, dN_ds, dN_dt);
        
        Matrix3d JinvT = evaluateJacobianInverse(natCoord);
        
        std::vector<double> dN_dx(kNodeCount), dN_dy(kNodeCount), dN_dz(kNodeCount);
        for (int i = 0; i < kNodeCount; ++i) {
            dN_dx[i] = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
            dN_dy[i] = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
            dN_dz[i] = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
        }
        
        const Matrix3ds& sigma = stressAtGauss[gp];
        
        // Kg += weight * B_geo^T * σ * B_geo
        for (int i = 0; i < kNodeCount; ++i) {
            for (int j = 0; j < kNodeCount; ++j) {
                // Contribution from stress terms
                double contrib_xx = dN_dx[i] * sigma.xx() * dN_dx[j];
                double contrib_yy = dN_dy[i] * sigma.yy() * dN_dy[j];
                double contrib_zz = dN_dz[i] * sigma.zz() * dN_dz[j];
                double contrib_xy = (dN_dx[i] * sigma.xy() * dN_dy[j] +
                                    dN_dy[i] * sigma.xy() * dN_dx[j]);
                double contrib_yz = (dN_dy[i] * sigma.yz() * dN_dz[j] +
                                    dN_dz[i] * sigma.yz() * dN_dy[j]);
                double contrib_xz = (dN_dx[i] * sigma.xz() * dN_dz[j] +
                                    dN_dz[i] * sigma.xz() * dN_dx[j]);
                
                double Kg_contrib = weight * (contrib_xx + contrib_yy + contrib_zz + 
                                             contrib_xy + contrib_yz + contrib_xz);
                
                // Add to all three diagonal blocks (one for each DOF direction)
                for (int d = 0; d < 3; ++d) {
                    Kg(3*i + d, 3*j + d) += Kg_contrib;
                }
            }
        }
    }
}

// ============================================================================
// Displacement Updates
// ============================================================================

void RgHex8GeomNLElement::updateCurrentDisplacement(const std::vector<double>& displacement)
{
    m_currentDisplacement = displacement;
    
    // Cache deformation gradient and stresses at all Gauss points
    int npts = getNumberOfGaussPoints();
    m_deformationGradient.resize(npts);
    m_strainAtGauss.resize(npts);
    m_stressAtGauss.resize(npts);
    
    // Extract displacement components
    std::vector<double> ux(kNodeCount), uy(kNodeCount), uz(kNodeCount);
    getNodalDisplacements(displacement, ux, uy, uz);
    
    // Compute at each Gauss point
    for (int gp = 0; gp < npts; ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        
        // Compute deformation gradient
        computeDeformationGradient(natCoord, ux, uy, uz, m_deformationGradient[gp]);
        
        // Compute Green-Lagrange strain
        computeGreenLagrangeStrain(m_deformationGradient[gp], m_strainAtGauss[gp]);
        
        // Compute second Piola-Kirchhoff stress
        computeSecondPiolaKirchhoffStress(m_strainAtGauss[gp], m_stressAtGauss[gp]);
    }
}

void RgHex8GeomNLElement::updatePreviousDisplacement(const std::vector<double>& displacement)
{
    m_previousDisplacement = displacement;
}

// ============================================================================
// Helper Methods
// ============================================================================

void RgHex8GeomNLElement::getNodalDisplacements(const std::vector<double>& u,
                                                 std::vector<double>& ux,
                                                 std::vector<double>& uy,
                                                 std::vector<double>& uz) const
{
    ux.resize(kNodeCount);
    uy.resize(kNodeCount);
    uz.resize(kNodeCount);
    
    for (int i = 0; i < kNodeCount; ++i) {
        ux[i] = (3*i < u.size()) ? u[3*i] : 0.0;
        uy[i] = (3*i + 1 < u.size()) ? u[3*i + 1] : 0.0;
        uz[i] = (3*i + 2 < u.size()) ? u[3*i + 2] : 0.0;
    }
}

// ============================================================================
// Serialization
// ============================================================================

void RgHex8GeomNLElement::Serialize(DumpStream& ar)
{
    RgHex8Element::Serialize(ar);
    
    if (!ar.IsShallow()) {
        ar & m_currentDisplacement & m_previousDisplacement;
    }
}

} // namespace RgFem
