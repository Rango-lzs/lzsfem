#include "RgHex8Element.h"
#include "elements/ElementShape/RgHex8Shape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "femcore/RgElementTraitsStore.h"
#include "../RgElementState.h"

namespace RgFem {

// ============================================================================
// Static trait initialization methods
// ============================================================================

RgElementTraits* RgHex8Element::fullIntTraits()
{
    static RgSolidElementTraits* traits = nullptr;
    if (!traits) {
        // Create traits for HEX8 element with 8 Gauss points (full integration)
        //traits = new RgSolidElementTraits(8, kNodeCount, ElementShape::ET_HEX8, ElementType::FE_HEX8G8);
        
        // Set up gauss-point coordinates and weights for hex8 element
        traits->gaussPoints.resize(8);
        
        const double rp = 0.577350269189626;  // 1/sqrt(3)
        const double rm = -rp;
        const double w = 1.0;                 // Weight for each point
        
        // Ordering: k fastest, then j, then i (i,j,k are indices in each direction)
        traits->gaussPoints[0] = RgGaussPoint(rm, rm, rm, w);
        traits->gaussPoints[1] = RgGaussPoint(rp, rm, rm, w);
        traits->gaussPoints[2] = RgGaussPoint(rp, rp, rm, w);
        traits->gaussPoints[3] = RgGaussPoint(rm, rp, rm, w);
        traits->gaussPoints[4] = RgGaussPoint(rm, rm, rp, w);
        traits->gaussPoints[5] = RgGaussPoint(rp, rm, rp, w);
        traits->gaussPoints[6] = RgGaussPoint(rp, rp, rp, w);
        traits->gaussPoints[7] = RgGaussPoint(rm, rp, rp, w);
        
        // Initialize the traits
        traits->init();
    }
    return traits;
}

RgElementTraits* RgHex8Element::reduceIntTraits()
{
    static RgSolidElementTraits* traits = nullptr;
    if (!traits) {
        // Create traits for HEX8 element with 1 Gauss point (reduced integration)
        //traits = new RgSolidElementTraits(1, kNodeCount, ElementShape::ET_HEX8, ElementType::FE_HEX8G1);
        
        // Set up gauss-point coordinates and weights for hex8 element
        traits->gaussPoints.resize(1);
        
        const double rp = 0.0;  // At center
        const double w = 8.0;   // Weight for entire element volume
        
        traits->gaussPoints[0] = RgGaussPoint(rp, rp, rp, w);
        
        // Initialize the traits
        traits->init();
    }
    return traits;
}

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8Element::RgHex8Element()
    : RgLinearSolid3dElement()
{
    // Default to full integration
    m_pTraits = RgHex8Element::fullIntTraits();
    m_node.resize(kNodeCount);
    m_loc_node.resize(kNodeCount);
}

RgHex8Element::RgHex8Element(bool fullInt)
    : RgLinearSolid3dElement()
{
    if (fullInt)
    {
        //m_pTraits = RgHex8Element::fullIntTraits();
        m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_HEX8G8);
    }
    else
    {
        //m_pTraits = RgHex8Element::reduceIntTraits();
        m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_HEX8G1);
    }
    m_node.resize(kNodeCount);
    m_loc_node.resize(kNodeCount);
}


RgHex8Element::RgHex8Element(const RgHex8Element& other)
    : RgLinearSolid3dElement(other)
{
    // m_pTraits is shared via static methods, no need to copy
    // m_node and m_loc_node will be handled by base class
}

RgHex8Element::~RgHex8Element()
{
    // Nothing specific to clean up - traits are static and managed automatically
}

RgHex8Element& RgHex8Element::operator=(const RgHex8Element& other)
{
    if (this != &other) {
        // Copy base class - this will handle nodes and other data
        RgLinearSolid3dElement::operator=(other);
        // Note: m_pTraits is not copied as it's a shared static object
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgHex8Element::elementType() const
{
    return ElementType::FE_HEX8G8;  // 8-node hex with 8 Gauss points (full integration)
}

ElementShape RgHex8Element::elementShape() const
{
    return ElementShape::ET_HEX8;
}

// ============================================================================ 
// FEM Matrix Calculations
// ============================================================================ 

void RgHex8Element::calculateStiffnessMatrix(Matrix& K) const
{
    // Get element DOFs (24 DOFs for 8 nodes with 3 DOF per node)
    int ndofs = NodeSize() * 3;
    K.resize(ndofs, ndofs);
    K.zero(); // Use setZero() for better performance with Eigen
    
    // Integrate stiffness matrix using Gauss quadrature
    int npts = GaussPointSize();
    
    for (int gp = 0; gp < npts; ++gp) {
        // Get material point to access constitutive properties at this Gauss point
        const RgMaterialPoint* matPt = getMaterialPoint(gp);
        if (!matPt || !matPt->m_pMat) continue;
        
        // Get material matrix (D matrix for isotropic linear elasticity)
        Matrix D;
        matPt->m_pMat->getConstitutiveMatrix(matPt, D);
        
        // Get gauss point coordinates
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        Vector3d natCoord(gaussPt.r, gaussPt.s, gaussPt.t);
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.w * jdet;
        
        // Compute B matrix at this Gauss point
        Matrix B;
        computeBMatrix(natCoord, B);
        
        // For linear elasticity: K += B^T * D * B * weight
        Matrix BTC = B.transpose() * D; 
        K += weight * (BTC * B); // K += B^T * D * B * weight
    }
}

void RgHex8Element::calculateMassMatrix(Matrix& M) const
{
    // Get element DOFs
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.setZero();
    
    // Get material point to access density
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    // Get material density from element material
    double rho = matPt->m_pMat->getDensity();
    
    int npts = GaussPointSize();
    
    for (int gp = 0; gp < npts; ++gp) {
        // Get gauss point coordinates from traits
        RgGaussPoint gaussPt = m_pTraits->gaussPoints[gp];
        Vector3d natCoord(gaussPt.r, gaussPt.s, gaussPt.t);
        
        // Evaluate shape functions
        std::vector<double> N;
        RgHex8Shape::evalH(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.w * jdet * rho;
        
        // Assemble consistent mass matrix: M = integral of rho * N^T * N dV
        for (int i = 0; i < kNodeCount; ++i) {
            for (int j = 0; j < kNodeCount; ++j) {
                double Mij = weight * N[i] * N[j];
                for (int d = 0; d < 3; ++d) {
                    M(3*i + d, 3*j + d) += Mij;
                }
            }
        }
    }
}

void RgHex8Element::calculateDampingMatrix(Matrix& C) const
{
    // Rayleigh damping: C = alpha*M + beta*K
    // For now, initialize as zero; actual implementation would use material damping parameters
    int ndofs = kNodeCount * 3;
    C.resize(ndofs, ndofs);
    C.zero();
}

void RgHex8Element::calculateInternalForceVector(std::vector<double>& F) const
{
    // Internal force: F = integral of B^T * sigma dV
    // This would require access to current stress state
    // Placeholder implementation
    F.clear();
    F.resize(kNodeCount * 3, 0.0);
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgHex8Element::calculateStress(FEMaterialPoint& matPt, StressTensor& stress) const
{
    // Calculate stress from strain using material model
    // Get the strain tensor from the material point
    StrainTensor strain = matPt.getStrain();
    
    // Use material model to compute stress from strain
    if (matPt.m_pMat) {
        matPt.m_pMat->calculateStress(matPt, stress);
    } else {
        // If no material model, set stress to zero
        stress.zero();
    }
}

void RgHex8Element::calculateStrain(FEMaterialPoint& matPt, StrainTensor& strain) const
{
    // Get gauss point index from material point
    int gp = matPt.m_index;
    
    // Check bounds
    if (gp < 0 || gp >= GaussPointSize()) {
        strain.zero();
        return;
    }
    
    // Get gauss point coordinates from traits
    RgGaussPoint gaussPt = m_pTraits->gaussPoints[gp];
    Vector3d natCoord(gaussPt.r, gaussPt.s, gaussPt.t);
    
    // Compute B matrix at this Gauss point
    Matrix B;
    computeBMatrix(natCoord, B);
    
    // Get element nodal displacements (assuming they are stored in the material point or accessible)
    // For now, we'll create a placeholder displacement vector
    std::vector<double> displacements(24, 0.0); // 8 nodes * 3 DOF per node
    
    // Calculate strain: ε = B * u
    // Note: In a real implementation, you would get actual displacements from the solution vector
    for (int i = 0; i < 6; ++i) { // 6 strain components
        double strain_component = 0.0;
        for (int j = 0; j < 24; ++j) { // 24 displacement DOFs
            strain_component += B(i, j) * displacements[j];
        }
        // Map to strain tensor (assuming Voigt notation: xx, yy, zz, xy, yz, xz)
        switch (i) {
            case 0: strain.xx = strain_component; break; // εxx
            case 1: strain.yy = strain_component; break; // εyy
            case 2: strain.zz = strain_component; break; // εzz
            case 3: strain.xy = strain_component; break; // γxy
            case 4: strain.yz = strain_component; break; // γyz
            case 5: strain.xz = strain_component; break; // γxz
        }
    }
}

double RgHex8Element::calculateStrainEnergy() const
{
    double strain_energy = 0.0;
    
    // Integrate strain energy over the element volume
    int npts = GaussPointSize();
    
    for (int gp = 0; gp < npts; ++gp) {
        // Get material point at this gauss point
        const FEMaterialPoint* matPt = getMaterialPoint(gp);
        if (!matPt) continue;
        
        // Get gauss point coordinates from traits
        RgGaussPoint gaussPt = m_pTraits->gaussPoints[gp];
        Vector3d natCoord(gaussPt.r, gaussPt.s, gaussPt.t);
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.w * jdet;
        
        // Get stress and strain tensors at this point
        StressTensor stress;
        StrainTensor strain;
        
        // Cast away const to call non-const methods (in real implementation, better approach needed)
        FEMaterialPoint& nonConstMatPt = const_cast<FEMaterialPoint&>(*matPt);
        calculateStress(nonConstMatPt, stress);
        calculateStrain(nonConstMatPt, strain);
        
        // Calculate strain energy density: 0.5 * σ : ε
        double strain_energy_density = 0.5 * (
            stress.xx * strain.xx +
            stress.yy * strain.yy +
            stress.zz * strain.zz +
            stress.xy * strain.xy +
            stress.yz * strain.yz +
            stress.xz * strain.xz
        );
        
        // Add to total strain energy
        strain_energy += strain_energy_density * weight;
    }
    
    return strain_energy;
}

double RgHex8Element::calculateKineticEnergy() const
{
    double kinetic_energy = 0.0;
    
    // Get material point to access density
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return 0.0;
    
    // Get material density from element material
    double rho = matPt->m_pMat->getDensity();
    
    // Integrate kinetic energy over the element volume
    int npts = GaussPointSize();
    
    for (int gp = 0; gp < npts; ++gp) {
        // Get gauss point coordinates from traits
        RgGaussPoint gaussPt = m_pTraits->gaussPoints[gp];
        Vector3d natCoord(gaussPt.r, gaussPt.s, gaussPt.t);
        
        // Evaluate shape functions
        std::vector<double> N;
        RgHex8Shape::evalH(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.w * jdet;
        
        // Get nodal velocities (assuming they are available)
        // For now, we'll create a placeholder velocity vector
        std::vector<double> velocities(24, 0.0); // 8 nodes * 3 DOF per node
        
        // Calculate velocity at gauss point by interpolation
        Vector3d velocity_gp(0.0, 0.0, 0.0);
        for (int i = 0; i < kNodeCount; ++i) {
            velocity_gp.x += N[i] * velocities[3*i + 0];
            velocity_gp.y += N[i] * velocities[3*i + 1];
            velocity_gp.z += N[i] * velocities[3*i + 2];
        }
        
        // Calculate kinetic energy density: 0.5 * rho * v^2
        double velocity_squared = velocity_gp.x*velocity_gp.x + velocity_gp.y*velocity_gp.y + velocity_gp.z*velocity_gp.z;
        double kinetic_energy_density = 0.5 * rho * velocity_squared;
        
        // Add to total kinetic energy
        kinetic_energy += kinetic_energy_density * weight;
    }
    
    return kinetic_energy;
}

// ============================================================================
// Serialization
// ============================================================================

void RgHex8Element::Serialize(DumpStream& ar)
{
    // Serialize base class
    RgLinearSolid3dElement::Serialize(ar);
}
===========================================================================

void RgHex8Element::computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const
{
    // Compute strain-displacement matrix B
    // Relates nodal displacements to strain: epsilon = B * u
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    RgHex8Shape::evalGradH(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                          dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    // B matrix has shape [6, 24] for 3D solid element
    // 6 strain components (Voigt): {exx, eyy, ezz, exy, eyz, exz}
    // 24 DOFs: 3 per node * 8 nodes
    
    B.resize(6, 24);
    B.setZero(); // Use setZero() for better performance with Eigen
    
    for (int i = 0; i < kNodeCount; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        // Using chain rule: dN/dx_i = dN/dξ_j * dξ_j/dx_i
        // where ξ_j are natural coordinates and x_i are physical coordinates
        const double dN_dx = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
        const double dN_dy = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
        const double dN_dz = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
        
        // Fill B matrix according to Voigt notation for strain components
        // ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz (engineering shear strains)
        B(0, 3*i + 0) = dN_dx;               // ε_xx = du/dx
        B(1, 3*i + 1) = dN_dy;               // ε_yy = dv/dy
        B(2, 3*i + 2) = dN_dz;               // ε_zz = dw/dz
        B(3, 3*i + 0) = dN_dy;               // γ_xy = du/dy + dv/dx
        B(3, 3*i + 1) = dN_dx;
        B(4, 3*i + 1) = dN_dz;               // γ_yz = dv/dz + dw/dy
        B(4, 3*i + 2) = dN_dy;
        B(5, 3*i + 0) = dN_dz;               // γ_xz = du/dz + dw/dx
        B(5, 3*i + 2) = dN_dx;
    }
}

