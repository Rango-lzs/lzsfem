#include "RgLinearSolid2dElement.h"
#include "materials/RgMaterial.h"
#include "datastructure/Vector3d.h"
#include "elements/NaturalCoord.h"
#include "materials/RgMaterialPointData.h"


// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgLinearSolid2dElement::calculateStiffnessMatrix(RgMatrix& K)
{
    // Get element DOFs (2 DOFs per node for 2D plane stress/plane strain)
    int ndofs = NodeSize() * 2;
    K.resize(ndofs, ndofs);
    K.zero();

    // Integrate stiffness matrix using Gauss quadrature
    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {
        // Get material matrix (D matrix for 2D plane stress/plane strain)
        Matrix D;
        RgMaterial* pMat = getMaterial();

        // Get material point to access constitutive properties at this Gauss point
        RgMaterialPoint* matPt = getMaterialPoint(gp);
        if (!matPt)
            continue;
        pMat->computeConstitutive(matPt, D);

        // Get gauss point coordinates
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        NaturalCoord natCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.getWeight() * jdet;

        // Compute B matrix at this Gauss point
        Matrix B;
        computeBMatrix(natCoord, B);

        // For linear elasticity: K += B^T * D * B * weight
        Matrix BTC = B.transpose() * D;
        K += (BTC * B) * weight;  // K += B^T * D * B * weight
    }
}

void RgLinearSolid2dElement::calculateMassMatrix(RgMatrix& M)
{
    // Get element DOFs
    int kNodeCount = NodeSize();
    int ndofs = NodeSize() * 2;
    M.resize(ndofs, ndofs);
    M.zero();

    // Get material density from the first material point
    const RgMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt)
        return;

    // Get material density from element material
    RgMaterial* pMat = getMaterial();
    double rho = pMat->getDensity();

    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {
        // Get gauss point coordinates from traits - using proper interface method
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        NaturalCoord natCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());

        // Evaluate shape functions
        std::vector<double> H = this->evalH(gp);

        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.getWeight() * jdet * rho;

        // Assemble consistent mass matrix: M = integral of rho * N^T * N dV
        for (int i = 0; i < kNodeCount; ++i)
        {
            for (int j = 0; j < kNodeCount; ++j)
            {
                double Mij = weight * H[i] * H[j];
                for (int d = 0; d < 2; ++d)  // 2 DOF per node in 2D
                {
                    M(2 * i + d, 2 * j + d) += Mij;
                }
            }
        }
    }
}

void RgLinearSolid2dElement::calculateDampingMatrix(RgMatrix& C)
{
    // Rayleigh damping: C = alpha*M + beta*K
    Matrix M, K;
    calculateMassMatrix(M);
    calculateStiffnessMatrix(K);

    // Get damping coefficients - using default values or from material if available
    double alpha = 0.0;  // Mass proportional damping coefficient
    double beta = 0.0;   // Stiffness proportional damping coefficient

    // Try to get damping coefficients from first material point if available
    const RgMaterialPoint* matPt = getMaterialPoint(0);
    if (matPt)
    {
        // For now using default values, in a real implementation these might come from
        // material properties or be calculated based on element characteristics
    }

    C = alpha * M + beta * K;
}

void RgLinearSolid2dElement::calculateInternalForceVector(std::vector<double>& F)
{
    // Internal force: F = integral of B^T * sigma dV
    F.clear();
    F.resize(NodeSize() * 2, 0.0);  // 2 DOF per node for 2D

    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {
        // Get gauss point coordinates
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        NaturalCoord natCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.getWeight() * jdet;

        // Compute B matrix at this Gauss point
        Matrix B;
        computeBMatrix(natCoord, B);

        // Get stress at this Gauss point
        StressTensor stress;
        RgMaterialPoint* matPt = getMaterialPoint(gp);
        if (!matPt)
            continue;

        // Calculate stress using material model
        RgMaterial* pMat = getMaterial();
        if (!pMat)
            continue;

        calculateStress(*matPt, stress);

        // Create a vector of stress components for calculation
        // (xx, yy, xy) - 2D stress state
        std::vector<double> stressVec = {stress.xx(), stress.yy(), stress.xy()};

        // Calculate internal force: F += B^T * sigma * weight
        for (int i = 0; i < B.columns(); ++i)
        {      // For each DOF
            for (int j = 0; j < B.rows(); ++j)
            {  // For each stress component
                F[i] += B(j, i) * stressVec[j] * weight;
            }
        }
    }
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgLinearSolid2dElement::getStress(RgMaterialPoint& matPt, StressTensor& stress)
{
    // Try to extract stress from the material point data if available
    // Look for specialized material point data that contains stress information
    auto* plasticData = matPt.ExtractData<SmallDef::SmallDefMaterialPointData>();
    if (plasticData)
    {
        // If we have plastic material point data, get the stress from there
        Matrix stress_matrix;
        // Convert Matrix to StressTensor (Matrix3ds)
        if (stress_matrix.rows() >= 3)  // 3 stress components in 2D: xx, yy, xy
        {
            stress.xx() = stress_matrix(0, 0);
            stress.yy() = stress_matrix(1, 0);
            stress.zz() = 0.0;  // Initialize for 2D case
            stress.xy() = stress_matrix(2, 0);
            stress.yz() = 0.0;  // Zero in 2D
            stress.xz() = 0.0;  // Zero in 2D
        }
        else
        {
            // If matrix is not in expected format, zero the stress
            stress.zero();
        }
        return;
    }

    // If no specialized data is available, calculate stress based on current state
    calculateStress(matPt, stress);
}

void RgLinearSolid2dElement::getStrain(RgMaterialPoint& matPt, StrainTensor& strain)
{
    // Try to extract strain from the material point data if available
    // Look for specialized material point data that contains strain information
    auto* plasticData = matPt.ExtractData<SmallDef::SmallDefMaterialPointData>();
    if (plasticData)
    {
        // If we have plastic material point data, get the strain from there
        Matrix strain_matrix;
        // Convert Matrix to StrainTensor (Matrix3ds)
        if (strain_matrix.rows() >= 3)  // 3 strain components in 2D: xx, yy, xy
        {
            strain.xx() = strain_matrix(0, 0);
            strain.yy() = strain_matrix(1, 0);
            strain.zz() = 0.0;  // Initialize for 2D case
            strain.xy() = strain_matrix(2, 0);
            strain.yz() = 0.0;  // Zero in 2D
            strain.xz() = 0.0;  // Zero in 2D
        }
        else
        {
            // If matrix is not in expected format, zero the strain
            strain.zero();
        }
        return;
    }

    // If no specialized data is available, calculate strain based on current state
    calculateStrain(matPt, strain);
}

void RgLinearSolid2dElement::updateStress(RgMaterialPoint& matPt)
{
    // Calculate the stress at the material point
    StressTensor stress;
    calculateStress(matPt, stress);

    // Note: In a complete implementation, we would store the stress in the material point
    // For now, we rely on the material model to update its own internal state
    RgMaterial* pMat = getMaterial();
    if (pMat)
    {
        pMat->commitState(&matPt);
    }
}

void RgLinearSolid2dElement::updateStrain(RgMaterialPoint& matPt)
{
    // Calculate the strain at the material point
    StrainTensor strain;
    calculateStrain(matPt, strain);

    // Note: In a complete implementation, we would store the strain in the material point
    // The strain is typically stored in the material point's state variables
    // For now, we rely on the material model to update its own internal state
    RgMaterial* pMat = getMaterial();
    if (pMat)
    {
        pMat->commitState(&matPt);
    }
}

double RgLinearSolid2dElement::calculateStrainEnergy()
{
    double strain_energy = 0.0;

    // Integrate strain energy over the element volume
    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {
        // Get material point at this gauss point
        const RgMaterialPoint* matPt = getMaterialPoint(gp);
        if (!matPt)
            continue;

        // Get gauss point coordinates from traits - using proper interface method
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        NaturalCoord natCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.getWeight() * jdet;

        // Get stress and strain tensors at this point
        StressTensor stress;
        StrainTensor strain;

        // Cast away const to call non-const methods (in real implementation, better approach needed)
        RgMaterialPoint& nonConstMatPt = const_cast<RgMaterialPoint&>(*matPt);
        calculateStress(nonConstMatPt, stress);
        calculateStrain(nonConstMatPt, strain);

        // Calculate strain energy density: 0.5 * σ : ε (for 2D: xx, yy, xy components)
        double strain_energy_density =
            0.5 * (stress.xx() * strain.xx() + stress.yy() * strain.yy() +
                   2.0 * stress.xy() * strain.xy());  // 2*τ_xy*γ_xy

        // Add to total strain energy
        strain_energy += strain_energy_density * weight;
    }

    return strain_energy;
}

double RgLinearSolid2dElement::calculateKineticEnergy()
{
    double kinetic_energy = 0.0;

    RgMaterial* pMat = getMaterial();
    double rho = pMat->getDensity();

    // Integrate kinetic energy over the element volume
    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {
        // Get gauss point coordinates from traits - using proper interface method
        RgGaussPoint gaussPt = m_pTraits->gaussPoint(gp);
        NaturalCoord natCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());

        // Evaluate shape functions
        std::vector<double> H = this->evalH(gp);

        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = gaussPt.getWeight() * jdet;

        // Get nodal velocities (assuming they are available)
        // For now, we'll create a placeholder velocity vector
        int numNodes = NodeSize();
        std::vector<double> velocities(numNodes * 2, 0.0);  // 2 DOF per node for 2D

        // Calculate velocity at gauss point by interpolation
        Vector3d velocity_gp(0.0, 0.0, 0.0);
        for (int i = 0; i < numNodes; ++i)
        {
            velocity_gp.x += H[i] * velocities[2 * i + 0];  // x-velocity
            velocity_gp.y += H[i] * velocities[2 * i + 1];  // y-velocity
            // z-velocity remains 0 in 2D
        }

        // Calculate kinetic energy density: 0.5 * rho * v^2
        double velocity_squared =
            velocity_gp.x * velocity_gp.x + velocity_gp.y * velocity_gp.y + velocity_gp.z * velocity_gp.z;
        double kinetic_energy_density = 0.5 * rho * velocity_squared;

        // Add to total kinetic energy
        kinetic_energy += kinetic_energy_density * weight;
    }

    return kinetic_energy;
}