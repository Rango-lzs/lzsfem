#include "RgLinearSolid3dElement.h"
#include "materials/RgMaterial.h"
#include "datastructure/Vector3d.h"
#include "elements/NaturalCoord.h"
#include "materials/RgMaterialPointData.h"

// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgLinearSolid3dElement::calculateStiffnessMatrix(Matrix& K)
{
    // Get element DOFs (24 DOFs for 8 nodes with 3 DOF per node)
    int ndofs = NodeSize() * 3;
    K.resize(ndofs, ndofs);
    K.zero();  // Use setZero() for better performance with Eigen

    // Integrate stiffness matrix using Gauss quadrature
    int npts = GaussPointSize();

    for (int gp = 0; gp < npts; ++gp)
    {

        // Get material matrix (D matrix for isotropic linear elasticity)
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

void RgLinearSolid3dElement::calculateMassMatrix(Matrix& M)
{
    // Get element DOFs
    int kNodeCount = NodeSize();
    int ndofs = NodeSize() * 3;
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
                for (int d = 0; d < 3; ++d)
                {
                    M(3 * i + d, 3 * j + d) += Mij;
                }
            }
        }
    }
}

void RgLinearSolid3dElement::calculateDampingMatrix(Matrix& C)
{
    // Rayleigh damping: C = alpha*M + beta*K
    Matrix M, K;
    calculateMassMatrix(M);
    calculateStiffnessMatrix(K);  // Now calling the const version

    // Get damping coefficients - using default values or from material if available
    // In a real implementation, these would typically come from material properties
    double alpha = 0.0;  // Mass proportional damping coefficient
    double beta = 0.0;   // Stiffness proportional damping coefficient

    // Try to get damping coefficients from first material point if available
    const RgMaterialPoint* matPt = getMaterialPoint(0);
    if (matPt)
    {
        // For now using default values, in a real implementation these might come from
        // material properties or be calculated based on element characteristics
        // Example: alpha and beta could be computed from target damping ratios at two frequencies
    }

    C = alpha * M + beta * K;
}

void RgLinearSolid3dElement::calculateInternalForceVector(std::vector<double>& F)
{
    // Internal force: F = integral of B^T * sigma dV
    F.clear();
    F.resize(NodeSize() * 3, 0.0);

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
        // (xx, yy, zz, xy, yz, xz)
        std::vector<double> stressVec = {stress.xx(), stress.yy(), stress.zz(), stress.xy(), stress.yz(), stress.xz()};

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

void RgLinearSolid3dElement::getStress(RgMaterialPoint& matPt, StressTensor& stress)
{
    // Try to extract stress from the material point data if available
    // Look for specialized material point data that contains stress information
    auto* plasticData = matPt.ExtractData<SmallDef::SmallDefMaterialPointData>();
    if (plasticData)
    {
        // If we have plastic material point data, get the stress from there
        Matrix stress_matrix; // = plasticData->getStress();
        // Convert Matrix to StressTensor (Matrix3ds)
        if (stress_matrix.rows() >= 6)
        {
            stress.xx() = stress_matrix(0, 0);
            stress.yy() = stress_matrix(1, 0);
            stress.zz() = stress_matrix(2, 0);
            stress.xy() = stress_matrix(3, 0);
            stress.yz() = stress_matrix(4, 0);
            stress.xz() = stress_matrix(5, 0);
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

void RgLinearSolid3dElement::getStrain(RgMaterialPoint& matPt, StrainTensor& strain)
{
    // Try to extract strain from the material point data if available
    // Look for specialized material point data that contains strain information
    auto* plasticData = matPt.ExtractData<SmallDef::SmallDefMaterialPointData>();
    if (plasticData)
    {
        // If we have plastic material point data, get the strain from there
        Matrix strain_matrix; // = plasticData->getTotalStrain();
        // Convert Matrix to StrainTensor (Matrix3ds)
        if (strain_matrix.rows() >= 6)
        {
            strain.xx() = strain_matrix(0, 0);
            strain.yy() = strain_matrix(1, 0);
            strain.zz() = strain_matrix(2, 0);
            strain.xy() = strain_matrix(3, 0);
            strain.yz() = strain_matrix(4, 0);
            strain.xz() = strain_matrix(5, 0);
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

void RgLinearSolid3dElement::updateStress(RgMaterialPoint& matPt)
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

void RgLinearSolid3dElement::updateStrain(RgMaterialPoint& matPt)
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

double RgLinearSolid3dElement::calculateStrainEnergy()
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

        // Calculate strain energy density: 0.5 * ¦Ò : ¦Å
        double strain_energy_density =
            0.5 * (stress.xx() * strain.xx() + stress.yy() * strain.yy() + stress.zz() * strain.zz() +
                   stress.xy() * strain.xy() + stress.yz() * strain.yz() + stress.xz() * strain.xz());

        // Add to total strain energy
        strain_energy += strain_energy_density * weight;
    }

    return strain_energy;
}

double RgLinearSolid3dElement::calculateKineticEnergy()
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
        std::vector<double> velocities(24, 0.0);  // 8 nodes * 3 DOF per node

        // Calculate velocity at gauss point by interpolation
        Vector3d velocity_gp(0.0, 0.0, 0.0);
        for (int i = 0; i < NodeSize(); ++i)
        {
            velocity_gp.x += H[i] * velocities[3 * i + 0];
            velocity_gp.y += H[i] * velocities[3 * i + 1];
            velocity_gp.z += H[i] * velocities[3 * i + 2];
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