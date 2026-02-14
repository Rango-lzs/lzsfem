#include "RgNLSolid3dElement.h"
#include "materials/RgMaterialPoint.h"
#include "datastructure/Matrix.h"
#include "datastructure/Matrix3d.h"
#include "elements/NaturalCoord.h"
#include <vector>
#include "materials/RgMaterial.h"



void RgNLSolid3dElement::calculateMassMatrix(Matrix& M)
{
    // Nonlinear 3D mass (constant in Lagrangian description)
    // M = integral of N^T * ρ * N dV (reference configuration)
    // Mass matrix does not change during nonlinear analysis
    // Implementation would integrate over element volume using shape functions
    // Placeholder for actual implementation
}

void RgNLSolid3dElement::calculateDampingMatrix(Matrix& C)
{
    // Nonlinear damping matrix implementation
    // Could be proportional damping: C = α*M + β*K where M is mass and K is stiffness
    // Or use algorithm-specific formulation
    // Placeholder for actual implementation
}

void RgNLSolid3dElement::calculateStiffnessMatrix(Matrix& Kt)
{
}

void RgNLSolid3dElement::calculateMaterialStiffnessMatrix(Matrix& Km)
{
    // Material stiffness matrix for nonlinear analysis
    // Km = ∫ B^T * D_mat * B dV where D_mat is the material tangent modulus
    // This is the constitutive component of the tangent stiffness matrix
    
    int numNodes = NodeSize();
    int dofPerNode = 3; // 3D solid element
    int totalDOF = numNodes * dofPerNode;
    
    Km.resize(totalDOF, totalDOF);
    Km.zero();
    
    // Iterate over all Gauss points
    for (int gp = 0; gp < GaussPointSize(); gp++) {
        RgGaussPoint gaussPt = gaussPoint(gp);
        NaturalCoord naturalCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getS());
        
        // Compute Jacobian determinant at this Gauss point
        Matrix3d jacobian = evaluateJacobian(naturalCoord);
        double jacobianDet = jacobian.det();
        
        // Compute B matrix at this Gauss point
        Matrix B(6, totalDOF); // 6 strain components, totalDOF displacement components
        computeBMatrix(naturalCoord, B);
        
        // Get material tangent stiffness matrix at this Gauss point
        Matrix D_mat(6, 6);
        // Get material from the first material point or from element property
        if (getMaterialPoint(gp)) {
            RgMaterial* material = getMaterial();  //getMaterialPoint(gp)->material();
            // Get the material tangent modulus based on current state
            // This is a placeholder - actual implementation depends on material model
            material->computeConstitutive(getMaterialPoint(gp), D_mat);
        }
        
        // Integration: B^T * D_mat * B * detJ * getWeight
        Matrix temp(totalDOF, totalDOF);
        //temp.mm(B.transpose(), D_mat);           // temp = B^T * D_mat
        //temp.mm(temp, B);                    // temp = B^T * D_mat * B
        //temp.sm(jacobianDet * gaussPt.getWeight()); // temp *= detJ * getWeight
        
        // Add contribution to total stiffness matrix
        for (int i = 0; i < totalDOF; i++) {
            for (int j = 0; j < totalDOF; j++) {
                Km(i, j) += temp(i, j);
            }
        }
    }
}


void RgNLSolid3dElement::calculateGeometricStiffnessMatrix(Matrix& Kg)
{
    // Geometric stiffness: Kg = integral of B_geo^T * σ * B_geo dV
    // Accounts for effect of initial/current stress on element stiffness
    // Important for stability and buckling analysis
    
    int numNodes = NodeSize();
    int dofPerNode = 3; // 3D solid element
    int totalDOF = numNodes * dofPerNode;
    
    Kg.resize(totalDOF, totalDOF);
    Kg.zero();
    
    // Iterate over all Gauss points
    for (int gp = 0; gp < GaussPointSize(); gp++) {
        RgGaussPoint gaussPt = gaussPoint(gp);
        NaturalCoord naturalCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());
        
        // Compute Jacobian determinant at this Gauss point
        Matrix3d jacobian = evaluateJacobian(naturalCoord);
        double jacobianDet = jacobian.det();
        
        // Get stress at this Gauss point
        StressTensor stress;
        if (getMaterialPoint(gp)) {
            calculateStress(*getMaterialPoint(gp), stress);
        } else {
            // If no material point, initialize to zero stress
            stress.zero();
        }
        
        // Compute shape function derivatives at this Gauss point
        std::vector<std::vector<double>> dN_dr_ds_dt = evalDeriv(naturalCoord);
        std::vector<double> dN_dr = dN_dr_ds_dt[0];
        std::vector<double> dN_ds = dN_dr_ds_dt[1];
        std::vector<double> dN_dt = dN_dr_ds_dt[2];
        
        // Compute Jacobian inverse transpose
        Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
        JinvT = JinvT.transpose();
        
        // Build geometric stiffness matrix contribution
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                // Compute physical derivatives: dN/dx, dN/dy, dN/dz
                // Using chain rule: dN/dx_i = dN/dξ_j * dξ_j/dx_i
                const double dN_dx_i = JinvT[0][0] * dN_dr[i] + JinvT[0][1] * dN_ds[i] + JinvT[0][2] * dN_dt[i];
                const double dN_dy_i = JinvT[1][0] * dN_dr[i] + JinvT[1][1] * dN_ds[i] + JinvT[1][2] * dN_dt[i];
                const double dN_dz_i = JinvT[2][0] * dN_dr[i] + JinvT[2][1] * dN_ds[i] + JinvT[2][2] * dN_dt[i];
                
                const double dN_dx_j = JinvT[0][0] * dN_dr[j] + JinvT[0][1] * dN_ds[j] + JinvT[0][2] * dN_dt[j];
                const double dN_dy_j = JinvT[1][0] * dN_dr[j] + JinvT[1][1] * dN_ds[j] + JinvT[1][2] * dN_dt[j];
                const double dN_dz_j = JinvT[2][0] * dN_dr[j] + JinvT[2][1] * dN_ds[j] + JinvT[2][2] * dN_dt[j];
                
                // Calculate submatrix contributions
                // Row i, column j of the geometric stiffness matrix
                // The geometric stiffness matrix contains contributions from current stress state
                Kg(3*i, 3*j) += (
                    stress.xx() * dN_dx_i * dN_dx_j +
                    stress.xy() * dN_dy_i * dN_dx_j +
                    stress.xz() * dN_dz_i * dN_dx_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i, 3*j+1) += (
                    stress.xy() * dN_dx_i * dN_dy_j +
                    stress.yy() * dN_dy_i * dN_dy_j +
                    stress.yz() * dN_dz_i * dN_dy_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i, 3*j+2) += (
                    stress.xz() * dN_dx_i * dN_dz_j +
                    stress.yz() * dN_dy_i * dN_dz_j +
                    stress.zz() * dN_dz_i * dN_dz_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+1, 3*j) += (
                    stress.xx() * dN_dx_i * dN_dx_j +
                    stress.xy() * dN_dy_i * dN_dx_j +
                    stress.xz() * dN_dz_i * dN_dx_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+1, 3*j+1) += (
                    stress.xy() * dN_dx_i * dN_dy_j +
                    stress.yy() * dN_dy_i * dN_dy_j +
                    stress.yz() * dN_dz_i * dN_dy_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+1, 3*j+2) += (
                    stress.xz() * dN_dx_i * dN_dz_j +
                    stress.yz() * dN_dy_i * dN_dz_j +
                    stress.zz() * dN_dz_i * dN_dz_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+2, 3*j) += (
                    stress.xx() * dN_dx_i * dN_dx_j +
                    stress.xy() * dN_dy_i * dN_dx_j +
                    stress.xz() * dN_dz_i * dN_dx_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+2, 3*j+1) += (
                    stress.xy() * dN_dx_i * dN_dy_j +
                    stress.yy() * dN_dy_i * dN_dy_j +
                    stress.yz() * dN_dz_i * dN_dy_j
                ) * jacobianDet * gaussPt.getWeight();
                
                Kg(3*i+2, 3*j+2) += (
                    stress.xz() * dN_dx_i * dN_dz_j +
                    stress.yz() * dN_dy_i * dN_dz_j +
                    stress.zz() * dN_dz_i * dN_dz_j
                ) * jacobianDet * gaussPt.getWeight();
            }
        }
    }
}


void RgNLSolid3dElement::calculateInternalForceVector(std::vector<double>& F)
{
    // Nonlinear 3D internal force: F_int = integral of B^T * σ dV
    // Uses updated Lagrangian formulation with current configuration
    // σ = Cauchy (true) stress in current configuration
    // B is the strain-displacement matrix at each Gauss point
    
    int numNodes = NodeSize();
    int dofPerNode = 3;
    int totalDOF = numNodes * dofPerNode;
    
    F.assign(totalDOF, 0.0);
    
    // Integrate over all Gauss points
    for (int gp = 0; gp < GaussPointSize(); gp++) {
        RgGaussPoint gaussPt = gaussPoint(gp);
        NaturalCoord naturalCoord(gaussPt.getR(), gaussPt.getS(), gaussPt.getT());
        
        // Compute Jacobian determinant
        Matrix3d jacobian = evaluateJacobian(naturalCoord);
        double jacobianDet = jacobian.det();
        
        // Compute B matrix
        Matrix B(6, totalDOF);
        computeLinearB0Matrix(naturalCoord, B);
        
        // Get stress at this Gauss point
        StressTensor stress;
        if (getMaterialPoint(gp)) {
            calculateStress(*getMaterialPoint(gp), stress);
        } else {
            stress.zero();
        }
        
        // Create stress vector in Voigt notation
        std::vector<double> stress_vec(6);
        stress_vec[0] = stress.xx(); // σ_xx
        stress_vec[1] = stress.yy(); // σ_yy
        stress_vec[2] = stress.zz(); // σ_zz
        stress_vec[3] = stress.xy(); // τ_xy
        stress_vec[4] = stress.yz(); // τ_yz
        stress_vec[5] = stress.xz(); // τ_xz
        
        // Integration: B^T * stress * detJ * getWeight
        for (int i = 0; i < totalDOF; i++) {
            double Fi = 0.0;
            for (int j = 0; j < 6; j++) {
                Fi += B(j, i) * stress_vec[j];
            }
            F[i] += Fi * jacobianDet * gaussPt.getWeight();
        }
    }
}

void RgNLSolid3dElement::calculateStress(RgMaterialPoint& matPt, StressTensor& stress)
{
    // Compute stress for the nonlinear element using finite strain approach
    // This would involve computing deformation gradient, then strain, then stress
    // Placeholder for actual implementation
}

void RgNLSolid3dElement::calculateStrain(RgMaterialPoint& matPt, StrainTensor& strain)
{
    // Compute strain for the nonlinear element using finite strain approach
    // This would typically compute Green-Lagrange strain
    // Placeholder for actual implementation
}

void RgNLSolid3dElement::updateDisplacementState(const std::vector<double>& displacement)
{
    // Update internal state with new displacement values
    // Caches deformation gradients, strains, and stresses at all Gauss points
    // Used for subsequent internal force and stiffness calculations
    
    for (int gp = 0; gp < GaussPointSize(); gp++) {
        // Compute deformation gradient at this Gauss point
        Matrix3d F;
        computeDeformationGradient(gp, displacement, F);
        
        // Compute Green-Lagrange strain
        Matrix3d E;
        computeGreenLagrangeStrain(F, E);
        
        // Calculate strain tensor
        StrainTensor strain;
        //strain.set(E[0][0], E[1][1], E[2][2],
        //           E[1][2], E[0][2], E[0][1]); // ε_yz, ε_xz, ε_xy
        
        // Update material point state
        if (getMaterialPoint(gp)) {
            RgMaterialPoint* matPt = getMaterialPoint(gp);
            //matPt->setStrain(strain);
            
            // Update deformation gradient
            // This might be stored in the material point or element state
            // For now, we'll assume it's handled by the material point
        }
    }
}

void RgNLSolid3dElement::computeDeformationGradient(
    int gaussPointIndex,
    const std::vector<double>& displacement,
    Matrix3d& F)
{
    // Compute F = I + ∂u/∂X
    // F[i][j] = δ_ij + ∂u_i/∂X_j
    // At specified Gauss point using physical derivatives
    // Placeholder for actual implementation
}

void RgNLSolid3dElement::computeGreenLagrangeStrain(
    const Matrix3d& F,
    Matrix3d& E)
{
    // Compute E = 0.5(C - I) where C = F^T * F
    // Green-Lagrange strain in material (Lagrangian) description
    // Symmetric 3×3 tensor
    
    Matrix3d C = F.transpose() * F;  // Right Cauchy-Green tensor: C = F^T * F
    
    // Green-Lagrange strain: E = 1/2(C - I)
    E = C;
    E[0][0] -= 1.0;  // Subtract identity matrix
    E[1][1] -= 1.0;
    E[2][2] -= 1.0;
    E *= 0.5;          // Multiply by 1/2
}

void RgNLSolid3dElement::computeCauchyStress(
    const Matrix3d& F,
    RgMaterial& material,
    StressTensor& sigma)
{
    // Compute Cauchy (true) stress: σ = (1/det(F)) * F * S * F^T
    // S is Second Piola-Kirchhoff stress from constitutive relation
    // σ is stress in current (spatial) description
    
    // First compute the Second Piola-Kirchhoff stress
    StressTensor S;
    computePK2Stress(F, material, S);
    
    // Then transform to Cauchy stress using the formula
    double det_F = F.det();
    Matrix3d inv_F_trans = F.inverse().transpose();
    
    // σ = (1/det(F)) * F * S * F^T
    //sigma = (1.0/det_F) * (F * S * F.transpose());
}

void RgNLSolid3dElement::computePK2Stress(
    const Matrix3d& F,
    RgMaterial& material,
    StressTensor& S)
{
    // Compute Second Piola-Kirchhoff stress
    // This stress measure is work-conjugate to the Green-Lagrange strain
    // Placeholder for actual implementation
}

double RgNLSolid3dElement::calculateStrainEnergy()
{
    // Calculate strain energy for the nonlinear element
    // In nonlinear analysis, this involves integrating stress power over deformation
    // Placeholder for actual implementation
    return 0.0;
}

double RgNLSolid3dElement::calculateKineticEnergy()
{
    // Calculate kinetic energy for the nonlinear element
    // Involves velocity-dependent terms and possibly deformation-dependent mass
    // Placeholder for actual implementation
    return 0.0;
}

void RgNLSolid3dElement::computeLinearBMatrix(const NaturalCoord& xi, Matrix& B_L)
{
    // Compute the linear (infinitesimal strain) B-matrix for nonlinear analysis
    // This includes small deformation terms and accounts for initial displacement effects
    // in the updated Lagrangian formulation
    
    int numNodes = NodeSize();
    int totalDOF = numNodes * 3; // 3 DOF per node for 3D solid
    
    B_L.resize(6, totalDOF); // 6 strain components (Voigt notation), totalDOF displacement components
    B_L.zero();
    
    // Compute shape function derivatives at natural coordinates
    std::vector<std::vector<double>> dN_dr_ds_dt = evalDeriv(xi);
    std::vector<double> dN_dr = dN_dr_ds_dt[0];
    std::vector<double> dN_ds = dN_dr_ds_dt[1];
    std::vector<double> dN_dt = dN_dr_ds_dt[2];
    
    // Compute Jacobian inverse
    Matrix3d JinvT = evaluateJacobianInverse(xi);
    JinvT = JinvT.transpose();
    
    for (int i = 0; i < numNodes; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        // Using chain rule: dN/dx_i = dN/dξ_j * dξ_j/dx_i
        const double dN_dx = JinvT[0][0] * dN_dr[i] + JinvT[0][1] * dN_ds[i] + JinvT[0][2] * dN_dt[i];
        const double dN_dy = JinvT[1][0] * dN_dr[i] + JinvT[1][1] * dN_ds[i] + JinvT[1][2] * dN_dt[i];
        const double dN_dz = JinvT[2][0] * dN_dr[i] + JinvT[2][1] * dN_ds[i] + JinvT[2][2] * dN_dt[i];
        
        // Fill B matrix according to Voigt notation for strain components
        // ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz (engineering shear strains)
        B_L(0, 3*i + 0) = dN_dx;               // ∂u/∂x (normal strain in x)
        B_L(1, 3*i + 1) = dN_dy;               // ∂v/∂y (normal strain in y)
        B_L(2, 3*i + 2) = dN_dz;               // ∂w/∂z (normal strain in z)
        B_L(3, 3*i + 0) = dN_dy;               // ∂u/∂y (shear strain xy)
        B_L(3, 3*i + 1) = dN_dx;               // ∂v/∂x (shear strain xy)
        B_L(4, 3*i + 1) = dN_dz;               // ∂v/∂z (shear strain yz)
        B_L(4, 3*i + 2) = dN_dy;               // ∂w/∂y (shear strain yz)
        B_L(5, 3*i + 0) = dN_dz;               // ∂u/∂z (shear strain xz)
        B_L(5, 3*i + 2) = dN_dx;               // ∂w/∂x (shear strain xz)
    }
}

void RgNLSolid3dElement::computeLinearB0Matrix(const NaturalCoord& xi, Matrix& B_L0)
{
}

void RgNLSolid3dElement::computeLinearB1Matrix(const NaturalCoord& xi, Matrix& B_L1)
{
}

void RgNLSolid3dElement::computeNonlinearBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_NL)
{
    // Compute the nonlinear B-matrix for finite strain analysis
    // This matrix captures higher-order terms needed for geometric nonlinearity
    // and accounts for the current deformation state represented by the deformation gradient F
    
    int numNodes = NodeSize();
    int totalDOF = numNodes * 3; // 3 DOF per node for 3D solid
    
    // The nonlinear B-matrix is typically larger or different in structure than the linear one
    // For finite strain analysis, we often need to consider the material or spatial variation
    // depending on the particular formulation being used
    
    // This is a simplified implementation - in practice, the nonlinear B-matrix
    // could be structured differently depending on the specific finite strain formulation
    
    B_NL.resize(9, totalDOF); // 9 components (3x3 gradient), totalDOF displacement components
    B_NL.zero();
    
    // Compute shape function derivatives at natural coordinates
    std::vector<std::vector<double>> dN_dr_ds_dt = evalDeriv(xi);
    std::vector<double> dN_dr = dN_dr_ds_dt[0];
    std::vector<double> dN_ds = dN_dr_ds_dt[1];
    std::vector<double> dN_dt = dN_dr_ds_dt[2];
    
    // Compute Jacobian inverse
    Matrix3d JinvT = evaluateJacobianInverse(xi);
    JinvT = JinvT.transpose();
    
    for (int i = 0; i < numNodes; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        const double dN_dx = JinvT[0][0] * dN_dr[i] + JinvT[0][1] * dN_ds[i] + JinvT[0][2] * dN_dt[i];
        const double dN_dy = JinvT[1][0] * dN_dr[i] + JinvT[1][1] * dN_ds[i] + JinvT[1][2] * dN_dt[i];
        const double dN_dz = JinvT[2][0] * dN_dr[i] + JinvT[2][1] * dN_ds[i] + JinvT[2][2] * dN_dt[i];
        
        // Fill the nonlinear B-matrix
        // For geometric nonlinearity, this relates displacement gradients to deformation
        // The structure depends on the particular finite strain formulation
        
        // Row for ∂u₁/∂X₁, ∂u₁/∂X₂, ∂u₁/∂X₃ (first row of displacement gradient)
        B_NL(0, 3*i + 0) = dN_dx;   // ∂u₁/∂X₁ term
        B_NL(1, 3*i + 0) = dN_dy;   // ∂u₁/∂X₂ term
        B_NL(2, 3*i + 0) = dN_dz;   // ∂u₁/∂X₃ term
        
        // Row for ∂u₂/∂X₁, ∂u₂/∂X₂, ∂u₂/∂X₃ (second row of displacement gradient)
        B_NL(3, 3*i + 1) = dN_dx;   // ∂u₂/∂X₁ term
        B_NL(4, 3*i + 1) = dN_dy;   // ∂u₂/∂X₂ term
        B_NL(5, 3*i + 1) = dN_dz;   // ∂u₂/∂X₃ term
        
        // Row for ∂u₃/∂X₁, ∂u₃/∂X₂, ∂u₃/∂X₃ (third row of displacement gradient)
        B_NL(6, 3*i + 2) = dN_dx;   // ∂u₃/∂X₁ term
        B_NL(7, 3*i + 2) = dN_dy;   // ∂u₃/∂X₂ term
        B_NL(8, 3*i + 2) = dN_dz;   // ∂u₃/∂X₃ term
    }
}

void RgNLSolid3dElement::computeTotalLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_TL)
{
}

void RgNLSolid3dElement::computeUpdatedLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_UL)
{
}

/* 
some reference

class NonlinearElement
{
private:
    int numNodes;          // 节点数
    int numDOFPerNode;     // 每个节点的自由度数
    Matrix coordinates;    // 初始构型的节点坐标
    Vector displacements;  // 节点位移
    Matrix D;              // 材料本构矩阵

    // 高斯积分点和权重
    std::vector<Vector> gaussPoints;
    std::vector<double> gaussgetWeights;

public:
    // 计算总刚度矩阵 (TL格式)
    virtual void calculateStiffnessMatrix(Matrix& Kt)
    {
        // Kt = Km + Kg (材料刚度矩阵 + 几何刚度矩阵)
        int dof = numNodes * numDOFPerNode;
        Kt.setZero(dof, dof);

        Matrix Km(dof, dof);
        Matrix Kg(dof, dof);

        calculateMaterialStiffnessMatrix(Km);
        calculateGeometricStiffnessMatrix(Kg);

        Kt = Km + Kg;
    }

    // 计算材料刚度矩阵 (TL格式)
    virtual void calculateMaterialStiffnessMatrix(Matrix& Km)
    {
        int dof = numNodes * numDOFPerNode;
        Km.setZero(dof, dof);

        // 在初始构型上进行高斯积分
        for (size_t gp = 0; gp < gaussPoints.size(); ++gp)
        {
            Vector xi = gaussPoints[gp];
            double getWeight = gaussgetWeights[gp];

            // 1. 计算形函数对参考坐标的导数 dN/dξ
            Matrix dN_dxi = calculateShapeFunctionDerivatives(xi);

            // 2. 计算雅可比矩阵 J0 = dX/dξ (初始构型)
            Matrix J0 = dN_dxi * coordinates;
            double detJ0 = J0.determinant();
            Matrix J0_inv = J0.inverse();

            // 3. 计算形函数对初始坐标的导数 dN/dX
            Matrix dN_dX = J0_inv * dN_dxi;

            // 4. 计算位移梯度张量 H = du/dX
            Matrix H = calculateDisplacementGradient(dN_dX);

            // 5. 计算Green-Lagrange应变张量 E = 0.5*(H + H^T + H^T*H)
            Matrix E = 0.5 * (H + H.transpose() + H.transpose() * H);

            // 6. 将应变张量转换为Voigt记号的应变向量
            Vector strainVoigt = tensorToVoigt(E);

            // 7. 计算第二Piola-Kirchhoff应力 S = D * E (Voigt记号)
            Vector stressVoigt = D * strainVoigt;

            // 8. 构造B矩阵 (应变-位移矩阵的非线性部分)
            Matrix BL = calculateLinearStrainDisplacementMatrix(dN_dX);
            Matrix BNL = calculateNonlinearStrainDisplacementMatrix(dN_dX, H);
            Matrix B = BL + BNL;

            // 9. 组装材料刚度矩阵
            // Km = ∫ B^T * D * B * dV0
            Km += B.transpose() * D * B * detJ0 * getWeight;
        }
    }

    // 计算几何刚度矩阵 (TL格式)
    virtual void calculateGeometricStiffnessMatrix(Matrix& Kg)
    {
        int dof = numNodes * numDOFPerNode;
        Kg.setZero(dof, dof);

        // 在初始构型上进行高斯积分
        for (size_t gp = 0; gp < gaussPoints.size(); ++gp)
        {
            Vector xi = gaussPoints[gp];
            double getWeight = gaussgetWeights[gp];

            // 1. 计算形函数对参考坐标的导数
            Matrix dN_dxi = calculateShapeFunctionDerivatives(xi);

            // 2. 计算雅可比矩阵和行列式
            Matrix J0 = dN_dxi * coordinates;
            double detJ0 = J0.determinant();
            Matrix J0_inv = J0.inverse();

            // 3. 计算形函数对初始坐标的导数
            Matrix dN_dX = J0_inv * dN_dxi;

            // 4. 计算位移梯度
            Matrix H = calculateDisplacementGradient(dN_dX);

            // 5. 计算Green-Lagrange应变
            Matrix E = 0.5 * (H + H.transpose() + H.transpose() * H);
            Vector strainVoigt = tensorToVoigt(E);

            // 6. 计算第二Piola-Kirchhoff应力
            Vector stressVoigt = D * strainVoigt;
            Matrix S = voigtToTensor(stressVoigt);

            // 7. 构造几何刚度矩阵
            // Kg = ∫ G^T * S * G * dV0
            // 其中 G 是形函数梯度矩阵
            for (int i = 0; i < numNodes; ++i)
            {
                for (int j = 0; j < numNodes; ++j)
                {
                    // 对于3D问题
                    Matrix Kg_ij(numDOFPerNode, numDOFPerNode);
                    Kg_ij.setZero();

                    for (int k = 0; k < numDOFPerNode; ++k)
                    {
                        for (int l = 0; l < numDOFPerNode; ++l)
                        {
                            // Kg_ij(k,l) = dNi/dX_m * S_mn * dNj/dX_n
                            for (int m = 0; m < numDOFPerNode; ++m)
                            {
                                for (int n = 0; n < numDOFPerNode; ++n)
                                {
                                    Kg_ij(k, l) += dN_dX(m, i) * S(m, n) * dN_dX(n, j);
                                }
                            }
                        }
                    }

                    // 组装到全局刚度矩阵
                    for (int k = 0; k < numDOFPerNode; ++k)
                    {
                        for (int l = 0; l < numDOFPerNode; ++l)
                        {
                            Kg(i * numDOFPerNode + k, j * numDOFPerNode + l) += Kg_ij(k, l) * detJ0 * getWeight;
                        }
                    }
                }
            }
        }
    }

private:
    // 辅助函数：计算形函数导数 (需要根据具体单元类型实现)
    Matrix calculateShapeFunctionDerivatives(const Vector& xi)
    {
        // 示例：4节点四边形单元
        // 返回 dN/dξ 矩阵 (2 x 4)
        Matrix dN(2, numNodes);
        // 实现形函数导数计算...
        return dN;
    }

    // 计算位移梯度 H = du/dX
    Matrix calculateDisplacementGradient(const Matrix& dN_dX)
    {
        int dim = numDOFPerNode;
        Matrix H(dim, dim);
        H.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            for (int k = 0; k < dim; ++k)
            {
                for (int l = 0; l < dim; ++l)
                {
                    H(k, l) += displacements(i * dim + k) * dN_dX(l, i);
                }
            }
        }
        return H;
    }

    // 线性应变-位移矩阵
    Matrix calculateLinearStrainDisplacementMatrix(const Matrix& dN_dX)
    {
        // 返回 BL 矩阵 (Voigt记号)
        // 对于3D: 6 x (numNodes * 3)
        int dim = numDOFPerNode;
        int strainSize = (dim == 3) ? 6 : ((dim == 2) ? 3 : 1);
        Matrix BL(strainSize, numNodes * dim);
        BL.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            if (dim == 3)
            {
                // 3D情况
                BL(0, i * 3 + 0) = dN_dX(0, i);
                BL(1, i * 3 + 1) = dN_dX(1, i);
                BL(2, i * 3 + 2) = dN_dX(2, i);
                BL(3, i * 3 + 0) = dN_dX(1, i);
                BL(3, i * 3 + 1) = dN_dX(0, i);
                BL(4, i * 3 + 1) = dN_dX(2, i);
                BL(4, i * 3 + 2) = dN_dX(1, i);
                BL(5, i * 3 + 0) = dN_dX(2, i);
                BL(5, i * 3 + 2) = dN_dX(0, i);
            }
        }
        return BL;
    }

    // 非线性应变-位移矩阵
    Matrix calculateNonlinearStrainDisplacementMatrix(const Matrix& dN_dX, const Matrix& H)
    {
        int dim = numDOFPerNode;
        int strainSize = (dim == 3) ? 6 : ((dim == 2) ? 3 : 1);
        Matrix BNL(strainSize, numNodes * dim);
        BNL.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            if (dim == 3)
            {
                // BNL 包含位移梯度的影响
                for (int k = 0; k < 3; ++k)
                {
                    BNL(0, i * 3 + k) += H(k, 0) * dN_dX(0, i);
                    BNL(1, i * 3 + k) += H(k, 1) * dN_dX(1, i);
                    BNL(2, i * 3 + k) += H(k, 2) * dN_dX(2, i);
                    BNL(3, i * 3 + k) += H(k, 0) * dN_dX(1, i) + H(k, 1) * dN_dX(0, i);
                    BNL(4, i * 3 + k) += H(k, 1) * dN_dX(2, i) + H(k, 2) * dN_dX(1, i);
                    BNL(5, i * 3 + k) += H(k, 0) * dN_dX(2, i) + H(k, 2) * dN_dX(0, i);
                }
            }
        }
        return BNL;
    }

    // 张量到Voigt记号的转换
    Vector tensorToVoigt(const Matrix& tensor)
    {
        int dim = tensor.rows();
        Vector voigt;

        if (dim == 3)
        {
            voigt.resize(6);
            voigt << tensor(0, 0), tensor(1, 1), tensor(2, 2), tensor(0, 1), tensor(1, 2), tensor(0, 2);
        }
        else if (dim == 2)
        {
            voigt.resize(3);
            voigt << tensor(0, 0), tensor(1, 1), tensor(0, 1);
        }
        return voigt;
    }

    // Voigt记号到张量的转换
    Matrix voigtToTensor(const Vector& voigt)
    {
        int size = voigt.size();
        Matrix tensor;

        if (size == 6)
        {
            tensor.resize(3, 3);
            tensor << voigt(0), voigt(3), voigt(5), voigt(3), voigt(1), voigt(4), voigt(5), voigt(4), voigt(2);
        }
        else if (size == 3)
        {
            tensor.resize(2, 2);
            tensor << voigt(0), voigt(2), voigt(2), voigt(1);
        }
        return tensor;
    }
};


class NonlinearElementUL
{
private:
    int numNodes;               // 节点数
    int numDOFPerNode;          // 每个节点的自由度数
    Matrix coordinates;         // 初始构型的节点坐标
    Matrix currentCoordinates;  // 当前构型的节点坐标 (t时刻)
    Vector displacements;       // 从初始构型到当前构型的总位移
    Vector incrementalDisp;     // 增量位移 (从t到t+Δt)
    Matrix D;                   // 材料本构矩阵 (Jaumann率形式或其他客观应力率)

    // 当前构型的应力状态 (Cauchy应力)
    std::vector<Matrix> cauchyStress;  // 每个积分点的Cauchy应力

    // 高斯积分点和权重
    std::vector<Vector> gaussPoints;
    std::vector<double> gaussgetWeights;

public:
    // 计算总刚度矩阵 (UL格式)
    virtual void calculateStiffnessMatrix(Matrix& Kt)
    {
        // Kt = Km + Kg (材料刚度矩阵 + 几何刚度矩阵)
        int dof = numNodes * numDOFPerNode;
        Kt.setZero(dof, dof);

        Matrix Km(dof, dof);
        Matrix Kg(dof, dof);

        calculateMaterialStiffnessMatrix(Km);
        calculateGeometricStiffnessMatrix(Kg);

        Kt = Km + Kg;
    }

    // 计算材料刚度矩阵 (UL格式)
    virtual void calculateMaterialStiffnessMatrix(Matrix& Km)
    {
        int dof = numNodes * numDOFPerNode;
        Km.setZero(dof, dof);

        // 在当前构型上进行高斯积分
        for (size_t gp = 0; gp < gaussPoints.size(); ++gp)
        {
            Vector xi = gaussPoints[gp];
            double getWeight = gaussgetWeights[gp];

            // 1. 计算形函数对参考坐标的导数 dN/dξ
            Matrix dN_dxi = calculateShapeFunctionDerivatives(xi);

            // 2. 计算当前构型的雅可比矩阵 J_t = dx/dξ
            Matrix J_t = dN_dxi * currentCoordinates;
            double detJ_t = J_t.determinant();
            Matrix J_t_inv = J_t.inverse();

            // 3. 计算形函数对当前坐标的导数 dN/dx
            Matrix dN_dx = J_t_inv * dN_dxi;

            // 4. 计算增量位移梯度 h = du/dx (从t到t+Δt)
            Matrix h = calculateIncrementalDisplacementGradient(dN_dx);

            // 5. 计算增量应变张量 (基于当前构型)
            // 线性化的应变增量: Δε = 0.5*(h + h^T)
            Matrix deltaEpsilon = 0.5 * (h + h.transpose());

            // 6. 将应变张量转换为Voigt记号
            Vector deltaStrainVoigt = tensorToVoigt(deltaEpsilon);

            // 7. 构造线性应变-位移矩阵 B_L
            Matrix BL = calculateLinearStrainDisplacementMatrix(dN_dx);
        
            // 10. 获取本构矩阵 (切线模量矩阵)
            // 对于UL格式,通常使用Jaumann率或Truesdell率的本构关系
            Matrix D_tangent = getMaterialTangentMatrix(gp);

            // 11. 组装材料刚度矩阵
            // Km = ∫ B^T * D * B * dv_t
            Km += BL.transpose() * D_tangent * BL * detJ_t * getWeight;
        }
    }

    // 计算几何刚度矩阵 (UL格式)
    virtual void calculateGeometricStiffnessMatrix(Matrix& Kg)
    {
        int dof = numNodes * numDOFPerNode;
        Kg.setZero(dof, dof);

        // 在当前构型上进行高斯积分
        for (size_t gp = 0; gp < gaussPoints.size(); ++gp)
        {
            Vector xi = gaussPoints[gp];
            double getWeight = gaussgetWeights[gp];

            // 1. 计算形函数对参考坐标的导数
            Matrix dN_dxi = calculateShapeFunctionDerivatives(xi);

            // 2. 计算当前构型的雅可比矩阵和行列式
            Matrix J_t = dN_dxi * currentCoordinates;
            double detJ_t = J_t.determinant();
            Matrix J_t_inv = J_t.inverse();

            // 3. 计算形函数对当前坐标的导数 dN/dx
            Matrix dN_dx = J_t_inv * dN_dxi;

            // 4. 获取当前构型的Cauchy应力张量
            Matrix sigma = cauchyStress[gp];

            // 5. 构造几何刚度矩阵
            // Kg = ∫ G^T * σ * G * dv_t
            // 其中 G 是形函数梯度矩阵, σ 是Cauchy应力
            for (int i = 0; i < numNodes; ++i)
            {
                for (int j = 0; j < numNodes; ++j)
                {
                    // 对于3D问题
                    Matrix Kg_ij(numDOFPerNode, numDOFPerNode);
                    Kg_ij.setZero();

                    for (int k = 0; k < numDOFPerNode; ++k)
                    {
                        for (int l = 0; l < numDOFPerNode; ++l)
                        {
                            // Kg_ij(k,l) = dNi/dx_m * σ_mn * dNj/dx_n
                            for (int m = 0; m < numDOFPerNode; ++m)
                            {
                                for (int n = 0; n < numDOFPerNode; ++n)
                                {
                                    Kg_ij(k, l) += dN_dx(m, i) * sigma(m, n) * dN_dx(n, j);
                                }
                            }
                        }
                    }

                    // 组装到全局刚度矩阵
                    for (int k = 0; k < numDOFPerNode; ++k)
                    {
                        for (int l = 0; l < numDOFPerNode; ++l)
                        {
                            Kg(i * numDOFPerNode + k, j * numDOFPerNode + l) += Kg_ij(k, l) * detJ_t * getWeight;
                        }
                    }
                }
            }
        }
    }

    // 更新当前构型 (在每个增量步之后调用)
    void updateConfiguration()
    {
        // 更新当前坐标: x_t+Δt = x_t + Δu
        currentCoordinates = currentCoordinates + reshapeDisplacement(incrementalDisp);

        // 更新总位移
        displacements += incrementalDisp;

        // 更新应力状态
        updateStressState();
    }

private:
    // 辅助函数：计算形函数导数 (需要根据具体单元类型实现)
    Matrix calculateShapeFunctionDerivatives(const Vector& xi)
    {
        // 示例：4节点四边形单元
        // 返回 dN/dξ 矩阵 (2 x 4)
        Matrix dN(2, numNodes);
        // 实现形函数导数计算...
        return dN;
    }

    // 计算增量位移梯度 h = Δu/dx (基于当前构型)
    Matrix calculateIncrementalDisplacementGradient(const Matrix& dN_dx)
    {
        int dim = numDOFPerNode;
        Matrix h(dim, dim);
        h.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            for (int k = 0; k < dim; ++k)
            {
                for (int l = 0; l < dim; ++l)
                {
                    h(k, l) += incrementalDisp(i * dim + k) * dN_dx(l, i);
                }
            }
        }
        return h;
    }

    // 线性应变-位移矩阵 (基于当前构型)
    Matrix calculateLinearStrainDisplacementMatrix(const Matrix& dN_dx)
    {
        // 返回 BL 矩阵 (Voigt记号)
        // 对于3D: 6 x (numNodes * 3)
        int dim = numDOFPerNode;
        int strainSize = (dim == 3) ? 6 : ((dim == 2) ? 3 : 1);
        Matrix BL(strainSize, numNodes * dim);
        BL.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            if (dim == 3)
            {
                // 3D情况 - 基于当前构型的导数
                BL(0, i * 3 + 0) = dN_dx(0, i);
                BL(1, i * 3 + 1) = dN_dx(1, i);
                BL(2, i * 3 + 2) = dN_dx(2, i);
                BL(3, i * 3 + 0) = dN_dx(1, i);
                BL(3, i * 3 + 1) = dN_dx(0, i);
                BL(4, i * 3 + 1) = dN_dx(2, i);
                BL(4, i * 3 + 2) = dN_dx(1, i);
                BL(5, i * 3 + 0) = dN_dx(2, i);
                BL(5, i * 3 + 2) = dN_dx(0, i);
            }
            else if (dim == 2)
            {
                // 2D情况
                BL(0, i * 2 + 0) = dN_dx(0, i);
                BL(1, i * 2 + 1) = dN_dx(1, i);
                BL(2, i * 2 + 0) = dN_dx(1, i);
                BL(2, i * 2 + 1) = dN_dx(0, i);
            }
        }
        return BL;
    }

    // 非线性应变-位移矩阵 (考虑增量位移梯度的影响)
    Matrix calculateNonlinearStrainDisplacementMatrix(const Matrix& dN_dx, const Matrix& h)
    {
        int dim = numDOFPerNode;
        int strainSize = (dim == 3) ? 6 : ((dim == 2) ? 3 : 1);
        Matrix BNL(strainSize, numNodes * dim);
        BNL.setZero();

        for (int i = 0; i < numNodes; ++i)
        {
            if (dim == 3)
            {
                // BNL 包含增量位移梯度的影响
                for (int k = 0; k < 3; ++k)
                {
                    BNL(0, i * 3 + k) += h(k, 0) * dN_dx(0, i);
                    BNL(1, i * 3 + k) += h(k, 1) * dN_dx(1, i);
                    BNL(2, i * 3 + k) += h(k, 2) * dN_dx(2, i);
                    BNL(3, i * 3 + k) += h(k, 0) * dN_dx(1, i) + h(k, 1) * dN_dx(0, i);
                    BNL(4, i * 3 + k) += h(k, 1) * dN_dx(2, i) + h(k, 2) * dN_dx(1, i);
                    BNL(5, i * 3 + k) += h(k, 0) * dN_dx(2, i) + h(k, 2) * dN_dx(0, i);
                }
            }
            else if (dim == 2)
            {
                // 2D情况
                for (int k = 0; k < 2; ++k)
                {
                    BNL(0, i * 2 + k) += h(k, 0) * dN_dx(0, i);
                    BNL(1, i * 2 + k) += h(k, 1) * dN_dx(1, i);
                    BNL(2, i * 2 + k) += h(k, 0) * dN_dx(1, i) + h(k, 1) * dN_dx(0, i);
                }
            }
        }
        return BNL;
    }

    // 获取材料切线模量矩阵 (Jaumann率或其他客观应力率)
    Matrix getMaterialTangentMatrix(int gaussPoint)
    {
        // 这里需要根据具体的本构模型实现
        // 对于UL格式,通常使用Jaumann率的本构关系
        // D_jaumann = D_material + correction terms

        // 简单情况:返回预设的D矩阵
        return D;

        // 复杂情况:需要考虑应力旋转的影响
        // 可能需要实现 Jaumann 校正或 Truesdell 校正
    }

    // 更新应力状态
    void updateStressState()
    {
        // 在每个积分点更新Cauchy应力
        for (size_t gp = 0; gp < gaussPoints.size(); ++gp)
        {
            Vector xi = gaussPoints[gp];

            // 计算形函数导数
            Matrix dN_dxi = calculateShapeFunctionDerivatives(xi);
            Matrix J_t = dN_dxi * currentCoordinates;
            Matrix J_t_inv = J_t.inverse();
            Matrix dN_dx = J_t_inv * dN_dxi;

            // 计算增量应变
            Matrix h = calculateIncrementalDisplacementGradient(dN_dx);
            Matrix deltaEpsilon = 0.5 * (h + h.transpose());
            Vector deltaStrainVoigt = tensorToVoigt(deltaEpsilon);

            // 计算应力增量
            Vector deltaStressVoigt = D * deltaStrainVoigt;
            Matrix deltaStress = voigtToTensor(deltaStressVoigt);

            // 更新Cauchy应力 (简化版本,实际应考虑应力旋转)
            // σ_t+Δt = σ_t + Δσ + rotation correction
            cauchyStress[gp] = cauchyStress[gp] + deltaStress;

            // 更精确的更新应该考虑Jaumann率:
            // σ_t+Δt = σ_t + Δσ + σ_t*W - W*σ_t
            // 其中 W = 0.5*(h - h^T) 是旋转张量
            Matrix W = 0.5 * (h - h.transpose());
            cauchyStress[gp] = cauchyStress[gp] + deltaStress + cauchyStress[gp] * W - W * cauchyStress[gp];
        }
    }

    // 张量到Voigt记号的转换
    Vector tensorToVoigt(const Matrix& tensor)
    {
        int dim = tensor.rows();
        Vector voigt;

        if (dim == 3)
        {
            voigt.resize(6);
            voigt << tensor(0, 0), tensor(1, 1), tensor(2, 2), tensor(0, 1), tensor(1, 2), tensor(0, 2);
        }
        else if (dim == 2)
        {
            voigt.resize(3);
            voigt << tensor(0, 0), tensor(1, 1), tensor(0, 1);
        }
        return voigt;
    }

    // Voigt记号到张量的转换
    Matrix voigtToTensor(const Vector& voigt)
    {
        int size = voigt.size();
        Matrix tensor;

        if (size == 6)
        {
            tensor.resize(3, 3);
            tensor << voigt(0), voigt(3), voigt(5), voigt(3), voigt(1), voigt(4), voigt(5), voigt(4), voigt(2);
        }
        else if (size == 3)
        {
            tensor.resize(2, 2);
            tensor << voigt(0), voigt(2), voigt(2), voigt(1);
        }
        return tensor;
    }

    // 将位移向量重塑为矩阵形式
    Matrix reshapeDisplacement(const Vector& disp)
    {
        Matrix dispMatrix(numNodes, numDOFPerNode);
        for (int i = 0; i < numNodes; ++i)
        {
            for (int j = 0; j < numDOFPerNode; ++j)
            {
                dispMatrix(i, j) = disp(i * numDOFPerNode + j);
            }
        }
        return dispMatrix;
    }
};

*/
