#pragma once

#include "RgSolid3dElement.h"

//! RgNLSolid3dElement - Nonlinear 3D solid element
//! Used for large strain, large displacement analysis (geometrically nonlinear)
//! Features:
//!   - Finite strain formulation
//!   - Updated Lagrangian description
//!   - Deformation gradient F = I + ∂u/∂X
//!   - Green-Lagrange strain E = 0.5(C - I)
//!   - Second Piola-Kirchhoff stress computation
//!   - Geometric stiffness from initial stress effects
//!   - Compatible with hyperelastic and elastoplastic materials
class FEM_EXPORT RgNLSolid3dElement : public RgSolid3dElement
{
public:
    RgNLSolid3dElement() = default;
    RgNLSolid3dElement(const RgNLSolid3dElement& el) = default;
    RgNLSolid3dElement& operator=(const RgNLSolid3dElement& el) = default;
    virtual ~RgNLSolid3dElement() = default;

    // Mass and Damping Matrices
    virtual void calculateMassMatrix(Matrix& M) override;
    virtual void calculateDampingMatrix(Matrix& C) override;

    virtual void calculateStiffnessMatrix(Matrix& Kt);

    // Stiffness matrices for nonlinear analysis
    virtual void calculateMaterialStiffnessMatrix(Matrix& Km);
    // Nonlinear-specific methods
    virtual void calculateGeometricStiffnessMatrix(Matrix& Kg);
    

    // Internal force vector (residual)
    virtual void calculateInternalForceVector(std::vector<double>& F) override;

    // Strain and Stress Calculations
    virtual void calculateStress(RgMaterialPoint& matPt, StressTensor& stress) override;
    virtual void calculateStrain(RgMaterialPoint& matPt, StrainTensor& strain) override;

   

    // Update displacement state for nonlinear iteration
    virtual void updateDisplacementState(const std::vector<double>& displacement);

    // Compute deformation gradient F = I + ∂u/∂X
    virtual void computeDeformationGradient(int gaussPointIndex, const std::vector<double>& displacement, Matrix3d& F);

    // Compute Green-Lagrange strain E = 0.5(C - I) where C = F^T * F
    virtual void computeGreenLagrangeStrain(const Matrix3d& F, Matrix3d& E);

    // Compute Cauchy stress from deformation and material
    virtual void computeCauchyStress(const Matrix3d& F, RgMaterial& material, StressTensor& sigma);

    // Compute Second Piola-Kirchhoff stress
    virtual void computePK2Stress(const Matrix3d& F, RgMaterial& material, StressTensor& S);

    // Strain energy calculation
    virtual double calculateStrainEnergy() override;

    // Kinetic energy calculation
    virtual double calculateKineticEnergy() override;

protected:
    // ========== B矩阵相关接口 ==========
    virtual void computeBMatrix(const NaturalCoord& xi, Matrix& B) const override = 0;

    //refer to Bathe  ： THe B matrix divided to linear and quadratic of the u incrementation
    virtual void computeLinearBMatrix(const NaturalCoord& xi, Matrix& B_L0);
    virtual void computeLinearB0Matrix(const NaturalCoord& xi, Matrix& B_L0);
    virtual void computeLinearB1Matrix(const NaturalCoord& xi, Matrix& B_L1);
    virtual void computeNonlinearBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_NL);

    virtual void computeTotalLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_TL);
    virtual void computeUpdatedLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_UL);
};


/*  reference

class FEM_EXPORT RgNLSolid3dElementAi : public RgSolid3dElement
{
public:
    // ========== 构造与析构 ==========
    RgNLSolid3dElementAi() = default;
    virtual ~RgNLSolid3dElementAi() = default;

    // ========== 状态管理接口 ==========
    virtual void initializeState();
    virtual void saveStateAtIncrementStart();
    virtual void restoreStateAtIncrementStart();
    virtual void commitState();
    virtual void revertToLastCommittedState();

    // ========== 位移和变形接口 ==========
    virtual void updateDisplacementState(const std::vector<double>& displacement) override;
    virtual void updateDisplacementIncrement(const std::vector<double>& displacementIncrement);

    virtual void computeDeformationGradient(int gaussPointIndex, const std::vector<double>& displacement, Matrix3d& F);

    virtual void computeIncrementalDeformationGradient(int gaussPointIndex,
                                                       const std::vector<double>& displacementIncrement,
                                                       Matrix3d& deltaF);

    // ========== 应变计算接口 ==========
    virtual void calculateStrain(RgMaterialPoint& matPt, StrainTensor& strain) override;
    virtual void calculateStrainIncrement(int gaussPointIndex, const std::vector<double>& displacementIncrement,
                                          StrainTensor& strainInc);

    virtual void computeGreenLagrangeStrain(const Matrix3d& F, Matrix3d& E);
    virtual void computeLogarithmicStrain(const Matrix3d& F, Matrix3d& eps);
    virtual void computeAlmansiStrain(const Matrix3d& F, Matrix3d& e);

    // ========== 应力计算接口 ==========
    virtual void calculateStress(RgMaterialPoint& matPt, StressTensor& stress) override;
    virtual void calculateStressIncrement(const StrainTensor& strainInc, RgMaterialPoint& matPt,
                                          StressTensor& stressInc);

    virtual void computeCauchyStress(const Matrix3d& F, RgMaterial& material, StressTensor& sigma);
    virtual void computePK2Stress(const Matrix3d& F, RgMaterial& material, StressTensor& S);
    virtual void computeFirstPiolaKirchhoffStress(const Matrix3d& F, RgMaterial& material, StressTensor& P);

    // ========== 刚度矩阵接口 ==========
    virtual void calculateStiffnessMatrix(Matrix& K) override;
    virtual void calculateTangentStiffnessMatrix(Matrix& Kt);
    virtual void calculateMaterialStiffnessMatrix(Matrix& Km);
    virtual void calculateGeometricStiffnessMatrix(Matrix& Kg);
    virtual void calculateConsistentTangentMatrix(Matrix& Ct);

    // ========== 内力向量接口 ==========
    virtual void calculateInternalForceVector(std::vector<double>& F) override;
    virtual void calculateInternalForceVectorTL(std::vector<double>& F);  // Total Lagrangian格式
    virtual void calculateInternalForceVectorUL(std::vector<double>& F);  // Updated Lagrangian格式

    // ========== 能量计算接口 ==========
    virtual double calculateStrainEnergy() override;
    virtual double calculateStrainEnergyTL();  // TL格式应变能
    virtual double calculateStrainEnergyUL();  // UL格式应变能
    virtual double calculateKineticEnergy() override;

    // ========== B矩阵相关接口 ==========
    virtual void computeBMatrix(const NaturalCoord& xi, Matrix& B) const override = 0;
    virtual void computeLinearBMatrix(const NaturalCoord& xi, Matrix& B_L) const;
    virtual void computeNonlinearBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_NL) const;
    virtual void computeTotalLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_TL) const;
    virtual void computeUpdatedLagrangianBMatrix(const NaturalCoord& xi, const Matrix3d& F, Matrix& B_UL) const;

    // ========== 几何和积分接口 ==========
    virtual void computeJacobian(const NaturalCoord& xi, Matrix3d& J) const;
    virtual void computeInverseJacobian(const NaturalCoord& xi, Matrix3d& invJ) const;
    virtual double computeJacobianDeterminant(const NaturalCoord& xi) const;

    virtual double computeReferenceVolume(int gaussPointIndex);
    virtual double computeCurrentVolume(int gaussPointIndex, const Matrix3d& F);

    // ========== 材料点管理接口 ==========
    virtual RgMaterialPoint& getMaterialPoint(int gaussPointIndex);
    virtual const RgMaterialPoint& getMaterialPoint(int gaussPointIndex) const;
    virtual void updateMaterialPointState(int gaussPointIndex, const RgMaterialPoint& mp);
    virtual int getNumberOfIntegrationPoints() const;

    // ========== 非线性格式设置 ==========
    enum NonlinearFormulation
    {
        TOTAL_LAGRANGIAN,    // TL格式
        UPDATED_LAGRANGIAN,  // UL格式
        COROTATIONAL         // 共旋格式
    };

    virtual void setNonlinearFormulation(NonlinearFormulation formulation);
    virtual NonlinearFormulation getNonlinearFormulation() const;

    // ========== 应变度量设置 ==========
    enum StrainMeasure
    {
        GREEN_LAGRANGE,  // 格林-拉格朗日应变
        LOGARITHMIC,     // 对数应变
        BIOT,            // Biot应变
        ALMANSI          // Almansi应变
    };

    virtual void setStrainMeasure(StrainMeasure measure);
    virtual StrainMeasure getStrainMeasure() const;

    // ========== 调试和验证接口 ==========
    virtual bool checkPatchTest() const;
    virtual bool checkEnergyConsistency() const;
    virtual void printState(int gaussPointIndex = -1) const;  // -1表示打印所有积分点

protected:
    // ========== 内部辅助方法 ==========
    virtual void assembleInternalForceFromPK2(const Matrix& B_TL, const StressTensor& S, double weight, double detJ0,
                                              std::vector<double>& F_elem);
    virtual void assembleInternalForceFromCauchy(const Matrix& B_UL, const StressTensor& sigma, double weight,
                                                 double detJ, std::vector<double>& F_elem);

    virtual void computeShapeFunctionDerivatives(const NaturalCoord& xi, std::vector<Vector3d>& dNdX) const;
    virtual void computeCurrentShapeFunctionDerivatives(const NaturalCoord& xi, const Matrix3d& F,
                                                        std::vector<Vector3d>& dNdx) const;

    // ========== 应力推前拉回转换 ==========
    virtual void pushForwardPK2toCauchy(const Matrix3d& F, const StressTensor& S, StressTensor& sigma);
    virtual void pullBackCauchytoPK2(const Matrix3d& F, const StressTensor& sigma, StressTensor& S);

private:
    // ========== 数据成员 ==========
    NonlinearFormulation m_formulation;
    StrainMeasure m_strainMeasure;

    // 状态变量（建议使用）
    std::vector<RgMaterialPoint> m_materialPoints;
    std::vector<Matrix3d> m_deformationGradients;
    // ... 其他必要的数据成员
};

*/