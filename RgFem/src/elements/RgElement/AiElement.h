// 层次1: 基础元素接口
class RgElement {
public:
    virtual ~RgElement() = default;
    virtual int getNumberOfNodes() const = 0;
    virtual int getNumberOfGaussPoints() const = 0;
    // 通用接口
};

// 层次2: 实体元素基类
class RgSolidElement : public RgElement {
public:
    virtual int dim() const = 0;  // 空间维度
};

// 层次3: 3D实体元素 - 定义几何和拓扑
class RgSolid3dElement : public RgSolidElement {
public:
    int dim() const override { return 3; }
    
    // 几何相关的共同接口（线性和非线性都需要）
    virtual void evaluateShapeFunctions(double r, double s, double t,
                                       std::vector<double>& N) const = 0;
    virtual void evaluateShapeDerivatives(double r, double s, double t,
                                         std::vector<double>& dNdr,
                                         std::vector<double>& dNds,
                                         std::vector<double>& dNdt) const = 0;
    virtual Matrix3d evaluateJacobian(const Vector3d& natCoord) const = 0;
    
    // 物理相关的共同接口
    virtual void calculateMassMatrix(Matrix& M) const = 0;
    virtual void applyBodyForce(const Vector3d& force, Vector& F) const = 0;
};

// 层次4a: 线性行为 - Mixin类
class LinearBehavior {
public:
    virtual ~LinearBehavior() = default;
    
    // 小变形线性分析接口
    virtual void calculateLinearStiffness(Matrix& K) const = 0;
    virtual void computeLinearStrain(const Vector& u, 
                                     std::vector<Matrix3ds>& strain) const = 0;
    virtual void computeLinearStress(const std::vector<Matrix3ds>& strain,
                                     std::vector<Matrix3ds>& stress) const = 0;
};

// 层次4b: 非线性行为 - Mixin类
class NonlinearBehavior {
public:
    virtual ~NonlinearBehavior() = default;
    
    // 大变形非线性分析接口
    virtual void calculateTangentStiffness(Matrix& Kt) const = 0;
    virtual void calculateGeometricStiffness(Matrix& Kg) const = 0;
    virtual void computeDeformationGradient(int gpId, const Vector& u,
                                           Matrix3d& F) const = 0;
    virtual void computeGreenLagrangeStrain(const Matrix3d& F,
                                           Matrix3ds& E) const = 0;
    virtual void updateNonlinearState(const Vector& u) = 0;
};

// 层次5: 具体拓扑 + 行为组合

// Hex8元素的几何和拓扑（线性和非线性共享）
class RgHex8Topology : public RgSolid3dElement {
public:
    static constexpr int kNodes = 8;
    static constexpr int kGaussPoints = 8;
    
    int getNumberOfNodes() const override { return kNodes; }
    int getNumberOfGaussPoints() const override { return kGaussPoints; }
    
    // 实现形函数（线性和非线性都用同一套）
    void evaluateShapeFunctions(double r, double s, double t,
                               std::vector<double>& N) const override {
        N.resize(kNodes);
        N[0] = 0.125 * (1-r) * (1-s) * (1-t);
        N[1] = 0.125 * (1+r) * (1-s) * (1-t);
        // ... 其他节点
    }
    
    void evaluateShapeDerivatives(double r, double s, double t,
                                 std::vector<double>& dNdr,
                                 std::vector<double>& dNds,
                                 std::vector<double>& dNdt) const override {
        dNdr.resize(kNodes);
        dNds.resize(kNodes);
        dNdt.resize(kNodes);
        
        dNdr[0] = -0.125 * (1-s) * (1-t);
        dNds[0] = -0.125 * (1-r) * (1-t);
        dNdt[0] = -0.125 * (1-r) * (1-s);
        // ... 其他节点
    }
    
    Matrix3d evaluateJacobian(const Vector3d& natCoord) const override {
        std::vector<double> dNdr, dNds, dNdt;
        evaluateShapeDerivatives(natCoord.x, natCoord.y, natCoord.z,
                                dNdr, dNds, dNdt);
        
        Matrix3d J;
        J.zero();
        for (int i = 0; i < kNodes; ++i) {
            Vector3d pos = getNode(i)->getPosition();
            J.m[0][0] += dNdr[i] * pos.x;
            J.m[0][1] += dNdr[i] * pos.y;
            J.m[0][2] += dNdr[i] * pos.z;
            // ... 其他分量
        }
        return J;
    }
    
    // 质量矩阵（线性和非线性都一样）
    void calculateMassMatrix(Matrix& M) const override {
        int ndofs = kNodes * 3;
        M.resize(ndofs, ndofs);
        M.zero();
        
        double rho = getMaterial()->getDensity();
        
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Vector3d natCoord = getGaussPoint(gp);
            std::vector<double> N;
            evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
            
            double detJ = evaluateJacobian(natCoord).det();
            double weight = getGaussWeight(gp) * detJ * rho;
            
            for (int i = 0; i < kNodes; ++i) {
                for (int j = 0; j < kNodes; ++j) {
                    double Mij = weight * N[i] * N[j];
                    for (int d = 0; d < 3; ++d) {
                        M(3*i+d, 3*j+d) += Mij;
                    }
                }
            }
        }
    }
    
protected:
    // 辅助函数供派生类使用
    void computeBMatrix(const Vector3d& natCoord, Matrix& B) const {
        // B矩阵计算（线性和非线性都需要，但使用方式不同）
    }
};

// 最终具体元素：Hex8 + 线性行为
class RgHex8Linear : public RgHex8Topology, public LinearBehavior {
public:
    void calculateLinearStiffness(Matrix& K) const override {
        int ndofs = kNodes * 3;
        K.resize(ndofs, ndofs);
        K.zero();
        
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Vector3d natCoord = getGaussPoint(gp);
            
            // 使用基类的几何计算
            Matrix B;
            computeBMatrix(natCoord, B);
            
            Matrix D = getMaterial()->getElasticityMatrix();
            double detJ = evaluateJacobian(natCoord).det();
            double weight = getGaussWeight(gp) * detJ;
            
            // K += B^T * D * B * weight
            K.addMatrixTripleProduct(B, D, weight);
        }
    }
    
    void computeLinearStrain(const Vector& u,
                            std::vector<Matrix3ds>& strain) const override {
        strain.resize(kGaussPoints);
        
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Vector3d natCoord = getGaussPoint(gp);
            Matrix B;
            computeBMatrix(natCoord, B);
            
            // ε = B * u (小变形线性关系)
            Vector strainVec = B * u;
            strain[gp] = Matrix3ds::fromVoigt(strainVec);
        }
    }
    
    void computeLinearStress(const std::vector<Matrix3ds>& strain,
                            std::vector<Matrix3ds>& stress) const override {
        stress.resize(kGaussPoints);
        Matrix D = getMaterial()->getElasticityMatrix();
        
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Vector strainVec = strain[gp].toVoigt();
            Vector stressVec = D * strainVec;
            stress[gp] = Matrix3ds::fromVoigt(stressVec);
        }
    }
};

// 最终具体元素：Hex8 + 非线性行为
class RgHex8Nonlinear : public RgHex8Topology, public NonlinearBehavior {
public:
    void calculateTangentStiffness(Matrix& Kt) const override {
        int ndofs = kNodes * 3;
        Kt.resize(ndofs, ndofs);
        Kt.zero();
        
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Vector3d natCoord = getGaussPoint(gp);
            
            // 材料刚度
            Matrix Km;
            computeMaterialStiffness(gp, natCoord, Km);
            
            // 几何刚度
            Matrix Kg;
            computeGeometricStiffness(gp, natCoord, Kg);
            
            double detJ = evaluateJacobian(natCoord).det();
            double weight = getGaussWeight(gp) * detJ;
            
            Kt += (Km + Kg) * weight;
        }
    }
    
    void calculateGeometricStiffness(Matrix& Kg) const override {
        // 几何刚度矩阵
    }
    
    void computeDeformationGradient(int gpId, const Vector& u,
                                   Matrix3d& F) const override {
        Vector3d natCoord = getGaussPoint(gpId);
        
        std::vector<double> dNdr, dNds, dNdt;
        evaluateShapeDerivatives(natCoord.x, natCoord.y, natCoord.z,
                                dNdr, dNds, dNdt);
        
        Matrix3d Jinv = evaluateJacobian(natCoord).inverse();
        
        // 计算位移梯度
        Matrix3d gradU;
        gradU.zero();
        
        for (int i = 0; i < kNodes; ++i) {
            // dN/dx = dN/dr * dr/dx
            Vector3d dNdx(
                dNdr[i] * Jinv.m[0][0] + dNds[i] * Jinv.m[1][0] + dNdt[i] * Jinv.m[2][0],
                dNdr[i] * Jinv.m[0][1] + dNds[i] * Jinv.m[1][1] + dNdt[i] * Jinv.m[2][1],
                dNdr[i] * Jinv.m[0][2] + dNds[i] * Jinv.m[1][2] + dNdt[i] * Jinv.m[2][2]
            );
            
            Vector3d ui(u[3*i], u[3*i+1], u[3*i+2]);
            
            // ∂u/∂x
            gradU.m[0][0] += ui.x * dNdx.x;
            gradU.m[0][1] += ui.x * dNdx.y;
            gradU.m[0][2] += ui.x * dNdx.z;
            // ... 其他分量
        }
        
        // F = I + ∂u/∂X
        F = Matrix3d::identity() + gradU;
    }
    
    void computeGreenLagrangeStrain(const Matrix3d& F,
                                   Matrix3ds& E) const override {
        // E = 0.5 * (F^T * F - I)
        Matrix3d C = F.transpose() * F;
        Matrix3d I = Matrix3d::identity();
        E = 0.5 * (C - I);
    }
    
    void updateNonlinearState(const Vector& u) override {
        // 更新所有高斯点的状态
        for (int gp = 0; gp < kGaussPoints; ++gp) {
            Matrix3d F;
            computeDeformationGradient(gp, u, F);
            
            Matrix3ds E;
            computeGreenLagrangeStrain(F, E);
            
            // 存储或更新材料点状态
            auto* mp = getMaterialPoint(gp);
            mp->setDeformationGradient(F);
            mp->setStrain(E);
        }
    }
    
private:
    void computeMaterialStiffness(int gp, const Vector3d& natCoord,
                                 Matrix& Km) const {
        // 非线性材料刚度
    }
    
    void computeGeometricStiffness(int gp, const Vector3d& natCoord,
                                  Matrix& Kg) const {
        // 考虑应力的几何刚度
    }
};