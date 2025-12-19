#include "RgSolidElement.h"
#include "materials/FEMaterialPoint.h"

RgSolidElement::RgSolidElement(const RgSolidElement& el) 
    : RgElement(el)
{
}

RgSolidElement& RgSolidElement::operator=(const RgSolidElement& el)
{
    if (this != &el) {
        RgElement::operator=(el);
    }
    return *this;
}

int RgSolidElement::dim()
{
    // 默认返回3，因为这是三维实体单元
    return 3;
}

RgFem::RgGaussPoint RgSolidElement::gaussPoint(int n) const
{
    static_cast<RgSolidElementTraits*>(m_pTraits)->gaussPoints(n);
}

// --- assembly interfaces (common to all solid elements) ---
void RgSolidElement::calculateStiffnessMatrix(Matrix& K) const
{
    // 默认实现，具体单元类型需要重写此方法
    K.setZero();
}

void RgSolidElement::calculateTangentStiffnessMatrix(Matrix& Kt) const
{
    // 默认实现，具体单元类型需要重写此方法
    Kt.setZero();
}

void RgSolidElement::calculateMassMatrix(Matrix& M) const
{
    // 默认实现，具体单元类型需要重写此方法
    M.setZero();
}

void RgSolidElement::calculateDampingMatrix(Matrix& C) const
{
    // 默认实现，具体单元类型需要重写此方法
    C.setZero();
}

void RgSolidElement::calculateInternalForceVector(Vector& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    F.setZero();
}

// --- geometric / nonlinear hooks ---
void RgSolidElement::computeGeometricStiffness(Matrix& Kg) const
{
    // 默认实现：零几何刚度矩阵
    Kg.setZero();
}

void RgSolidElement::computeDeformationGradient(int gp, Matrix3d& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    F = Matrix3d(); // Zero matrix
}

// --- integration / material point ---
FEMaterialPoint* RgSolidElement::getMaterialPoint(int gp)
{
    // 默认实现，具体单元类型需要重写此方法
    return nullptr;
}

const FEMaterialPoint* RgSolidElement::getMaterialPoint(int gp) const
{
    // 默认实现，具体单元类型需要重写此方法
    return nullptr;
}

// --- utilities ---
double RgSolidElement::elementVolume() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0.0;
}

void RgSolidElement::updateState(double dt)
{
    // override: 更新材料点状态，具体单元类型需要重写此方法
}

void RgSolidElement::commitState()
{
    // override: 提交当前状态至下一时间步，具体单元类型需要重写此方法
}

void RgSolidElement::resetState()
{
    // override: 重置状态到初始条件，具体单元类型需要重写此方法
}

// Protected helper methods
void RgSolidElement::computeBMatrixAtGauss(int gp, Matrix& B) const
{
    // 默认实现，具体单元类型需要重写此方法
    B.setZero();
}

void RgSolidElement::getConstitutiveTangentAtGauss(int gp, Matrix& D) const
{
    // 默认实现，具体单元类型需要重写此方法
    D.setZero();
}

int RgSolidElement::getNumberOfFaces() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0;
}

int RgSolidElement::getNumberOfEdges() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0;
}

void RgSolidElement::getFaceNodeIds(int faceId, std::vector<int>& faceNodes) const
{
    // 默认实现，具体单元类型需要重写此方法
    faceNodes.clear();
}

void RgSolidElement::getEdgeNodeIds(int edgeId, std::vector<int>& edgeNodes) const
{
    // 默认实现，具体单元类型需要重写此方法
    edgeNodes.clear();
}

Vector3d RgSolidElement::evaluateCoordinates(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return Vector3d(0, 0, 0);
}

Matrix3d RgSolidElement::evaluateJacobian(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return Matrix3d();
}

double RgSolidElement::evaluateJacobianDeterminant(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0.0;
}

Matrix3d RgSolidElement::evaluateJacobianInverse(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return Matrix3d();
}

bool RgSolidElement::isValidNaturalCoordinate(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return false;
}

double RgSolidElement::getCharacteristicLength() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0.0;
}

void RgSolidElement::applyBodyForce(const Vector3d& force, Vector& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    // F向量大小应该在派生类中调整
}

void RgSolidElement::applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    // F向量大小应该在派生类中调整
}

void RgSolidElement::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    // F向量大小应该在派生类中调整
}

double RgSolidElement::calculateStrainEnergy() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0.0;
}

double RgSolidElement::calculateKineticEnergy() const
{
    // 默认实现，具体单元类型需要重写此方法
    return 0.0;
}