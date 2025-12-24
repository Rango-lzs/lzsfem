#include "RgSolidElement.h"
#include "femcore/RgElementTraitsStore.h"
#include "femcore/FEMesh.h"
#include "femcore/FENode.h"
#include "datastructure/Vector3d.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/tens4d.h"
#include "femcore/Matrix/DenseMatrix.h"
#include "elements/NaturalCoord.h"
#include <vector>
#include <stdexcept>
#include "../ElementTraits/RgSolidElementTraits.h"

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
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            return solidTraits->gaussPoints[n];
        }
    }
    return RgFem::RgGaussPoint(); // Return default if traits not found
}

// --- assembly interfaces (common to all solid elements) ---
void RgSolidElement::calculateStiffnessMatrix(Matrix& K) const
{
    // 默认实现，具体单元类型需要重写此方法
    K.zero();
}

void RgSolidElement::calculateMassMatrix(Matrix& M) const
{
    // 默认实现，具体单元类型需要重写此方法
    M.zero();
}

void RgSolidElement::calculateDampingMatrix(Matrix& C) const
{
    // 默认实现，具体单元类型需要重写此方法
    C.zero();
}

void RgSolidElement::calculateInternalForceVector(Vector& F) const
{
    // 默认实现，具体单元类型需要重写此方法
    //F.zero();
}

void RgSolidElement::calculateStress(FEMaterialPoint& matPt, StressTensor& stress)
{
    // 默认实现，具体单元类型需要重写此方法
    stress.zero();
}

void RgSolidElement::calculateStrain(FEMaterialPoint& matPt, StrainTensor& strain)
{
    // 默认实现，具体单元类型需要重写此方法
    strain.zero();
}

// Protected helper methods
void RgSolidElement::computeBMatrixAtGauss(int gp, Matrix& B) const
{
    // 默认实现，具体单元类型需要重写此方法
    B.zero();
}

void RgSolidElement::getConstitutiveTangentAtGauss(int gp, Matrix& D) const
{
    // 默认实现，具体单元类型需要重写此方法
    D.zero();
}

Matrix3d RgSolidElement::evaluateJacobian(const Vector3d& naturalCoord) const
{
    // 默认实现，具体单元类型需要重写此方法
    return Matrix3d();
}

double RgSolidElement::evaluateJacobianDeterminant(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.det();
}

Matrix3d RgSolidElement::evaluateJacobianInverse(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.inverse();
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

std::vector<double> RgSolidElement::evalH(int n)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits && n < traits->m_nint) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            std::vector<double> result(solidTraits->m_neln);
            for (int i = 0; i < solidTraits->m_neln; ++i) {
                result[i] = solidTraits->m_H[n][i];
            }
            return result;
        }
    }
    return std::vector<double>(); // Return empty vector if traits not found or invalid index
}

std::vector<std::vector<double>> RgSolidElement::evalDeriv(int n)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits && n < traits->m_nint) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            std::vector<std::vector<double>> result(3, std::vector<double>(solidTraits->m_neln));
            for (int i = 0; i < solidTraits->m_neln; ++i) {
                result[0][i] = solidTraits->m_Gr[n][i];  // dN/dr
                result[1][i] = solidTraits->m_Gs[n][i];  // dN/ds
                result[2][i] = solidTraits->m_Gt[n][i];  // dN/dt
            }
            return result;
        }
    }
    return std::vector<std::vector<double>>(); // Return empty vector if traits not found or invalid index
}

std::vector<std::vector<double>> RgSolidElement::evalDeriv2(int n)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits && n < traits->m_nint) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            std::vector<std::vector<double>> result(6, std::vector<double>(solidTraits->m_neln));
            for (int i = 0; i < solidTraits->m_neln; ++i) {
                result[0][i] = solidTraits->Grr[n][i];  // d²N/dr²
                result[1][i] = solidTraits->Gss[n][i];  // d²N/ds²
                result[2][i] = solidTraits->Gtt[n][i];  // d²N/dt²
                result[3][i] = solidTraits->Grs[n][i];  // d²N/drds
                result[4][i] = solidTraits->Gst[n][i];  // d²N/dsdt
                result[5][i] = solidTraits->Grt[n][i];  // d²N/drdt
            }
            return result;
        }
    }
    return std::vector<std::vector<double>>(); // Return empty vector if traits not found or invalid index
}

std::vector<double> RgSolidElement::evalH(const RgFem::NaturalCoord& coord)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            // Use the virtual method from RgSolidElementTraits
            return solidTraits->evalH(coord);
        }
    }
    return std::vector<double>(); // Return empty vector if traits not found
}

std::vector<std::vector<double>> RgSolidElement::evalDeriv(const RgFem::NaturalCoord& coord)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            // Use the virtual method from RgSolidElementTraits
            return solidTraits->evalDeriv(coord);
        }
    }
    return std::vector<std::vector<double>>(); // Return empty vector if traits not found
}

std::vector<std::vector<double>> RgSolidElement::evalDeriv2(const RgFem::NaturalCoord& coord)
{
    // Get the element traits for this element
    RgElementTraits* traits = RgElementTraitsStore::GetInstance()->GetElementTraits(elementType());
    if (traits) {
        // Cast to the appropriate traits type
        RgSolidElementTraits* solidTraits = dynamic_cast<RgSolidElementTraits*>(traits);
        if (solidTraits) {
            // Use the virtual method from RgSolidElementTraits
            return solidTraits->evalDeriv2(coord);
        }
    }
    return std::vector<std::vector<double>>(); // Return empty vector if traits not found
}

void RgSolidElement::Serialize(DumpStream& ar)
{
    // Default implementation - to be overridden by derived classes
}