#include "RgBeam2dElement.h"
#include "elements/ElementShape/RgElementShape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <array>



// ============================================================================
// Constructor and Destructor
// ============================================================================

RgBeam2dElement::RgBeam2dElement()
    : RgBeamElement()
    , m_I(1.0)
    , m_A(1.0)
    , m_shearFactor(5.0 / 6.0)  // Default for rectangular cross-section
{
    initTraits();
}

RgBeam2dElement::RgBeam2dElement(const std::array<int, kNodeCount>& nodeIds)
    : RgStructureElement()
    , m_I(1.0)
    , m_A(1.0)
    , m_shearFactor(5.0 / 6.0)
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
    initTraits();
}

RgBeam2dElement::RgBeam2dElement(const RgBeam2dElement& other)
    : RgStructureElement(other)
    , m_I(other.m_I)
    , m_A(other.m_A)
    , m_shearFactor(other.m_shearFactor)
{
    m_gaussR = other.m_gaussR;
    m_gaussW = other.m_gaussW;
}

RgBeam2dElement::~RgBeam2dElement()
{
}

RgBeam2dElement& RgBeam2dElement::operator=(const RgBeam2dElement& other)
{
    if (this != &other) {
        RgBeamElement::operator=(other);
        m_I = other.m_I;
        m_A = other.m_A;
        m_shearFactor = other.m_shearFactor;
        m_gaussR = other.m_gaussR;
        m_gaussW = other.m_gaussW;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgBeam2dElement::elementType() const
{
    return ElementType::FE_BEAM2DG2;
}

ElementShape RgBeam2dElement::elementShape() const
{
    return ElementShape::ET_BEAM2D;
}

ElementCategory RgBeam2dElement::elementClass() const
{
    return ElementCategory::FE_ELEM_BEAM;
}

// ============================================================================
// Element Properties & Initialization
// ============================================================================

void RgBeam2dElement::initTraits()
{
    RgSolidElementTraits* traits = new RgSolidElementTraits(kNumGaussPoints, kNodeCount, ET_BEAM2D, FE_BEAM2DG2);
    traits->init();
    
    if (m_pTraits == nullptr) {
        m_pTraits = traits;
        m_node.resize(kNodeCount);
        m_loc_node.resize(kNodeCount);
    }
    
    initializeGaussPoints();
}

void RgBeam2dElement::initializeGaussPoints()
{
    m_gaussR.resize(kNumGaussPoints);
    m_gaussW.resize(kNumGaussPoints);
    
    // 2-point Gauss quadrature
    const double a = 1.0 / std::sqrt(3.0);  // ±0.5773502691896...
    
    m_gaussR[0] = -a;
    m_gaussR[1] = a;
    m_gaussW[0] = 1.0;
    m_gaussW[1] = 1.0;
}

// ============================================================================
// Geometric Operations
// ============================================================================

Vector3d RgBeam2dElement::evaluateCoordinates(double naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord, N);
    
    Vector3d coords(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            coords += N[i] * node->getPosition();
        }
    }
    return coords;
}

double RgBeam2dElement::evaluateLength() const
{
    FENode* n0 = getNode(0);
    FENode* n1 = getNode(1);
    
    if (!n0 || !n1) return 0.0;
    
    Vector3d v = n1->getPosition() - n0->getPosition();
    return v.Length();
}

Vector3d RgBeam2dElement::getAxis() const
{
    FENode* n0 = getNode(0);
    FENode* n1 = getNode(1);
    
    if (!n0 || !n1) return Vector3d(1, 0, 0);
    
    Vector3d v = n1->getPosition() - n0->getPosition();
    double len = v.Length();
    if (len > 1e-10) {
        v.Normalize();
    }
    return v;
}

// ============================================================================
// Shape Function Evaluations
// ============================================================================

void RgBeam2dElement::evaluateShapeFunctions(double r, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    
    // Linear shape functions
    // N0 = (1 - r) / 2
    // N1 = (1 + r) / 2
    N[0] = (1.0 - r) / 2.0;
    N[1] = (1.0 + r) / 2.0;
}

void RgBeam2dElement::evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const
{
    dN_dr.assign(kNodeCount, 0.0);
    
    // dN0/dr = -1/2
    // dN1/dr = 1/2
    dN_dr[0] = -0.5;
    dN_dr[1] = 0.5;
}

// ============================================================================
// Transformation Matrix
// ============================================================================

Matrix3d RgBeam2dElement::computeTransformationMatrix() const
{
    Vector3d axis = getAxis();
    
    // Transformation from local to global coordinates
    // Local x-axis along beam axis
    Vector3d ex = axis;
    
    // Local y-axis perpendicular to beam axis (in xy-plane)
    Vector3d ey(0, 1, 0);
    if (std::fabs(ex.y) < 0.99) {
        ey = (Vector3d(0, 1, 0) - ex.y * ex).Normalized();
    } else {
        ey = (Vector3d(1, 0, 0) - ex.x * ex).Normalized();
    }
    
    // Local z-axis
    Vector3d ez = ex.Cross(ey);
    
    Matrix3d T;
    T.m[0][0] = ex.x; T.m[0][1] = ex.y; T.m[0][2] = ex.z;
    T.m[1][0] = ey.x; T.m[1][1] = ey.y; T.m[1][2] = ey.z;
    T.m[2][0] = ez.x; T.m[2][1] = ez.y; T.m[2][2] = ez.z;
    
    return T;
}

void RgBeam2dElement::computeLocalStiffness(Matrix& K_local) const
{
    int ndofs = kTotalDofs;
    K_local.resize(ndofs, ndofs);
    K_local.zero();
    
    double L = evaluateLength();
    if (L < 1e-10) return;
    
    // Get material properties
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    double E = matPt->m_pMat->getElasticModulus();
    double G = matPt->m_pMat->getShearModulus();
    
    // Timoshenko beam stiffness matrix (local coordinates)
    double EI = E * m_I;
    double EA = E * m_A;
    double GAeff = G * m_shearFactor * m_A;
    double phi = 12.0 * EI / (GAeff * L * L);  // Shear parameter
    
    // Element stiffness in local frame (6x6 matrix, 3 DOFs per node)
    // DOFs: [u1, v1, r1, u2, v2, r2]
    
    double k1 = EA / L;
    double k2 = 12.0 * EI / (L * L * L * (1.0 + phi));
    double k3 = 6.0 * EI / (L * L * (1.0 + phi));
    double k4 = (4.0 + phi) * EI / (L * (1.0 + phi));
    double k5 = (2.0 - phi) * EI / (L * (1.0 + phi));
    
    // Axial stiffness (DOFs 0 and 3)
    K_local[0][0] = k1;
    K_local[3][3] = k1;
    K_local[0][3] = -k1;
    K_local[3][0] = -k1;
    
    // Bending stiffness
    K_local[1][1] = k2;
    K_local[2][2] = k4;
    K_local[4][4] = k2;
    K_local[5][5] = k4;
    
    K_local[1][2] = k3;
    K_local[2][1] = k3;
    K_local[1][4] = -k2;
    K_local[4][1] = -k2;
    K_local[1][5] = k3;
    K_local[5][1] = k3;
    K_local[2][4] = -k3;
    K_local[4][2] = -k3;
    K_local[2][5] = k5;
    K_local[5][2] = k5;
    K_local[4][5] = -k3;
    K_local[5][4] = -k3;
}

// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgBeam2dElement::calculateStiffnessMatrix(Matrix& K) const
{
    int ndofs = kTotalDofs;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Compute local stiffness matrix
    Matrix K_local;
    computeLocalStiffness(K_local);
    
    // Compute transformation matrix
    Matrix3d T = computeTransformationMatrix();
    
    // For 2D beam with 3 DOFs per node, create full transformation
    // T_full is 6x6 block diagonal with T_2d blocks
    Matrix T_full(ndofs, ndofs);
    T_full.zero();
    
    // First node transformation (2D)
    T_full[0][0] = T.m[0][0]; T_full[0][1] = T.m[0][1];
    T_full[1][0] = T.m[1][0]; T_full[1][1] = T.m[1][1];
    T_full[2][2] = 1.0;  // Rotation about z
    
    // Second node transformation (2D)
    T_full[3][3] = T.m[0][0]; T_full[3][4] = T.m[0][1];
    T_full[4][3] = T.m[1][0]; T_full[4][4] = T.m[1][1];
    T_full[5][5] = 1.0;  // Rotation about z
    
    // Transform to global: K = T_full^T * K_local * T_full
    Matrix KT(ndofs, ndofs);
    KT.zero();
    
    for (int i = 0; i < ndofs; ++i) {
        for (int j = 0; j < ndofs; ++j) {
            for (int k = 0; k < ndofs; ++k) {
                KT[i][j] += K_local[i][k] * T_full[k][j];
            }
        }
    }
    
    for (int i = 0; i < ndofs; ++i) {
        for (int j = 0; j < ndofs; ++j) {
            for (int k = 0; k < ndofs; ++k) {
                K[i][j] += T_full[k][i] * KT[k][j];
            }
        }
    }
}

void RgBeam2dElement::calculateMassMatrix(Matrix& M) const
{
    int ndofs = kTotalDofs;
    M.resize(ndofs, ndofs);
    M.zero();
    
    double L = evaluateLength();
    if (L < 1e-10) return;
    
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    double rho = matPt->m_pMat->getDensity();
    
    // Consistent mass matrix for linear beam
    // Local coordinates
    double m = rho * m_A * L;  // Total mass
    double I_rot = rho * m_I * L;  // Rotational inertia
    
    // Lumped mass matrix (simplified)
    M[0][0] = m / 2.0;  // Node 0, axial
    M[1][1] = m / 2.0;  // Node 0, transverse
    M[2][2] = I_rot / 12.0;  // Node 0, rotation
    M[3][3] = m / 2.0;  // Node 1, axial
    M[4][4] = m / 2.0;  // Node 1, transverse
    M[5][5] = I_rot / 12.0;  // Node 1, rotation
}

void RgBeam2dElement::calculateDampingMatrix(Matrix& C) const
{
    int ndofs = kTotalDofs;
    C.resize(ndofs, ndofs);
    C.zero();
    
    // Rayleigh damping: C = a0*M + a1*K
    Matrix M, K;
    calculateMassMatrix(M);
    calculateStiffnessMatrix(K);
    
    // Default Rayleigh coefficients
    double a0 = 0.0;
    double a1 = 0.0;
    
    for (int i = 0; i < ndofs; ++i) {
        for (int j = 0; j < ndofs; ++j) {
            C[i][j] = a0 * M[i][j] + a1 * K[i][j];
        }
    }
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgBeam2dElement::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    if (!matPt.m_pMat) return;
    
    // Axial stress: sigma = E * strain
    // where strain is computed from displacement gradient
}

void RgBeam2dElement::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    strain.zero();
    
    // Small strain: ε = du/dx
    // Computed from displacement field
}

// ============================================================================
// Loading Operations
// ============================================================================

void RgBeam2dElement::applyBodyForce(const Vector3d& force, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    // Distributed body force over element length
    double L = evaluateLength();
    double f_mag = force.Length();
    
    // Distribute force to nodes
    for (int i = 0; i < kNodeCount; ++i) {
        F[3*i]   += force.x * L / 2.0;  // Axial
        F[3*i+1] += force.y * L / 2.0;  // Transverse
    }
}

void RgBeam2dElement::applyDistributedLoad(const Vector3d& load, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    // Distributed load per unit length
    double L = evaluateLength();
    
    // Equivalent nodal forces for distributed load
    // For uniform load w: F = w*L/2 at each end
    F[1] = load.y * L / 2.0;  // Node 0, transverse
    F[4] = load.y * L / 2.0;  // Node 1, transverse
    
    // Moments due to distributed load
    F[2] = load.y * L * L / 12.0;   // Node 0, moment
    F[5] = -load.y * L * L / 12.0;  // Node 1, moment
}

void RgBeam2dElement::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    if (nodeId >= 0 && nodeId < kNodeCount) {
        F[3*nodeId]   = force.x;
        F[3*nodeId+1] = force.y;
    }
}

// ============================================================================
// Material Point Access
// ============================================================================

FEMaterialPoint* RgBeam2dElement::getMaterialPoint(int gaussPtId)
{
    if (gaussPtId >= 0 && gaussPtId < kNumGaussPoints && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

const FEMaterialPoint* RgBeam2dElement::getMaterialPoint(int gaussPtId) const
{
    if (gaussPtId >= 0 && gaussPtId < kNumGaussPoints && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

// ============================================================================
// Serialization
// ============================================================================

void RgBeam2dElement::Serialize(DumpStream& ar)
{
    RgBeamElement::Serialize(ar);
    
    if (ar.IsSaving()) {
        ar << m_I << m_A << m_shearFactor;
        ar << (int)m_gaussR.size();
        for (double v : m_gaussR) ar << v;
        for (double v : m_gaussW) ar << v;
    } else {
        ar >> m_I >> m_A >> m_shearFactor;
        int size = 0;
        ar >> size;
        m_gaussR.resize(size);
        m_gaussW.resize(size);
        for (int i = 0; i < size; ++i) {
            ar >> m_gaussR[i] >> m_gaussW[i];
        }
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

double RgBeam2dElement::getCharacteristicLength() const
{
    return evaluateLength();
}

double RgBeam2dElement::getVolume() const
{
    return m_A * evaluateLength();
}


