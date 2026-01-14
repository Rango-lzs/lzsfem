#include "RgBeam3dElement.h"
#include "elements/ElementShape/RgElementShape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/RgMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <array>



// ============================================================================
// Constructor and Destructor
// ============================================================================

RgBeam3dElement::RgBeam3dElement()
    : RgBeamElement()
    , m_Iy(1.0)
    , m_Iz(1.0)
    , m_Ip(1.0)
    , m_A(1.0)
    , m_shearFactor(5.0 / 6.0)
    , m_orient(0, 1, 0)
{
    initTraits();
}

RgBeam3dElement::RgBeam3dElement(const std::array<int, kNodeCount>& nodeIds)
    : RgStructureElement()
    , m_Iy(1.0)
    , m_Iz(1.0)
    , m_Ip(1.0)
    , m_A(1.0)
    , m_shearFactor(5.0 / 6.0)
    , m_orient(0, 1, 0)
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
    initTraits();
}

RgBeam3dElement::RgBeam3dElement(const RgBeam3dElement& other)
    : RgStructureElement(other)
    , m_Iy(other.m_Iy)
    , m_Iz(other.m_Iz)
    , m_Ip(other.m_Ip)
    , m_A(other.m_A)
    , m_shearFactor(other.m_shearFactor)
    , m_orient(other.m_orient)
{
    m_gaussR = other.m_gaussR;
    m_gaussW = other.m_gaussW;
    m_rotation = other.m_rotation;
}

RgBeam3dElement::~RgBeam3dElement()
{
}

RgBeam3dElement& RgBeam3dElement::operator=(const RgBeam3dElement& other)
{
    if (this != &other) {
        RgBeamElement::operator=(other);
        m_Iy = other.m_Iy;
        m_Iz = other.m_Iz;
        m_Ip = other.m_Ip;
        m_A = other.m_A;
        m_shearFactor = other.m_shearFactor;
        m_orient = other.m_orient;
        m_gaussR = other.m_gaussR;
        m_gaussW = other.m_gaussW;
        m_rotation = other.m_rotation;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgBeam3dElement::elementType() const
{
    return ElementType::FE_BEAM3DG2;
}

ElementShape RgBeam3dElement::elementShape() const
{
    return ElementShape::ET_BEAM3D;
}

ElementCategory RgBeam3dElement::elementClass() const
{
    return ElementCategory::FE_ELEM_BEAM;
}

// ============================================================================
// Element Properties & Initialization
// ============================================================================

void RgBeam3dElement::initTraits()
{
    RgSolidElementTraits* traits = new RgSolidElementTraits(kNumGaussPoints, kNodeCount, ET_BEAM3D, FE_BEAM3DG2);
    traits->init();
    
    if (m_pTraits == nullptr) {
        m_pTraits = traits;
        m_node.resize(kNodeCount);
        m_loc_node.resize(kNodeCount);
    }
    
    initializeGaussPoints();
}

void RgBeam3dElement::initializeGaussPoints()
{
    m_gaussR.resize(kNumGaussPoints);
    m_gaussW.resize(kNumGaussPoints);
    m_rotation.resize(kNumGaussPoints, Matrix3d());
    
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

Vector3d RgBeam3dElement::evaluateCoordinates(double naturalCoord) const
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

double RgBeam3dElement::evaluateLength() const
{
    FENode* n0 = getNode(0);
    FENode* n1 = getNode(1);
    
    if (!n0 || !n1) return 0.0;
    
    Vector3d v = n1->getPosition() - n0->getPosition();
    return v.Length();
}

Vector3d RgBeam3dElement::getAxis() const
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

void RgBeam3dElement::evaluateShapeFunctions(double r, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    
    // Linear shape functions
    N[0] = (1.0 - r) / 2.0;
    N[1] = (1.0 + r) / 2.0;
}

void RgBeam3dElement::evaluateShapeDerivatives(double r, std::vector<double>& dN_dr) const
{
    dN_dr.assign(kNodeCount, 0.0);
    
    dN_dr[0] = -0.5;
    dN_dr[1] = 0.5;
}

// ============================================================================
// Transformation Matrix (Local to Global)
// ============================================================================

Matrix3d RgBeam3dElement::computeTransformationMatrix() const
{
    Vector3d e1 = getAxis();  // Beam axis
    
    // Create orthogonal triad e1, e2, e3
    Vector3d e2, e3;
    
    // e2 perpendicular to e1, preferably in direction of m_orient
    Vector3d up = m_orient;
    Vector3d temp = e1.Cross(up);
    
    if (temp.Length() < 0.01) {
        // up is parallel to e1, use different reference
        if (std::fabs(e1.z) < 0.99) {
            up = Vector3d(0, 0, 1);
        } else {
            up = Vector3d(1, 0, 0);
        }
        temp = e1.Cross(up);
    }
    
    e3 = temp.Normalized();  // e3 perpendicular to e1
    e2 = e3.Cross(e1);       // e2 perpendicular to both
    e2.Normalize();
    
    Matrix3d T;
    T.m[0][0] = e1.x; T.m[0][1] = e1.y; T.m[0][2] = e1.z;
    T.m[1][0] = e2.x; T.m[1][1] = e2.y; T.m[1][2] = e2.z;
    T.m[2][0] = e3.x; T.m[2][1] = e3.y; T.m[2][2] = e3.z;
    
    return T;
}

Matrix3d RgBeam3dElement::computeRotationMatrix() const
{
    // For small rotation analysis, rotation matrix ≈ I
    // For geometric nonlinear: implement exponential map or quaternion
    return Matrix3d::Identity();
}

void RgBeam3dElement::updateRotations(const double* displacements)
{
    // Update rotations based on displacement field
    // For large rotation analysis
}

void RgBeam3dElement::computeLocalStiffness(Matrix& K_local) const
{
    int ndofs = kTotalDofs;
    K_local.resize(ndofs, ndofs);
    K_local.zero();
    
    double L = evaluateLength();
    if (L < 1e-10) return;
    
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    double E = matPt->m_pMat->getElasticModulus();
    double G = matPt->m_pMat->getShearModulus();
    
    // Timoshenko beam stiffness in local coordinates
    double EIy = E * m_Iy;
    double EIz = E * m_Iz;
    double GIp = G * m_Ip;
    double EA = E * m_A;
    double GAeff_y = G * m_shearFactor * m_A;  // Shear in y direction
    double GAeff_z = G * m_shearFactor * m_A;  // Shear in z direction
    
    // Shear parameters
    double phi_y = 12.0 * EIz / (GAeff_y * L * L);
    double phi_z = 12.0 * EIy / (GAeff_z * L * L);
    
    // 12x12 local stiffness matrix
    // DOFs: [u1, v1, w1, rx1, ry1, rz1, u2, v2, w2, rx2, ry2, rz2]
    
    // Axial (compression/tension) - terms for DOF 0 and 6
    double k_ax = EA / L;
    K_local[0][0] = k_ax;
    K_local[6][6] = k_ax;
    K_local[0][6] = -k_ax;
    K_local[6][0] = -k_ax;
    
    // Torsion (rotation about x-axis) - terms for DOF 3 and 9
    double k_tors = GIp / L;
    K_local[3][3] = k_tors;
    K_local[9][9] = k_tors;
    K_local[3][9] = -k_tors;
    K_local[9][3] = -k_tors;
    
    // Bending about y-axis and shear in z - terms for DOF 2, 5, 8, 11
    double k2y = 12.0 * EIz / (L * L * L * (1.0 + phi_y));
    double k3y = 6.0 * EIz / (L * L * (1.0 + phi_y));
    double k4y = (4.0 + phi_y) * EIz / (L * (1.0 + phi_y));
    double k5y = (2.0 - phi_y) * EIz / (L * (1.0 + phi_y));
    
    K_local[2][2] = k2y;
    K_local[5][5] = k4y;
    K_local[8][8] = k2y;
    K_local[11][11] = k4y;
    K_local[2][5] = k3y;
    K_local[5][2] = k3y;
    K_local[2][8] = -k2y;
    K_local[8][2] = -k2y;
    K_local[2][11] = k3y;
    K_local[11][2] = k3y;
    K_local[5][8] = -k3y;
    K_local[8][5] = -k3y;
    K_local[5][11] = k5y;
    K_local[11][5] = k5y;
    K_local[8][11] = -k3y;
    K_local[11][8] = -k3y;
    
    // Bending about z-axis and shear in y - terms for DOF 1, 4, 7, 10
    double k2z = 12.0 * EIy / (L * L * L * (1.0 + phi_z));
    double k3z = 6.0 * EIy / (L * L * (1.0 + phi_z));
    double k4z = (4.0 + phi_z) * EIy / (L * (1.0 + phi_z));
    double k5z = (2.0 - phi_z) * EIy / (L * (1.0 + phi_z));
    
    K_local[1][1] = k2z;
    K_local[4][4] = k4z;
    K_local[7][7] = k2z;
    K_local[10][10] = k4z;
    K_local[1][4] = -k3z;
    K_local[4][1] = -k3z;
    K_local[1][7] = -k2z;
    K_local[7][1] = -k2z;
    K_local[1][10] = -k3z;
    K_local[10][1] = -k3z;
    K_local[4][7] = k3z;
    K_local[7][4] = k3z;
    K_local[4][10] = k5z;
    K_local[10][4] = k5z;
    K_local[7][10] = k3z;
    K_local[10][7] = k3z;
}

// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgBeam3dElement::calculateStiffnessMatrix(Matrix& K) const
{
    int ndofs = kTotalDofs;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Compute local stiffness
    Matrix K_local;
    computeLocalStiffness(K_local);
    
    // Compute transformation matrix (3x3)
    Matrix3d T = computeTransformationMatrix();
    
    // Create full 12x12 transformation matrix (block diagonal)
    Matrix T_full(ndofs, ndofs);
    T_full.zero();
    
    // First node (6 DOFs)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_full[i][j] = T.m[i][j];
            T_full[i+3][j+3] = T.m[i][j];  // Rotation part
        }
    }
    
    // Second node (6 DOFs)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_full[i+6][j+6] = T.m[i][j];
            T_full[i+9][j+9] = T.m[i][j];  // Rotation part
        }
    }
    
    // Transform: K = T_full^T * K_local * T_full
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

void RgBeam3dElement::calculateMassMatrix(Matrix& M) const
{
    int ndofs = kTotalDofs;
    M.resize(ndofs, ndofs);
    M.zero();
    
    double L = evaluateLength();
    if (L < 1e-10) return;
    
    const FEMaterialPoint* matPt = getMaterialPoint(0);
    if (!matPt || !matPt->m_pMat) return;
    
    double rho = matPt->m_pMat->getDensity();
    
    // Lumped mass matrix
    double m = rho * m_A * L / 2.0;  // Half mass per node
    
    // Approximate rotational inertia
    double I_rot = rho * (m_Iy + m_Iz) * L / 12.0;
    
    // Translational mass
    M[0][0] = m;
    M[1][1] = m;
    M[2][2] = m;
    
    M[6][6] = m;
    M[7][7] = m;
    M[8][8] = m;
    
    // Rotational inertia
    M[3][3] = I_rot;
    M[4][4] = I_rot;
    M[5][5] = I_rot;
    
    M[9][9] = I_rot;
    M[10][10] = I_rot;
    M[11][11] = I_rot;
}

void RgBeam3dElement::calculateDampingMatrix(Matrix& C) const
{
    int ndofs = kTotalDofs;
    C.resize(ndofs, ndofs);
    C.zero();
    
    // Rayleigh damping: C = a0*M + a1*K
    Matrix M, K;
    calculateMassMatrix(M);
    calculateStiffnessMatrix(K);
    
    double a0 = 0.0;  // Mass damping coefficient
    double a1 = 0.0;  // Stiffness damping coefficient
    
    for (int i = 0; i < ndofs; ++i) {
        for (int j = 0; j < ndofs; ++j) {
            C[i][j] = a0 * M[i][j] + a1 * K[i][j];
        }
    }
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgBeam3dElement::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    if (!matPt.m_pMat) return;
    
    Matrix3ds strain;
    calculateStrain(matPt, strain);
    matPt.m_pMat->getStress(matPt, strain, stress);
}

void RgBeam3dElement::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    strain.zero();
}

// ============================================================================
// Loading Operations
// ============================================================================

void RgBeam3dElement::applyBodyForce(const Vector3d& force, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    double L = evaluateLength();
    
    // Distribute body force to nodes
    for (int i = 0; i < kNodeCount; ++i) {
        F[6*i]   += force.x * L / 2.0;  // fx
        F[6*i+1] += force.y * L / 2.0;  // fy
        F[6*i+2] += force.z * L / 2.0;  // fz
    }
}

void RgBeam3dElement::applyDistributedLoad(const Vector3d& load, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    double L = evaluateLength();
    
    // Distributed load per unit length
    F[1] = load.y * L / 2.0;   // Node 0
    F[7] = load.y * L / 2.0;   // Node 1
    
    F[2] = load.z * L / 2.0;   // Node 0
    F[8] = load.z * L / 2.0;   // Node 1
}

void RgBeam3dElement::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    if (nodeId >= 0 && nodeId < kNodeCount) {
        F[6*nodeId]   = force.x;
        F[6*nodeId+1] = force.y;
        F[6*nodeId+2] = force.z;
    }
}

void RgBeam3dElement::applyMoment(int nodeId, const Vector3d& moment, Vector& F) const
{
    int ndofs = kTotalDofs;
    F.resize(ndofs);
    F.zero();
    
    if (nodeId >= 0 && nodeId < kNodeCount) {
        F[6*nodeId+3] = moment.x;
        F[6*nodeId+4] = moment.y;
        F[6*nodeId+5] = moment.z;
    }
}

// ============================================================================
// Material Point Access
// ============================================================================

RgMaterialPoint* RgBeam3dElement::getMaterialPoint(int gaussPtId)
{
    if (gaussPtId >= 0 && gaussPtId < kNumGaussPoints && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

const RgMaterialPoint* RgBeam3dElement::getMaterialPoint(int gaussPtId) const
{
    if (gaussPtId >= 0 && gaussPtId < kNumGaussPoints && m_MaterialPoint) {
        return m_MaterialPoint + gaussPtId;
    }
    return nullptr;
}

// ============================================================================
// Serialization
// ============================================================================

void RgBeam3dElement::Serialize(DumpStream& ar)
{
    RgBeamElement::Serialize(ar);
    
    if (ar.IsSaving()) {
        ar << m_Iy << m_Iz << m_Ip << m_A << m_shearFactor;
        ar << m_orient.x << m_orient.y << m_orient.z;
        ar << (int)m_gaussR.size();
        for (double v : m_gaussR) ar << v;
        for (double v : m_gaussW) ar << v;
    } else {
        ar >> m_Iy >> m_Iz >> m_Ip >> m_A >> m_shearFactor;
        ar >> m_orient.x >> m_orient.y >> m_orient.z;
        int size = 0;
        ar >> size;
        m_gaussR.resize(size);
        m_gaussW.resize(size);
        m_rotation.resize(size);
        for (int i = 0; i < size; ++i) {
            ar >> m_gaussR[i] >> m_gaussW[i];
        }
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

double RgBeam3dElement::getCharacteristicLength() const
{
    return evaluateLength();
}

double RgBeam3dElement::getVolume() const
{
    return m_A * evaluateLength();
}


