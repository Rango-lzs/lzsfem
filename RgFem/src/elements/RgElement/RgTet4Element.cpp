#include "RgTet4Element.h"
#include "elements/ElementShape/RgElementShape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>

namespace RgFem {

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgTet4Element::RgTet4Element()
    : RgSolid3dElement()
{
}

RgTet4Element::RgTet4Element(const std::array<int, kNodeCount>& nodeIds)
    : RgSolid3dElement()
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
}

RgTet4Element::RgTet4Element(const RgTet4Element& other)
    : RgSolid3dElement(other)
{
    m_gaussR = other.m_gaussR;
    m_gaussS = other.m_gaussS;
    m_gaussT = other.m_gaussT;
    m_gaussW = other.m_gaussW;
    m_jacobianInverse = other.m_jacobianInverse;
    m_jacobianDet = other.m_jacobianDet;
}

RgTet4Element::~RgTet4Element()
{
}

RgTet4Element& RgTet4Element::operator=(const RgTet4Element& other)
{
    if (this != &other) {
        RgSolid3dElement::operator=(other);
        m_gaussR = other.m_gaussR;
        m_gaussS = other.m_gaussS;
        m_gaussT = other.m_gaussT;
        m_gaussW = other.m_gaussW;
        m_jacobianInverse = other.m_jacobianInverse;
        m_jacobianDet = other.m_jacobianDet;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgTet4Element::elementType() const
{
    return ElementType::FE_TET4G1;
}

ElementShape RgTet4Element::elementShape() const
{
    return ElementShape::ET_TET4;
}

ElementCategory RgTet4Element::elementClass() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

// ============================================================================
// Element Properties
// ============================================================================

int RgTet4Element::getNumberOfGaussPoints() const
{
    return m_pTraits ? m_pTraits->m_nint : 1;
}

// ============================================================================
// Initialize Element Traits
// ============================================================================

void RgTet4Element::initTraits()
{
    RgSolidElementTraits* traits = new RgSolidElementTraits(1, kNodeCount, ET_TET4, FE_TET4G1);
    traits->init();
    
    if (m_pTraits == nullptr) {
        m_pTraits = traits;
        m_node.resize(kNodeCount);
        m_loc_node.resize(kNodeCount);
    }
    
    // Single Gauss point at centroid
    m_gaussR.assign(1, 0.25);
    m_gaussS.assign(1, 0.25);
    m_gaussT.assign(1, 0.25);
    m_gaussW.assign(1, 1.0);
    
    m_jacobianInverse.resize(1);
    m_jacobianDet.resize(1);
}

// ============================================================================
// Shape Function Evaluations
// ============================================================================

void RgTet4Element::evaluateShapeFunctions(double r, double s, double t, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    
    // Linear tetrahedral shape functions
    // Node numbering: 0:(0,0,0), 1:(1,0,0), 2:(0,1,0), 3:(0,0,1)
    N[0] = 1.0 - r - s - t;  // L1
    N[1] = r;                 // L2
    N[2] = s;                 // L3
    N[3] = t;                 // L4
}

void RgTet4Element::evaluateShapeDerivatives(double r, double s, double t,
                                             std::vector<double>& dN_dr,
                                             std::vector<double>& dN_ds,
                                             std::vector<double>& dN_dt) const
{
    dN_dr.assign(kNodeCount, 0.0);
    dN_ds.assign(kNodeCount, 0.0);
    dN_dt.assign(kNodeCount, 0.0);
    
    // dN/dr
    dN_dr[0] = -1.0;
    dN_dr[1] = 1.0;
    dN_dr[2] = 0.0;
    dN_dr[3] = 0.0;
    
    // dN/ds
    dN_ds[0] = -1.0;
    dN_ds[1] = 0.0;
    dN_ds[2] = 1.0;
    dN_ds[3] = 0.0;
    
    // dN/dt
    dN_dt[0] = -1.0;
    dN_dt[1] = 0.0;
    dN_dt[2] = 0.0;
    dN_dt[3] = 1.0;
}

// ============================================================================
// Geometric Operations
// ============================================================================

Vector3d RgTet4Element::evaluateCoordinates(const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d coords(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            coords += N[i] * node->getPosition();
        }
    }
    return coords;
}

Matrix3d RgTet4Element::evaluateJacobian(const Vector3d& naturalCoord) const
{
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d J(0, 0, 0, 0, 0, 0, 0, 0, 0);
    
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            Vector3d pos = node->getPosition();
            J.m[0][0] += dN_dr[i] * pos.x;
            J.m[0][1] += dN_dr[i] * pos.y;
            J.m[0][2] += dN_dr[i] * pos.z;
            J.m[1][0] += dN_ds[i] * pos.x;
            J.m[1][1] += dN_ds[i] * pos.y;
            J.m[1][2] += dN_ds[i] * pos.z;
            J.m[2][0] += dN_dt[i] * pos.x;
            J.m[2][1] += dN_dt[i] * pos.y;
            J.m[2][2] += dN_dt[i] * pos.z;
        }
    }
    return J;
}

double RgTet4Element::evaluateJacobianDeterminant(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.det();
}

Matrix3d RgTet4Element::evaluateJacobianInverse(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.inverse();
}

// ============================================================================
// Physical Field Evaluations
// ============================================================================

Vector3d RgTet4Element::evaluateField(const Vector3d* nodeValues, const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d field(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        field += N[i] * nodeValues[i];
    }
    return field;
}

double RgTet4Element::evaluateScalarField(const double* nodeValues, const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    double field = 0.0;
    for (int i = 0; i < kNodeCount; ++i) {
        field += N[i] * nodeValues[i];
    }
    return field;
}

// ============================================================================
// FEM Matrix Calculations
// ============================================================================

void RgTet4Element::calculateStiffnessMatrix(Matrix& K) const
{
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Integration over single Gauss point
    Vector3d natCoord(m_gaussR[0], m_gaussS[0], m_gaussT[0]);
    double jdet = evaluateJacobianDeterminant(natCoord);
    double weight = m_gaussW[0] * jdet;
    
    // Compute B matrix
    Matrix B;
    computeBMatrix(natCoord, B);
    
    // K = B^T * D * B * weight
    // (Material matrix D would come from material model)
}

void RgTet4Element::calculateMassMatrix(Matrix& M) const
{
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    double rho = 1.0;  // From material
    
    Vector3d natCoord(m_gaussR[0], m_gaussS[0], m_gaussT[0]);
    std::vector<double> N;
    evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
    
    double jdet = evaluateJacobianDeterminant(natCoord);
    double weight = m_gaussW[0] * jdet * rho;
    
    for (int i = 0; i < kNodeCount; ++i) {
        for (int j = 0; j < kNodeCount; ++j) {
            double Mij = weight * N[i] * N[j];
            for (int d = 0; d < 3; ++d) {
                M(3*i + d, 3*j + d) += Mij;
            }
        }
    }
}

void RgTet4Element::calculateDampingMatrix(Matrix& C) const
{
    int ndofs = kNodeCount * 3;
    C.resize(ndofs, ndofs);
    C.zero();
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgTet4Element::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    // Placeholder for stress calculation
}

void RgTet4Element::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    // Placeholder for strain calculation
}

// ============================================================================
// Face and Edge Operations
// ============================================================================

void RgTet4Element::getFaceNodeIds(int faceId, std::array<int, 3>& faceNodes) const
{
    // Four triangular faces
    static const int faces[4][3] = {
        {0, 1, 2},  // Face 0
        {0, 3, 1},  // Face 1
        {1, 3, 2},  // Face 2
        {0, 2, 3}   // Face 3
    };
    
    if (faceId >= 0 && faceId < 4) {
        std::copy(faces[faceId], faces[faceId] + 3, faceNodes.begin());
    }
}

void RgTet4Element::getEdgeNodeIds(int edgeId, std::array<int, 2>& edgeNodes) const
{
    // Six edges
    static const int edges[6][2] = {
        {0, 1}, {1, 2}, {2, 0},  // Base triangle
        {0, 3}, {1, 3}, {2, 3}   // Edges to apex
    };
    
    if (edgeId >= 0 && edgeId < 6) {
        std::copy(edges[edgeId], edges[edgeId] + 2, edgeNodes.begin());
    }
}

// ============================================================================
// Loading Operations
// ============================================================================

void RgTet4Element::applyBodyForce(const Vector3d& force, Vector& F) const
{
    int ndofs = kNodeCount * 3;
    if (F.size() != ndofs) F.resize(ndofs);
    
    Vector3d natCoord(m_gaussR[0], m_gaussS[0], m_gaussT[0]);
    std::vector<double> N;
    evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
    
    double jdet = evaluateJacobianDeterminant(natCoord);
    double weight = m_gaussW[0] * jdet;
    
    for (int i = 0; i < kNodeCount; ++i) {
        Vector3d fNode = weight * N[i] * force;
        F(3*i + 0) += fNode.x;
        F(3*i + 1) += fNode.y;
        F(3*i + 2) += fNode.z;
    }
}

void RgTet4Element::applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const
{
    int ndofs = kNodeCount * 3;
    if (F.size() != ndofs) F.resize(ndofs);
    // Placeholder for face integration
}

void RgTet4Element::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    int ndofs = kNodeCount * 3;
    if (F.size() != ndofs) F.resize(ndofs);
    
    if (nodeId >= 0 && nodeId < kNodeCount) {
        F(3*nodeId + 0) += force.x;
        F(3*nodeId + 1) += force.y;
        F(3*nodeId + 2) += force.z;
    }
}

// ============================================================================
// Material Point Access
// ============================================================================

FEMaterialPoint* RgTet4Element::getMaterialPoint(int gaussPtId)
{
    return nullptr;  // Placeholder
}

const FEMaterialPoint* RgTet4Element::getMaterialPoint(int gaussPtId) const
{
    return nullptr;  // Placeholder
}

// ============================================================================
// Serialization
// ============================================================================

void RgTet4Element::Serialize(DumpStream& ar)
{
    RgSolid3dElement::Serialize(ar);
    
    if (!ar.IsShallow()) {
        ar & m_gaussR & m_gaussS & m_gaussT & m_gaussW;
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

bool RgTet4Element::isValidNaturalCoordinate(const Vector3d& naturalCoord) const
{
    double sum = naturalCoord.x + naturalCoord.y + naturalCoord.z;
    return (naturalCoord.x >= -1e-10 && naturalCoord.y >= -1e-10 && 
            naturalCoord.z >= -1e-10 && sum <= 1.0 + 1e-10);
}

double RgTet4Element::getCharacteristicLength() const
{
    std::array<int, 2> edgeNodes;
    double length = 0.0;
    
    for (int i = 0; i < kNumEdges; ++i) {
        getEdgeNodeIds(i, edgeNodes);
        FENode* n0 = getNode(edgeNodes[0]);
        FENode* n1 = getNode(edgeNodes[1]);
        
        if (n0 && n1) {
            Vector3d v = n1->getPosition() - n0->getPosition();
            length += v.Length();
        }
    }
    
    return length / kNumEdges;
}

double RgTet4Element::getVolume() const
{
    // For tetrahedral element, compute volume from coordinates
    std::vector<Vector3d> coords(kNodeCount);
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            coords[i] = node->getPosition();
        }
    }
    
    // Volume = |det(v1, v2, v3)| / 6 where v_i = p_i - p_0
    Vector3d v1 = coords[1] - coords[0];
    Vector3d v2 = coords[2] - coords[0];
    Vector3d v3 = coords[3] - coords[0];
    
    double det = v1.x * (v2.y * v3.z - v2.z * v3.y) -
                 v1.y * (v2.x * v3.z - v2.z * v3.x) +
                 v1.z * (v2.x * v3.y - v2.y * v3.x);
    
    return std::abs(det) / 6.0;
}

// ============================================================================
// Protected Methods
// ============================================================================

void RgTet4Element::computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const
{
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    int ndofs = kNodeCount * 3;
    B.resize(6, ndofs);
    B.zero();
    
    for (int i = 0; i < kNodeCount; ++i) {
        double dN_dx = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
        double dN_dy = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
        double dN_dz = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
        
        B(0, 3*i + 0) = dN_dx;
        B(1, 3*i + 1) = dN_dy;
        B(2, 3*i + 2) = dN_dz;
        B(3, 3*i + 0) = dN_dy;
        B(3, 3*i + 1) = dN_dx;
        B(4, 3*i + 1) = dN_dz;
        B(4, 3*i + 2) = dN_dy;
        B(5, 3*i + 0) = dN_dz;
        B(5, 3*i + 2) = dN_dx;
    }
}

} // namespace RgFem
