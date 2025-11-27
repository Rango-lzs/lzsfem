#include "RgHex8Element.h"
#include "elements/ElementShape/RgHex8Shape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>

namespace RgFem {

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8Element::RgHex8Element()
    : RgSolid3dElement()
{
    // Node and Gauss point vectors will be resized when traits are initialized
}

RgHex8Element::RgHex8Element(const std::array<int, kNodeCount>& nodeIds)
    : RgSolid3dElement()
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
}

RgHex8Element::RgHex8Element(const RgHex8Element& other)
    : RgSolid3dElement(other)
{
    // Copy Gauss point data
    m_gaussR = other.m_gaussR;
    m_gaussS = other.m_gaussS;
    m_gaussT = other.m_gaussT;
    m_gaussW = other.m_gaussW;
    
    // Copy cached Jacobian data
    m_jacobianInverse = other.m_jacobianInverse;
    m_jacobianDet = other.m_jacobianDet;
}

RgHex8Element::~RgHex8Element()
{
    // Clean up if needed
}

RgHex8Element& RgHex8Element::operator=(const RgHex8Element& other)
{
    if (this != &other) {
        // Copy base class
        RgSolid3dElement::operator=(other);
        
        // Copy Gauss point data
        m_gaussR = other.m_gaussR;
        m_gaussS = other.m_gaussS;
        m_gaussT = other.m_gaussT;
        m_gaussW = other.m_gaussW;
        
        // Copy cached Jacobian data
        m_jacobianInverse = other.m_jacobianInverse;
        m_jacobianDet = other.m_jacobianDet;
    }
    return *this;
}

// ============================================================================
// Element Type Identification
// ============================================================================

ElementType RgHex8Element::elementType() const
{
    return ElementType::FE_HEX8G8;  // 8-node hex with 8 Gauss points (full integration)
}

ElementShape RgHex8Element::elementShape() const
{
    return ElementShape::ET_HEX8;
}

ElementCategory RgHex8Element::elementClass() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

// ============================================================================
// Element Properties
// ============================================================================

int RgHex8Element::getNumberOfGaussPoints() const
{
    // Full integration: 2x2x2 = 8 Gauss points
    // Reduced integration: 1x1x1 = 1 Gauss point
    return m_pTraits ? m_pTraits->m_nint : 8;
}

// ============================================================================
// Initialize Element Traits
// ============================================================================

void RgHex8Element::initTraits()
{
    // Create traits for HEX8 element with 8 Gauss points (full integration)
    RgSolidElementTraits* traits = new RgSolidElementTraits(8, kNodeCount, ET_HEX8, FE_HEX8G8);
    traits->init();
    
    // Set element traits (inherited from RgElement)
    // Note: This is typically called during element construction or by the library
    if (m_pTraits == nullptr) {
        m_pTraits = traits;
        m_node.resize(kNodeCount);
        m_loc_node.resize(kNodeCount);
    }
    
    // Initialize Gauss points
    initializeGaussPoints();
    
    // Resize Jacobian cache
    m_jacobianInverse.resize(getNumberOfGaussPoints());
    m_jacobianDet.resize(getNumberOfGaussPoints());
}

void RgHex8Element::initializeGaussPoints()
{
    int npts = getNumberOfGaussPoints();
    
    if (npts == 1) {
        // Reduced integration: 1 point at origin
        m_gaussR.assign(1, 0.0);
        m_gaussS.assign(1, 0.0);
        m_gaussT.assign(1, 0.0);
        m_gaussW.assign(1, 8.0);  // Weight = 2^3 = 8
    }
    else if (npts == 8) {
        // Full integration: 2x2x2 Gauss points
        m_gaussR.resize(8);
        m_gaussS.resize(8);
        m_gaussT.resize(8);
        m_gaussW.resize(8);
        
        double gp = 1.0 / std::sqrt(3.0);  // ±1/√3 ≈ ±0.57735
        double w = 1.0;  // Weight for each point
        
        int idx = 0;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 2; ++k) {
                    m_gaussR[idx] = (i == 0) ? -gp : gp;
                    m_gaussS[idx] = (j == 0) ? -gp : gp;
                    m_gaussT[idx] = (k == 0) ? -gp : gp;
                    m_gaussW[idx] = w;
                    ++idx;
                }
            }
        }
    }
}

// ============================================================================
// Shape Function Evaluations
// ============================================================================

void RgHex8Element::evaluateShapeFunctions(double r, double s, double t, std::vector<double>& N) const
{
    N.resize(kNodeCount);
    
    // Trilinear shape functions for 8-node hexahedral element
    // Node numbering:
    // 0: (-1, -1, -1), 1: (+1, -1, -1), 2: (+1, +1, -1), 3: (-1, +1, -1)
    // 4: (-1, -1, +1), 5: (+1, -1, +1), 6: (+1, +1, +1), 7: (-1, +1, +1)
    
    double r1 = 1.0 - r;
    double r2 = 1.0 + r;
    double s1 = 1.0 - s;
    double s2 = 1.0 + s;
    double t1 = 1.0 - t;
    double t2 = 1.0 + t;
    
    N[0] = 0.125 * r1 * s1 * t1;
    N[1] = 0.125 * r2 * s1 * t1;
    N[2] = 0.125 * r2 * s2 * t1;
    N[3] = 0.125 * r1 * s2 * t1;
    N[4] = 0.125 * r1 * s1 * t2;
    N[5] = 0.125 * r2 * s1 * t2;
    N[6] = 0.125 * r2 * s2 * t2;
    N[7] = 0.125 * r1 * s2 * t2;
}

void RgHex8Element::evaluateShapeDerivatives(double r, double s, double t,
                                             std::vector<double>& dN_dr,
                                             std::vector<double>& dN_ds,
                                             std::vector<double>& dN_dt) const
{
    dN_dr.resize(kNodeCount);
    dN_ds.resize(kNodeCount);
    dN_dt.resize(kNodeCount);
    
    double s1 = 1.0 - s;
    double s2 = 1.0 + s;
    double t1 = 1.0 - t;
    double t2 = 1.0 + t;
    
    double r1 = 1.0 - r;
    double r2 = 1.0 + r;
    double r1t = 1.0 - r;
    double r2t = 1.0 + r;
    
    double t1s = 1.0 - t;
    double t2s = 1.0 + t;
    
    // dN/dr
    dN_dr[0] = -0.125 * s1 * t1;
    dN_dr[1] =  0.125 * s1 * t1;
    dN_dr[2] =  0.125 * s2 * t1;
    dN_dr[3] = -0.125 * s2 * t1;
    dN_dr[4] = -0.125 * s1 * t2;
    dN_dr[5] =  0.125 * s1 * t2;
    dN_dr[6] =  0.125 * s2 * t2;
    dN_dr[7] = -0.125 * s2 * t2;
    
    // dN/ds
    dN_ds[0] = -0.125 * r1 * t1;
    dN_ds[1] = -0.125 * r2 * t1;
    dN_ds[2] =  0.125 * r2 * t1;
    dN_ds[3] =  0.125 * r1 * t1;
    dN_ds[4] = -0.125 * r1 * t2;
    dN_ds[5] = -0.125 * r2 * t2;
    dN_ds[6] =  0.125 * r2 * t2;
    dN_ds[7] =  0.125 * r1 * t2;
    
    // dN/dt
    dN_dt[0] = -0.125 * r1 * s1;
    dN_dt[1] = -0.125 * r2 * s1;
    dN_dt[2] = -0.125 * r2 * s2;
    dN_dt[3] = -0.125 * r1 * s2;
    dN_dt[4] =  0.125 * r1 * s1;
    dN_dt[5] =  0.125 * r2 * s1;
    dN_dt[6] =  0.125 * r2 * s2;
    dN_dt[7] =  0.125 * r1 * s2;
}

void RgHex8Element::evaluateShapeDerivatives2(double r, double s, double t,
                                              std::vector<double>& d2N_drr,
                                              std::vector<double>& d2N_dss,
                                              std::vector<double>& d2N_dtt,
                                              std::vector<double>& d2N_drs,
                                              std::vector<double>& d2N_dst,
                                              std::vector<double>& d2N_drt) const
{
    // For trilinear shape functions, second derivatives are zero except for mixed derivatives
    d2N_drr.assign(kNodeCount, 0.0);
    d2N_dss.assign(kNodeCount, 0.0);
    d2N_dtt.assign(kNodeCount, 0.0);
    
    d2N_drs.resize(kNodeCount);
    d2N_dst.resize(kNodeCount);
    d2N_drt.resize(kNodeCount);
    
    double t1 = 1.0 - t;
    double t2 = 1.0 + t;
    double s1 = 1.0 - s;
    double s2 = 1.0 + s;
    double r1 = 1.0 - r;
    double r2 = 1.0 + r;
    
    // d2N/drds
    d2N_drs[0] =  0.125 * t1;
    d2N_drs[1] = -0.125 * t1;
    d2N_drs[2] =  0.125 * t1;
    d2N_drs[3] = -0.125 * t1;
    d2N_drs[4] =  0.125 * t2;
    d2N_drs[5] = -0.125 * t2;
    d2N_drs[6] =  0.125 * t2;
    d2N_drs[7] = -0.125 * t2;
    
    // d2N/dsdt
    d2N_dst[0] =  0.125 * r1;
    d2N_dst[1] =  0.125 * r2;
    d2N_dst[2] = -0.125 * r2;
    d2N_dst[3] = -0.125 * r1;
    d2N_dst[4] = -0.125 * r1;
    d2N_dst[5] = -0.125 * r2;
    d2N_dst[6] =  0.125 * r2;
    d2N_dst[7] =  0.125 * r1;
    
    // d2N/drdt
    d2N_drt[0] =  0.125 * s1;
    d2N_drt[1] = -0.125 * s1;
    d2N_drt[2] = -0.125 * s2;
    d2N_drt[3] =  0.125 * s2;
    d2N_drt[4] = -0.125 * s1;
    d2N_drt[5] =  0.125 * s1;
    d2N_drt[6] =  0.125 * s2;
    d2N_drt[7] = -0.125 * s2;
}

// ============================================================================
// Geometric Operations
// ============================================================================

Vector3d RgHex8Element::evaluateCoordinates(const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d coords(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            Vector3d nodePos = node->getPosition();
            coords += N[i] * nodePos;
        }
    }
    return coords;
}

Matrix3d RgHex8Element::evaluateJacobian(const Vector3d& naturalCoord) const
{
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z, 
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d J(0, 0, 0, 0, 0, 0, 0, 0, 0);
    
    for (int i = 0; i < kNodeCount; ++i) {
        FENode* node = getNode(i);
        if (node) {
            Vector3d pos = node->getPosition();
            
            // J[0][0] = dx/dr, J[0][1] = dy/dr, J[0][2] = dz/dr
            J.m[0][0] += dN_dr[i] * pos.x;
            J.m[0][1] += dN_dr[i] * pos.y;
            J.m[0][2] += dN_dr[i] * pos.z;
            
            // J[1][0] = dx/ds, J[1][1] = dy/ds, J[1][2] = dz/ds
            J.m[1][0] += dN_ds[i] * pos.x;
            J.m[1][1] += dN_ds[i] * pos.y;
            J.m[1][2] += dN_ds[i] * pos.z;
            
            // J[2][0] = dx/dt, J[2][1] = dy/dt, J[2][2] = dz/dt
            J.m[2][0] += dN_dt[i] * pos.x;
            J.m[2][1] += dN_dt[i] * pos.y;
            J.m[2][2] += dN_dt[i] * pos.z;
        }
    }
    return J;
}

double RgHex8Element::evaluateJacobianDeterminant(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.det();
}

Matrix3d RgHex8Element::evaluateJacobianInverse(const Vector3d& naturalCoord) const
{
    Matrix3d J = evaluateJacobian(naturalCoord);
    return J.inverse();
}

// ============================================================================
// Physical Field Evaluations
// ============================================================================

Vector3d RgHex8Element::evaluateField(const Vector3d* nodeValues, const Vector3d& naturalCoord) const
{
    std::vector<double> N;
    evaluateShapeFunctions(naturalCoord.x, naturalCoord.y, naturalCoord.z, N);
    
    Vector3d field(0, 0, 0);
    for (int i = 0; i < kNodeCount; ++i) {
        field += N[i] * nodeValues[i];
    }
    return field;
}

double RgHex8Element::evaluateScalarField(const double* nodeValues, const Vector3d& naturalCoord) const
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

void RgHex8Element::calculateStiffnessMatrix(Matrix& K) const
{
    // Get element DOFs (24 DOFs for 8 nodes with 3 DOF per node)
    int ndofs = kNodeCount * 3;
    K.resize(ndofs, ndofs);
    K.zero();
    
    // Integrate stiffness matrix using Gauss quadrature
    int npts = getNumberOfGaussPoints();
    
    for (int gp = 0; gp < npts; ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[gp] * jdet;
        
        // Compute B matrix at this Gauss point
        Matrix B;
        computeBMatrix(natCoord, B);
        
        // K += B^T * D * B * weight
        // (This is simplified - actual implementation depends on material model)
        // For now, we provide the structure; material model would be integrated here
    }
}

void RgHex8Element::calculateMassMatrix(Matrix& M) const
{
    // Get element DOFs
    int ndofs = kNodeCount * 3;
    M.resize(ndofs, ndofs);
    M.zero();
    
    // Get material density from element material
    // (Actual implementation depends on FEMaterial interface)
    double rho = 1.0;  // Default density
    
    int npts = getNumberOfGaussPoints();
    
    for (int gp = 0; gp < npts; ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        std::vector<double> N;
        evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[gp] * jdet * rho;
        
        // Assemble consistent mass matrix: M = integral of rho * N^T * N dV
        for (int i = 0; i < kNodeCount; ++i) {
            for (int j = 0; j < kNodeCount; ++j) {
                double Mij = weight * N[i] * N[j];
                for (int d = 0; d < 3; ++d) {
                    M(3*i + d, 3*j + d) += Mij;
                }
            }
        }
    }
}

void RgHex8Element::calculateDampingMatrix(Matrix& C) const
{
    // Rayleigh damping: C = alpha*M + beta*K
    // For now, initialize as zero; actual implementation would use material damping parameters
    int ndofs = kNodeCount * 3;
    C.resize(ndofs, ndofs);
    C.zero();
}

void RgHex8Element::calculateTangentStiffnessMatrix(Matrix& Kt) const
{
    // Tangent stiffness for nonlinear analysis
    // Kt = K + Kg (material stiffness + geometric stiffness)
    int ndofs = kNodeCount * 3;
    Kt.resize(ndofs, ndofs);
    Kt.zero();
    
    // Calculate standard stiffness
    Matrix K;
    calculateStiffnessMatrix(K);
    Kt = K;
    
    // Add geometric stiffness for large deformations
    Matrix Kg;
    computeGeometricStiffness(Kg);
    Kt += Kg;
}

void RgHex8Element::calculateInternalForceVector(Vector& F) const
{
    // Internal force: F = integral of B^T * sigma dV
    // This would require access to current stress state
    // Placeholder implementation
}

// ============================================================================
// Strain and Stress Calculations
// ============================================================================

void RgHex8Element::calculateStress(FEMaterialPoint& matPt, Matrix3ds& stress)
{
    // Calculate stress from strain using material model
    // (Implementation depends on material interface)
    // stress = D * strain (elasticity) or more complex for nonlinear materials
}

void RgHex8Element::calculateStrain(FEMaterialPoint& matPt, Matrix3ds& strain)
{
    // Calculate strain from nodal displacements
    // epsilon = B * u, where B is strain-displacement matrix
    // (Implementation depends on displacement storage)
}

double RgHex8Element::calculateStrainEnergy() const
{
    // U = 0.5 * integral of sigma : epsilon dV
    return 0.0;  // Placeholder
}

double RgHex8Element::calculateKineticEnergy() const
{
    // KE = 0.5 * integral of rho * v^2 dV
    return 0.0;  // Placeholder
}

// ============================================================================
// Face and Edge Operations
// ============================================================================

void RgHex8Element::getFaceNodeIds(int faceId, std::array<int, 4>& faceNodes) const
{
    // Hexahedron faces (quad, 4 nodes each)
    // Face numbering follows standard FEA convention
    static const int faces[6][4] = {
        {0, 1, 5, 4},  // Face 0: bottom face (z=-1)
        {1, 2, 6, 5},  // Face 1: right face (r=+1)
        {2, 3, 7, 6},  // Face 2: top face (z=+1)
        {3, 0, 4, 7},  // Face 3: left face (r=-1)
        {0, 3, 2, 1},  // Face 4: front face (s=-1)
        {4, 5, 6, 7}   // Face 5: back face (s=+1)
    };
    
    if (faceId >= 0 && faceId < 6) {
        std::copy(faces[faceId], faces[faceId] + 4, faceNodes.begin());
    }
}

void RgHex8Element::getEdgeNodeIds(int edgeId, std::array<int, 2>& edgeNodes) const
{
    // Hexahedron edges (2 nodes each)
    static const int edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},  // Bottom face edges
        {4, 5}, {5, 6}, {6, 7}, {7, 4},  // Top face edges
        {0, 4}, {1, 5}, {2, 6}, {3, 7}   // Vertical edges
    };
    
    if (edgeId >= 0 && edgeId < 12) {
        std::copy(edges[edgeId], edges[edgeId] + 2, edgeNodes.begin());
    }
}

// ============================================================================
// Loading Operations
// ============================================================================

void RgHex8Element::applyBodyForce(const Vector3d& force, Vector& F) const
{
    // Distribute body force to nodes using shape functions
    int ndofs = kNodeCount * 3;
    if (F.size() != ndofs) F.resize(ndofs);
    
    int npts = getNumberOfGaussPoints();
    
    for (int gp = 0; gp < npts; ++gp) {
        Vector3d natCoord(m_gaussR[gp], m_gaussS[gp], m_gaussT[gp]);
        std::vector<double> N;
        evaluateShapeFunctions(natCoord.x, natCoord.y, natCoord.z, N);
        
        double jdet = evaluateJacobianDeterminant(natCoord);
        double weight = m_gaussW[gp] * jdet;
        
        for (int i = 0; i < kNodeCount; ++i) {
            Vector3d fNode = weight * N[i] * force;
            F(3*i + 0) += fNode.x;
            F(3*i + 1) += fNode.y;
            F(3*i + 2) += fNode.z;
        }
    }
}

void RgHex8Element::applyDistributedLoad(int faceId, const Vector3d& traction, Vector& F) const
{
    // Apply distributed load on element face using surface integration
    // (Implementation requires surface Gauss quadrature)
    int ndofs = kNodeCount * 3;
    if (F.size() != ndofs) F.resize(ndofs);
    
    // This would require 2D Gauss integration over the face
    // Placeholder for actual implementation
}

void RgHex8Element::applyPointLoad(int nodeId, const Vector3d& force, Vector& F) const
{
    // Apply concentrated load at node
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

FEMaterialPoint* RgHex8Element::getMaterialPoint(int gaussPtId)
{
    if (gaussPtId >= 0 && gaussPtId < getNumberOfGaussPoints()) {
        return m_pTraits ? nullptr : nullptr;  // Return material point at Gauss point
    }
    return nullptr;
}

const FEMaterialPoint* RgHex8Element::getMaterialPoint(int gaussPtId) const
{
    if (gaussPtId >= 0 && gaussPtId < getNumberOfGaussPoints()) {
        return nullptr;  // Return const material point at Gauss point
    }
    return nullptr;
}

// ============================================================================
// State Management
// ============================================================================

void RgHex8Element::updateState(double timeStep)
{
    // Update material state at all Gauss points
    // (Implementation depends on material model)
}

void RgHex8Element::commitState()
{
    // Commit state variables at all Gauss points
    // (Implementation depends on material model)
}

void RgHex8Element::resetState()
{
    // Reset state variables at all Gauss points
    // (Implementation depends on material model)
}

// ============================================================================
// Serialization
// ============================================================================

void RgHex8Element::Serialize(DumpStream& ar)
{
    // Serialize base class
    RgSolid3dElement::Serialize(ar);
    
    if (!ar.IsShallow()) {
        ar & m_gaussR & m_gaussS & m_gaussT & m_gaussW;
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

bool RgHex8Element::isValidNaturalCoordinate(const Vector3d& naturalCoord) const
{
    // Check if point is within reference hexahedron [-1, 1]^3
    const double tolerance = 1.0 + 1e-10;
    return (std::abs(naturalCoord.x) <= tolerance &&
            std::abs(naturalCoord.y) <= tolerance &&
            std::abs(naturalCoord.z) <= tolerance);
}

double RgHex8Element::getCharacteristicLength() const
{
    // Compute characteristic length (e.g., mean edge length)
    double length = 0.0;
    std::array<int, 2> edgeNodes;
    
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

// ============================================================================
// Protected Methods
// ============================================================================

void RgHex8Element::computeBMatrix(const Vector3d& naturalCoord, Matrix& B) const
{
    // Compute strain-displacement matrix B
    // Relates nodal displacements to strain: epsilon = B * u
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    evaluateShapeDerivatives(naturalCoord.x, naturalCoord.y, naturalCoord.z,
                            dN_dr, dN_ds, dN_dt);
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    // B matrix has shape [6, 24] for 3D solid element
    // 6 strain components (Voigt): {exx, eyy, ezz, exy, eyz, exz}
    // 24 DOFs: 3 per node * 8 nodes
    
    B.resize(6, 24);
    B.zero();
    
    for (int i = 0; i < kNodeCount; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        double dN_dx = dN_dr[i] * JinvT.m[0][0] + dN_ds[i] * JinvT.m[1][0] + dN_dt[i] * JinvT.m[2][0];
        double dN_dy = dN_dr[i] * JinvT.m[0][1] + dN_ds[i] * JinvT.m[1][1] + dN_dt[i] * JinvT.m[2][1];
        double dN_dz = dN_dr[i] * JinvT.m[0][2] + dN_ds[i] * JinvT.m[1][2] + dN_dt[i] * JinvT.m[2][2];
        
        // Fill B matrix
        B(0, 3*i + 0) = dN_dx;  // exx
        B(1, 3*i + 1) = dN_dy;  // eyy
        B(2, 3*i + 2) = dN_dz;  // ezz
        B(3, 3*i + 0) = dN_dy;  // exy (shear)
        B(3, 3*i + 1) = dN_dx;
        B(4, 3*i + 1) = dN_dz;  // eyz (shear)
        B(4, 3*i + 2) = dN_dy;
        B(5, 3*i + 0) = dN_dz;  // exz (shear)
        B(5, 3*i + 2) = dN_dx;
    }
}

void RgHex8Element::computeGeometricStiffness(Matrix& Kg) const
{
    // Compute geometric (initial stress) stiffness matrix
    // Kg = integral of B_geo^T * sigma * B_geo dV
    // This accounts for large deformation effects
    
    int ndofs = kNodeCount * 3;
    Kg.resize(ndofs, ndofs);
    Kg.zero();
    
    // (Actual implementation requires current stress state)
}

} // namespace RgFem
