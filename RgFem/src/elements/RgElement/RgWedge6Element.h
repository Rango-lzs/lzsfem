#pragma once

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>

namespace RgFem {

//! RgWedge6Element - 6-node linear wedge (triangular prism) element
//! 3D solid element with triangular cross-section
//! Node numbering:
//!   Bottom triangle: 0-1-2
//!   Top triangle:    3-4-5
class FEM_EXPORT RgWedge6Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 6;
    static constexpr int kNumFaces = 5;
    static constexpr int kNumEdges = 9;

    // Constructors and Destructors
    RgWedge6Element();
    explicit RgWedge6Element(const std::array<int, kNodeCount>& nodeIds);
    RgWedge6Element(const RgWedge6Element& other);
    RgWedge6Element& operator=(const RgWedge6Element& other);
    virtual ~RgWedge6Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementClass() const;

    // Element Properties
    virtual int getNumberOfNodes() const override { return kNodeCount; }
    virtual int getNumberOfGaussPoints() const;
    virtual int getNumberOfFaces() const { return kNumFaces; }
    virtual int getNumberOfEdges() const { return kNumEdges; }

    // Initialize element traits
    virtual void initTraits() override;

    // Shape functions and derivatives
    virtual double shapeFunction(int nodeId, double r, double s, double t) const override;
    virtual void shapeDerivatives(int nodeId, double r, double s, double t,
                                  double& dNdr, double& dNds, double& dNdt) const override;

    // Evaluate coordinates and Jacobian
    virtual void evaluateCoordinates(double r, double s, double t,
                                     std::array<double, 3>& coord) const override;
    virtual void evaluateJacobian(double r, double s, double t,
                                  std::array<std::array<double, 3>, 3>& J) const override;
    virtual double evaluateJacobianDeterminant(double r, double s, double t) const override;
    virtual void evaluateJacobianInverse(double r, double s, double t,
                                        std::array<std::array<double, 3>, 3>& Jinv) const override;

    // Matrix calculations
    virtual void calculateStiffnessMatrix(RgMatrix& K) const override;
    virtual void calculateMassMatrix(RgMatrix& M) const override;
    virtual void calculateInternalForceVector(RgVector& F) const override;

private:
    // Shape function helpers (triangular x linear in z)
    double N_tri(int triNodeId, double r, double s) const;
    void dN_tri_dr(int triNodeId, double r, double s, double& dNdr, double& dNds) const;
};

} // namespace RgFem
