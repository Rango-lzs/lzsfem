#pragma once

#include "RgLinearSolid2dElement.h"
#include <array>
#include <vector>

namespace RgFem {

//! RgTri3Element - 3-node linear triangular 2D element
//! 2D solid element with linear shape functions
//! Node numbering (counterclockwise):
//!   0: (0, 0)
//!   1: (1, 0)
//!   2: (0, 1)
class FEM_EXPORT RgTri3Element : public RgLinearSolid2dElement
{
public:
    static constexpr int kNodeCount = 3;
    static constexpr int kNumEdges = 3;

    // Constructors and Destructors
    RgTri3Element();
    explicit RgTri3Element(const std::array<int, kNodeCount>& nodeIds);
    RgTri3Element(const RgTri3Element& other);
    RgTri3Element& operator=(const RgTri3Element& other);
    virtual ~RgTri3Element();

    // Element Type Identification
    virtual ElementType elementType() const override;
    virtual ElementShape elementShape() const;
    virtual ElementCategory elementClass() const;

    // Element Properties
    virtual int getNumberOfNodes() const override { return kNodeCount; }
    virtual int getNumberOfGaussPoints() const;
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
    // Area calculation for Jacobian
    double getElementArea() const;
};

} // namespace RgFem
