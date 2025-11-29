#pragma once

#include "RgLinearSolid3dElement.h"
#include <array>
#include <vector>

namespace RgFem {

//! RgPyramid5Element - 5-node linear pyramid element
//! 3D solid element with quadrilateral base and apex
//! Node numbering:
//!   Base quad: 0-1-2-3 (counterclockwise when viewed from apex)
//!   Apex:      4
class FEM_EXPORT RgPyramid5Element : public RgLinearSolid3dElement
{
public:
    static constexpr int kNodeCount = 5;
    static constexpr int kNumFaces = 5;
    static constexpr int kNumEdges = 8;

    // Constructors and Destructors
    RgPyramid5Element();
    explicit RgPyramid5Element(const std::array<int, kNodeCount>& nodeIds);
    RgPyramid5Element(const RgPyramid5Element& other);
    RgPyramid5Element& operator=(const RgPyramid5Element& other);
    virtual ~RgPyramid5Element();

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
    // Helper for shape function computation
    double computeHeightFactor(double t) const;
};

} // namespace RgFem
