#pragma once
#include "elements/RgSolidElement.h"
#include <array>
#include <vector>

// 20-node serendipity/hexahedral element
class RgHex20Element : public RgSolidElement
{
public:
    static constexpr int kNumNodes = 20;

    RgHex20Element();
    RgHex20Element(int id, const std::array<int, kNumNodes> &nodeIds);
    virtual ~RgHex20Element();

    // Basic info
    int nodeCount() const override { return kNumNodes; }
    void setNodeIds(const std::array<int,kNumNodes> &ids) { m_nodeIds = ids; }
    const std::array<int,kNumNodes>& nodeIds() const { return m_nodeIds; }

    // Natural coords xi = (r,s,t) in [-1,1]
    // N will be resized/filled to kNumNodes
    void shapeFunctions(const std::array<double,3> &xi, std::array<double,kNumNodes> &N) const override;
    // dN/dxi: for each node a 3-vector (dN/dr, dN/ds, dN/dt)
    void shapeFunctionDerivatives(const std::array<double,3> &xi,
                                  std::array<std::array<double,3>,kNumNodes> &dN) const override;

    // Return standard natural coordinates of the 20 nodes (order used by this element)
    static const std::array<std::array<double,3>, kNumNodes>& naturalNodeCoords();

protected:
    std::array<int,kNumNodes> m_nodeIds;

private:
    // Declaration only; define values in the .cpp file
    static const std::array<std::array<double,3>, kNumNodes> s_nodeCoords;
};