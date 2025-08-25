#ifndef RGHEX8ELEMENT_H
#define RGHEX8ELEMENT_H

#include "RgSolid3dElement.h"
#include <array>
#include <vector>
#include <string>

namespace RgFem {

// forward declarations for project types (adjust to real names in your codebase)
class RgMatrix;
class RgVector;
class RgMaterial;

/*
    RgHex8Element
    - 8-node linear hexahedral solid element
    - derives from RgSolid3dElement
    - provides shape functions, jacobian, B-matrix helpers and standard element interfaces:
            * stiffness and mass computation (Gauss integration)
            * shape function evaluation and derivatives (dN/dxi)
            * conversion dN/dxi -> dN/dx
            * integration point access
    - All method signatures are kept generic (use your project's matrix/vector types).
*/
class RgHex8Element : public RgSolid3dElement {
public:
        static constexpr int kNodeCount = 8;

        RgHex8Element();
        explicit RgHex8Element(const std::array<int, kNodeCount>& nodeIds);
        RgHex8Element(const RgHex8Element& other);
        virtual ~RgHex8Element();

        // RgElement / RgSolid3dElement overrides
        virtual RgElement* clone() const override;
        virtual std::string typeName() const override;
        virtual int nodeCount() const override { return kNodeCount; }

        // N_i(xi,eta,zeta) for xi,eta,zeta in [-1,1]
        // N is filled with length kNodeCount
        virtual void shapeFunctions(double xi, double eta, double zeta, std::array<double, kNodeCount>& N) const;

        // dN/dxi (derivatives with respect to local coordinates xi,eta,zeta)
        // dN_dxi[i] = {dNi/dxi, dNi/deta, dNi/dzeta}
        virtual void shapeFunctionDerivatives(double xi, double eta, double zeta,
                                                                                    std::array<std::array<double,3>, kNodeCount>& dN_dxi) const;

        // Compute Jacobian J = [dx/dxi] (3x3), returns det(J)
        // coords: nodal coordinates as coords[node][{x,y,z}]
        virtual double jacobian(const std::array<std::array<double,3>, kNodeCount>& coords,
                                                        const std::array<std::array<double,3>, kNodeCount>& dN_dxi,
                                                        std::array<std::array<double,3>,3>& J,
                                                        std::array<std::array<double,3>,3>* invJ = nullptr) const;

        // Convert dN/dxi -> dN/dx using precomputed invJ (or compute inside if invJ==nullptr).
        // Fills dN_dx[node][{d/dx,d/dy,d/dz}] and returns detJ.
        virtual double computeDNdx(const std::array<std::array<double,3>, kNodeCount>& coords,
                                                             const std::array<std::array<double,3>, kNodeCount>& dN_dxi,
                                                             std::array<std::array<double,3>, kNodeCount>& dN_dx) const;

        // Build strain-displacement matrix B (3D small-strain linear elasticity)
        // B is expected to be a 6 x (3*kNodeCount) matrix (project type: RgMatrix)
        virtual void buildBMatrix(const std::array<std::array<double,3>, kNodeCount>& dN_dx, RgMatrix& B) const;

        // Element stiffness and consistent mass computation using Gauss integration
        // Ke and Me must be of appropriate size (3*kNodeCount square)
        virtual void computeStiffness(RgMatrix& Ke, const RgMaterial& material, int integrationOrder = 2) const override;
        virtual void computeMass(RgMatrix& Me, double density, int integrationOrder = 2) const override;

        // Given nodal displacement vector (length 3*kNodeCount), compute strains (6) and stresses (6) at one gauss point
        virtual void computeStrainStressAtGauss(const std::array<std::array<double,3>, kNodeCount>& coords,
                                                                                        const std::array<std::array<double,3>, kNodeCount>& dN_dxi,
                                                                                        const std::array<double, 3*kNodeCount>& nodalDisp,
                                                                                        const RgMaterial& material,
                                                                                        std::array<double,6>& strain,
                                                                                        std::array<double,6>& stress) const;

        // Integration point description
        struct GaussPoint {
                double xi, eta, zeta;
                double weight;
        };
        // Returns standard 2x2x2 Gauss points when order==2. Higher orders optional.
        static std::vector<GaussPoint> gaussPoints(int order = 2);

        // Accessors for node ids stored by element
        const std::array<int, kNodeCount>& nodeIds() const { return m_nodeIds; }
        void setNodeIds(const std::array<int, kNodeCount>& ids) { m_nodeIds = ids; }

protected:
        std::array<int, kNodeCount> m_nodeIds;

        // small helpers that may be reused in cpp implementation
        static double shape1D(double xi, int nodeLocalIndex); // returns shape factor for one direction
};

} // namespace RgFem

#endif // RGHEX8ELEMENT_H