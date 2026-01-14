#include "RgHex8Element.h"
#include "elements/ElementShape/RgHex8Shape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/RgMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include "femcore/RgElementTraitsStore.h"
#include "../RgElementState.h"
#include "materials/RgMaterial.h"
#include "materials/SmallDeformation/RgElastoPlasticMaterialPoint.h"
#include "elements/NaturalCoord.h"

//应力应变存储在哪儿，由谁计算

// ============================================================================
// Constructor and Destructor
// ============================================================================

RgHex8Element::RgHex8Element()
    : RgLinearSolid3dElement()
{
    // Default to full integration
    m_node.resize(kNodeCount);
    m_loc_node.resize(kNodeCount);
}

RgHex8Element::RgHex8Element(bool fullInt)
    : RgLinearSolid3dElement()
{
    if (fullInt)
    {
        //m_pTraits = RgHex8Element::fullIntTraits();
        m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_HEX8G8);
    }
    else
    {
        //m_pTraits = RgHex8Element::reduceIntTraits();
        m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_HEX8G1);
    }
    m_node.resize(kNodeCount);
    m_loc_node.resize(kNodeCount);
}


RgHex8Element::RgHex8Element(const RgHex8Element& other)
    : RgLinearSolid3dElement(other)
{
    // m_pTraits is shared via static methods, no need to copy
    // m_node and m_loc_node will be handled by base class
}

RgHex8Element::~RgHex8Element()
{
    // Nothing specific to clean up - traits are static and managed automatically
}

RgHex8Element& RgHex8Element::operator=(const RgHex8Element& other)
{
    if (this != &other) {
        // Copy base class - this will handle nodes and other data
        RgLinearSolid3dElement::operator=(other);
        // Note: m_pTraits is not copied as it's a shared static object
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

ElementCategory RgHex8Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

// ============================================================================
// Serialization
// ============================================================================

void RgHex8Element::Serialize(DumpStream& ar)
{
    // Serialize base class
    RgLinearSolid3dElement::Serialize(ar);
}

// ============================================================================

void RgHex8Element::computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B)
{
    // Compute strain-displacement matrix B
    // Relates nodal displacements to strain: epsilon = B * u
    
    std::vector<double> dN_dr, dN_ds, dN_dt;
    std::vector<std::vector<double>> H_Deriv = evalDeriv(NaturalCoord(naturalCoord.getR(), naturalCoord.getS(), naturalCoord.getT()));
    dN_dr = H_Deriv[0];
    dN_ds = H_Deriv[1];
    dN_ds = H_Deriv[2];
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    
    // B matrix has shape [6, 24] for 3D solid element
    // 6 strain components (Voigt): {exx, eyy, ezz, exy, eyz, exz}
    // 24 DOFs: 3 per node * 8 nodes
    
    B.resize(6, 24);
    B.zero(); // Use setZero() for better performance with Eigen
    
    for (int i = 0; i < kNodeCount; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        // Using chain rule: dN/dx_i = dN/dξ_j * dξ_j/dx_i
        // where ξ_j are natural coordinates and x_i are physical coordinates
        const double dN_dx = dN_dr[i] * JinvT[0][0] + dN_ds[i] * JinvT[1][0] + dN_dt[i] * JinvT[2][0];
        const double dN_dy = dN_dr[i] * JinvT[0][1] + dN_ds[i] * JinvT[1][1] + dN_dt[i] * JinvT[2][1];
        const double dN_dz = dN_dr[i] * JinvT[0][2] + dN_ds[i] * JinvT[1][2] + dN_dt[i] * JinvT[2][2];
        
        // Fill B matrix according to Voigt notation for strain components
        // ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz (engineering shear strains)
        B(0, 3*i + 0) = dN_dx;               // ε_xx = du/dx
        B(1, 3*i + 1) = dN_dy;               // ε_yy = dv/dy
        B(2, 3*i + 2) = dN_dz;               // ε_zz = dw/dz
        B(3, 3*i + 0) = dN_dy;               // γ_xy = du/dy + dv/dx
        B(3, 3*i + 1) = dN_dx;
        B(4, 3*i + 1) = dN_dz;               // γ_yz = dv/dz + dw/dy
        B(4, 3*i + 2) = dN_dy;
        B(5, 3*i + 0) = dN_dz;               // γ_xz = du/dz + dw/dx
        B(5, 3*i + 2) = dN_dx;
    }
}

