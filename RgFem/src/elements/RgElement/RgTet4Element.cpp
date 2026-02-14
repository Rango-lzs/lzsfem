#include "RgTet4Element.h"
#include "elements/ElementShape/RgTet4Shape.h"
#include "elements/ElementTraits/RgSolidElementTraits.h"
#include "materials/RgMaterialPoint.h"
#include "materials/RgMaterial.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "basicio/DumpStream.h"
#include "femcore/RgElementTraitsStore.h"
#include "elements/NaturalCoord.h"
#include <cmath>
#include <array>



// ============================================================================
// Constructor and Destructor
// ============================================================================

RgTet4Element::RgTet4Element()
    : RgLinearSolid3dElement()
{
    // Initialize traits for Tet4 element
    m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_TET4G1);
    m_node.resize(kNodeCount);
    m_loc_node.resize(kNodeCount);
}

RgTet4Element::RgTet4Element(const std::array<int, kNodeCount>& nodeIds)
    : RgLinearSolid3dElement()
{
    m_node.assign(nodeIds.begin(), nodeIds.end());
    m_pTraits = RgElementTraitsStore::GetInstance()->GetElementTraits(FE_TET4G1);
    m_loc_node.resize(kNodeCount);
}

RgTet4Element::RgTet4Element(const RgTet4Element& other)
    : RgLinearSolid3dElement(other)
{
  
}

RgTet4Element::~RgTet4Element()
{
}

RgTet4Element& RgTet4Element::operator=(const RgTet4Element& other)
{
    if (this != &other) {
        RgSolid3dElement::operator=(other);      
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

ElementCategory RgTet4Element::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}


void RgTet4Element::Serialize(DumpStream& ar)
{
}

void RgTet4Element::computeBMatrix(const NaturalCoord& naturalCoord, Matrix& B)
{
    // Compute strain-displacement matrix B
    // Relates nodal displacements to strain: epsilon = B * u
    
    std::vector<std::vector<double>> dN_dr_ds_dt = evalDeriv(NaturalCoord(naturalCoord.getR(), naturalCoord.getS(), naturalCoord.getT()));
    std::vector<double> dN_dr = dN_dr_ds_dt[0];
    std::vector<double> dN_ds = dN_dr_ds_dt[1];
    std::vector<double> dN_dt = dN_dr_ds_dt[2];
    
    Matrix3d JinvT = evaluateJacobianInverse(naturalCoord);
    JinvT = JinvT.transpose();
    
    // B matrix has shape [6, 12] for 3D solid element
    // 6 strain components (Voigt): {exx, eyy, ezz, exy, eyz, exz}
    // 12 DOFs: 3 per node * 4 nodes
    
    B.resize(6, 12);
    B.zero();
    
    for (int i = 0; i < kNodeCount; ++i) {
        // Compute physical derivatives: dN/dx, dN/dy, dN/dz
        // Using chain rule: dN/dx_i = dN/dξ_j * dξ_j/dx_i
        const double dN_dx = JinvT[0][0] * dN_dr[i] + JinvT[0][1] * dN_ds[i] + JinvT[0][2] * dN_dt[i];
        const double dN_dy = JinvT[1][0] * dN_dr[i] + JinvT[1][1] * dN_ds[i] + JinvT[1][2] * dN_dt[i];
        const double dN_dz = JinvT[2][0] * dN_dr[i] + JinvT[2][1] * dN_ds[i] + JinvT[2][2] * dN_dt[i];
        
        // Strain-displacement matrix for node i (3 DOFs: u, v, w)
        // Voigt notation: [εxx, εyy, εzz, γxy, γyz, γzx]
        B(0, 3*i)     = dN_dx;   // εxx = ∂u/∂x
        B(1, 3*i+1)   = dN_dy;   // εyy = ∂v/∂y
        B(2, 3*i+2)   = dN_dz;   // εzz = ∂w/∂z
        B(3, 3*i)     = dN_dy;   // γxy = ∂u/∂y
        B(3, 3*i+1)   = dN_dx;   // γxy = ∂v/∂x
        B(4, 3*i+1)   = dN_dz;   // γyz = ∂v/∂z
        B(4, 3*i+2)   = dN_dy;   // γyz = ∂w/∂y
        B(5, 3*i)     = dN_dz;   // γzx = ∂u/∂z
        B(5, 3*i+2)   = dN_dx;   // γzx = ∂w/∂x
    }
}


