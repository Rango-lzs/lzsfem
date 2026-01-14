#include "elements/RgElement/RgElement.h"
#include "basicio/DumpStream.h"
#include "materials/RgMaterialPoint.h"
#include "femcore/Domain/RgDomain.h"
#include <math.h>
#include "materials/RgMaterial.h"

RgElement::RgElement() : m_pTraits(0) 
{ 
	static int n = 1;
	m_id = n++;
	m_part = nullptr;
	m_mat_id = -1;
	m_loc_id = -1;
}

// --- Element Identification ---
int RgElement::getId() const { return m_id; }

void RgElement::setId(int n) { m_id = n; }

int RgElement::getMatId() const { return m_mat_id; }

void RgElement::setMatId(int id) { m_mat_id = id; }

void RgElement::setLocalId(int lid)
{
    m_loc_id = lid;
}

int RgElement::getLocalId() const
{
    return m_loc_id;
}

// --- Mesh Partition ---
RgDomain* RgElement::getDomain() const
{
    return m_part;
}

void RgElement::setDomain(RgDomain* dom)
{
    m_part = dom;
}

RgMaterial* RgElement::getMaterial() const
{
    return getDomain()->GetMaterial();
}

// --- Node Connectivity ---
const std::vector<NodeId>& RgElement::getNodeIds() const
{
    return m_node;
}

NodeId RgElement::getNodeId(int idx) const
{
    if (idx >= 0 && idx < m_node.size())
    {
        return m_node[idx];
	}
    return -1;
}

NodeId RgElement::getLocNodeId(int idx) const
{
    if (idx >= 0 && idx < m_loc_node.size())
    {
        return m_loc_node[idx];
    }
    return -1;
}

void RgElement::setNodeId(int idx, int id)
{
    if (idx >= 0 && idx < m_node.size()) {
        m_node[idx] = id;
    }
}

void RgElement::setNode(FENode* n, int i)
{
    // Default implementation does nothing
}

FENode* RgElement::getNode(int idx) const
{
    return nullptr;
}

// --- Element Properties ---
ElementType RgElement::elementType() const
{
    return ElementType::FE_ELEM_INVALID_TYPE;
}

ElementCategory RgElement::elementCategory() const
{
    return ElementCategory::FE_ELEM_SOLID;
}

ElementShape RgElement::elementShape() const
{
    return m_pTraits ? m_pTraits->Shape() : ElementShape::ET_INVALID;
}

// --- Traits Management ---
void RgElement::initTraits()
{
    // Default implementation does nothing
}

RgElementTraits* RgElement::getTraits()
{
    return m_pTraits;
}

int RgElement::NodeSize() const
{
    return m_node.size();
}

int RgElement::GaussPointSize() const
{
    return m_pTraits ? m_pTraits->guassSize() : 0;
}

int RgElement::ShapeFunctions() const
{
    return m_pTraits ? m_pTraits->shapeSize() : 0;
}

// --- Material Point Data ---
RgMaterialPoint* RgElement::getMaterialPoint(int n)
{
    return (n >= 0 && n < m_state.size()) ? m_state[n] : nullptr;
}

// --- Material Point Data ---
 const RgMaterialPoint* RgElement::getMaterialPoint(int n) const
{
    //return (n >= 0 && n < m_state.size()) ? m_state[n] : nullptr;
     return nullptr;
}

void RgElement::setRgMaterialPointData(RgMaterialPoint* pmp, int n)
{
    if (pmp && n >= 0 && n < m_state.size()) {
        /* pmp->m_elem = this;
         pmp->m_index = n;*/
        m_state[n] = pmp;
    }
}

// --- Serialization ---
void RgElement::Serialize(DumpStream& ar)
{
    // Default implementation does nothing
}

// --- Field Evaluation ---
Vector3d RgElement::Evaluate(Vector3d* value, int iGauss)
{
    // Default implementation returns zero vector
    return Vector3d{0, 0, 0};
}

// --- Other Utilities ---
bool RgElement::isActive() const
{
    return true;
}

void RgElement::ClearData()
{
    // Default implementation does nothing
}
