#include "RgElementShapeStore.h"
#include "elements/RgElement/RgElement.h"
#include "elements/ElementShape/RgSolidElementShape.h"
#include "elements/ElementShape/RgSurfaceElementShape.h"

RgElementShapeStore* RgElementShapeStore::m_pThis = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//! initialize library
void RgElementShapeStore::Initialize()
{
    // Calling GetInstance will initialize the static pointer
    if (m_pThis == 0) GetInstance();
}

RgElementShapeStore* RgElementShapeStore::GetInstance()
{
    if (m_pThis == 0)
    {
        m_pThis = new RgElementShapeStore;

        int n;
        // register element shapes (must be done before types!)
        m_pThis->RegisterShape(ET_TET4, new FETet4);
        m_pThis->RegisterShape(ET_TET10, new FETet10);
        m_pThis->RegisterShape(ET_PENTA6, new FEPenta6);
        m_pThis->RegisterShape(ET_PENTA15, new FEPenta15);
        m_pThis->RegisterShape(ET_HEX8, new FEHex8);
        m_pThis->RegisterShape(ET_HEX20, new FEHex20);
        m_pThis->RegisterShape(ET_HEX27, new FEHex27);
        m_pThis->RegisterShape(ET_PYRA5, new FEPyra5);
        m_pThis->RegisterShape(ET_PYRA13, new FEPyra13);
        m_pThis->RegisterShape(ET_QUAD4, new FEQuad4);
        m_pThis->RegisterShape(ET_QUAD8, new FEQuad8);
        m_pThis->RegisterShape(ET_QUAD9, new FEQuad9);
        m_pThis->RegisterShape(ET_TRI3, new FETri3);
        m_pThis->RegisterShape(ET_TRI6, new FETri6);
        m_pThis->RegisterShape(ET_TRI7, new FETri7);
        m_pThis->RegisterShape(ET_TRI10, new FETri10);
    }
    return m_pThis;
}

RgElementShapeStore::~RgElementShapeStore()
{
    for (auto& pair : m_ShapeMap) {
        delete pair.second;
    }
    m_ShapeMap.clear();
}

void RgElementShapeStore::RegisterShape(ElementShape shapeType, RgElementShape* pshape)
{
    m_ShapeMap[shapeType] = pshape;
}

//! return element shape class
RgElementShape* RgElementShapeStore::GetElementShape(ElementShape eshape)
{
    auto it = m_ShapeMap.find(eshape);
    if (it != m_ShapeMap.end()) {
        return it->second;
    }
    return nullptr;
}
