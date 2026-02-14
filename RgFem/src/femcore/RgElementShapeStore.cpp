#include "RgElementShapeStore.h"
#include "elements/RgElement/RgElement.h"
#include "elements/ElementShape/RgSolidElementShape.h"
#include "elements/ElementShape/RgTet4Shape.h"
#include "elements/ElementShape/RgTet10Shape.h"
#include "elements/ElementShape/RgPenta6Shape.h"
#include "elements/ElementShape/RgPenta15Shape.h"
#include "elements/ElementShape/RgHex8Shape.h"
#include "elements/ElementShape/RgHex20Shape.h"
#include "elements/ElementShape/RgHex27Shape.h"
#include "elements/ElementShape/RgPyra5Shape.h"
#include "elements/ElementShape/RgPyra13Shape.h"
//#include "elements/ElementShape/RgSurfaceElementShape.h"

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
        m_pThis->RegisterShape(ET_TET4, new RgTet4Shape);
        m_pThis->RegisterShape(ET_TET10, new RgTet10Shape);
        m_pThis->RegisterShape(ET_PENTA6, new RgPenta6Shape);
        m_pThis->RegisterShape(ET_PENTA15, new RgPenta15Shape);
        m_pThis->RegisterShape(ET_HEX8, new RgHex8Shape);
        m_pThis->RegisterShape(ET_HEX20, new RgHex20Shape);
        m_pThis->RegisterShape(ET_HEX27, new RgHex27Shape);
        m_pThis->RegisterShape(ET_PYRA5, new RgPyra5Shape);
        m_pThis->RegisterShape(ET_PYRA13, new RgPyra13Shape);
        /* m_pThis->RegisterShape(ET_QUAD4, new FEQuad4);
         m_pThis->RegisterShape(ET_QUAD8, new FEQuad8);
         m_pThis->RegisterShape(ET_QUAD9, new FEQuad9);
         m_pThis->RegisterShape(ET_TRI3, new FETri3);
         m_pThis->RegisterShape(ET_TRI6, new FETri6);
         m_pThis->RegisterShape(ET_TRI7, new FETri7);
         m_pThis->RegisterShape(ET_TRI10, new FETri10);*/
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
