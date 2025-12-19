/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "stdafx.h"
#include "RgElementShapeStore.h"
#include "FEElement.h"
#include "FESolidElementShape.h"
#include "FESurfaceElementShape.h"

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
        m_pThis->RegisterShape(ET_TET5, new FETet5);
        m_pThis->RegisterShape(ET_TET10, new FETet10);
        m_pThis->RegisterShape(ET_TET15, new FETet15);
        m_pThis->RegisterShape(ET_TET20, new FETet20);
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

void RgElementShapeStore::RegisterShape(FE_Element_Shape shapeType, FEElementShape* pshape)
{
    m_ShapeMap[shapeType] = pshape;
}

//! return element shape class
FEElementShape* RgElementShapeStore::GetElementShapeClass(FE_Element_Shape eshape)
{
    auto it = m_ShapeMap.find(eshape);
    if (it != m_ShapeMap.end()) {
        return it->second;
    }
    return nullptr;
}

FE_Element_Shape RgElementShapeStore::GetElementShape(int ntype)
{
    // This function doesn't make sense in the new design since we're separating shapes from types
    // We'll provide a stub implementation for compatibility
    return FE_ELEM_INVALID_SHAPE;
}