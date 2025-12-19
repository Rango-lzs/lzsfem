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
#include "RgElementTraitsStore.h"
#include "FEElement.h"
#include "FESolidElementShape.h"
#include "FESurfaceElementShape.h"

RgElementTraitsStore* RgElementTraitsStore::m_pThis = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//! initialize library
void RgElementTraitsStore::Initialize()
{
    // Calling GetInstance will initialize the static pointer
    if (m_pThis == 0) GetInstance();
}

RgElementTraitsStore* RgElementTraitsStore::GetInstance()
{
    if (m_pThis == 0)
    {
        m_pThis = new RgElementTraitsStore;

        // register element types
        m_pThis->RegisterTraits(FE_HEX8G8, new FEHex8G8);
        m_pThis->RegisterTraits(FE_HEX8RI, new FEHex8RI);
        m_pThis->RegisterTraits(FE_HEX8G1, new FEHex8G1);
        m_pThis->RegisterTraits(FE_TET4G1, new FETet4G1);
        m_pThis->RegisterTraits(FE_TET4G4, new FETet4G4);
        m_pThis->RegisterTraits(FE_TET5G4, new FETet5G4);
        m_pThis->RegisterTraits(FE_PENTA6G6, new FEPenta6G6);
        m_pThis->RegisterTraits(FE_TET10G1, new FETet10G1);
        m_pThis->RegisterTraits(FE_TET10G4, new FETet10G4);
        m_pThis->RegisterTraits(FE_TET10G8, new FETet10G8);
        m_pThis->RegisterTraits(FE_TET10GL11, new FETet10GL11);
        m_pThis->RegisterTraits(FE_TET10G4RI1, new FETet10G4RI1);
        m_pThis->RegisterTraits(FE_TET10G8RI4, new FETet10G8RI4);
        m_pThis->RegisterTraits(FE_TET15G4, new FETet15G4);
        m_pThis->RegisterTraits(FE_TET15G8, new FETet15G8);
        m_pThis->RegisterTraits(FE_TET15G11, new FETet15G11);
        m_pThis->RegisterTraits(FE_TET15G15, new FETet15G15);
        m_pThis->RegisterTraits(FE_TET15G15RI4, new FETet15G15RI4);
        m_pThis->RegisterTraits(FE_TET20G15, new FETet20G15);
        m_pThis->RegisterTraits(FE_HEX20G8, new FEHex20G8);
        m_pThis->RegisterTraits(FE_HEX20G27, new FEHex20G27);
        m_pThis->RegisterTraits(FE_HEX27G27, new FEHex27G27);
        m_pThis->RegisterTraits(FE_PENTA15G8, new FEPenta15G8);
        m_pThis->RegisterTraits(FE_PENTA15G21, new FEPenta15G21);
        m_pThis->RegisterTraits(FE_PYRA5G8, new FEPyra5G8);
        m_pThis->RegisterTraits(FE_PYRA13G8, new FEPyra13G8);
        m_pThis->RegisterTraits(FE_QUAD4G4, new FEQuad4G4);
        m_pThis->RegisterTraits(FE_QUAD4NI, new FEQuad4NI);
        m_pThis->RegisterTraits(FE_TRI3G1, new FETri3G1);
        m_pThis->RegisterTraits(FE_TRI3G3, new FETri3G3);
        m_pThis->RegisterTraits(FE_TRI3G7, new FETri3G7);
        m_pThis->RegisterTraits(FE_TRI3NI, new FETri3NI);
        m_pThis->RegisterTraits(FE_TRI6G3, new FETri6G3);
        m_pThis->RegisterTraits(FE_TRI6G4, new FETri6G4);
        m_pThis->RegisterTraits(FE_TRI6G7, new FETri6G7);
        m_pThis->RegisterTraits(FE_TRI6GL7, new FETri6GL7);
        m_pThis->RegisterTraits(FE_TRI6NI, new FETri6NI);
        m_pThis->RegisterTraits(FE_TRI7G3, new FETri7G3);
        m_pThis->RegisterTraits(FE_TRI7G4, new FETri7G4);
        m_pThis->RegisterTraits(FE_TRI7G7, new FETri7G7);
        m_pThis->RegisterTraits(FE_TRI7GL7, new FETri7GL7);
        m_pThis->RegisterTraits(FE_TRI10G7, new FETri10G7);
        m_pThis->RegisterTraits(FE_TRI10G12, new FETri10G12);
        m_pThis->RegisterTraits(FE_QUAD8G9, new FEQuad8G9);
        m_pThis->RegisterTraits(FE_QUAD8NI, new FEQuad8NI);
        m_pThis->RegisterTraits(FE_QUAD9G9, new FEQuad9G9);
        m_pThis->RegisterTraits(FE_QUAD9NI, new FEQuad9NI);
        m_pThis->RegisterTraits(FE_SHELL_QUAD4G8, new FEShellQuad4G8);
        m_pThis->RegisterTraits(FE_SHELL_QUAD4G12, new FEShellQuad4G12);
        m_pThis->RegisterTraits(FE_SHELL_QUAD8G18, new FEShellQuad8G18);
        m_pThis->RegisterTraits(FE_SHELL_QUAD8G27, new FEShellQuad8G27);
        m_pThis->RegisterTraits(FE_SHELL_TRI3G6, new FEShellTri3G6);
        m_pThis->RegisterTraits(FE_SHELL_TRI3G9, new FEShellTri3G9);
        m_pThis->RegisterTraits(FE_SHELL_TRI6G14, new FEShellTri6G14);
        m_pThis->RegisterTraits(FE_SHELL_TRI6G21, new FEShellTri6G21);
        m_pThis->RegisterTraits(FE_TRUSS, new FETrussElementTraits);
        m_pThis->RegisterTraits(FE_DISCRETE, new FEDiscreteElementTraits);
        m_pThis->RegisterTraits(FE2D_TRI3G1, new FE2DTri3G1);
        m_pThis->RegisterTraits(FE2D_TRI6G3, new FE2DTri6G3);
        m_pThis->RegisterTraits(FE2D_QUAD4G4, new FE2DQuad4G4);
        m_pThis->RegisterTraits(FE2D_QUAD8G9, new FE2DQuad8G9);
        m_pThis->RegisterTraits(FE2D_QUAD9G9, new FE2DQuad9G9);
        m_pThis->RegisterTraits(FE_LINE2G1, new FELine2G1);
    }
    return m_pThis;
}

RgElementTraitsStore::~RgElementTraitsStore()
{
    for (auto& pair : m_TraitsMap) {
        delete pair.second;
    }
    m_TraitsMap.clear();
}

void RgElementTraitsStore::RegisterTraits(FE_Element_Type traitType, FEElementTraits* ptrait)
{
    m_TraitsMap[traitType] = ptrait;
}

void RgElementTraitsStore::SetElementTraits(FEElement& el, int nid)
{
    // In the new design, we need to find the trait by index
    // This is a bit awkward with map-only storage, so we'll iterate
    if (nid >= 0 && nid < (int)m_TraitsMap.size()) {
        auto it = m_TraitsMap.begin();
        std::advance(it, nid);
        el.SetTraits(it->second);
    }
}

//! return element traits data
FEElementTraits* RgElementTraitsStore::GetElementTraits(int ntype)
{
    auto it = m_TraitsMap.find((FE_Element_Type)ntype);
    if (it != m_TraitsMap.end()) {
        return it->second;
    }
    return nullptr;
}

//! return the element class of a given element type
FE_Element_Class RgElementTraitsStore::GetElementClass(int ntype)
{
    auto it = m_TraitsMap.find((FE_Element_Type)ntype);
    if (it != m_TraitsMap.end()) {
        return it->second->Class();
    }
    return FE_ELEM_INVALID_CLASS;
}

bool RgElementTraitsStore::IsValid(const FE_Element_Spec& c)
{
    if (c.eclass == FE_ELEM_INVALID_CLASS) return false;
    if (c.eshape == FE_ELEM_INVALID_SHAPE) return false;
    if (c.etype == FE_ELEM_INVALID_TYPE) return false;
    
    auto it = m_TraitsMap.find(c.etype);
    if (it == m_TraitsMap.end()) return false;
    
    if (c.eclass != it->second->Class()) return false;
    if (c.eshape != it->second->Shape()) return false;
    return true;
}

//! get the element spec from the type
FE_Element_Spec RgElementTraitsStore::GetElementSpecFromType(FE_Element_Type elemType)
{
    FE_Element_Spec espec;
    espec.etype = elemType;
    if (elemType != FE_ELEM_INVALID_TYPE)
    {
        auto it = m_TraitsMap.find(elemType);
        if (it != m_TraitsMap.end()) {
            espec.eclass = it->second->Class();
            espec.eshape = it->second->Shape();
        }
    }
    return espec;
}