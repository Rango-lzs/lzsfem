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
#include "FESurfaceMap.h"
#include "FESurface.h"
#include "DumpStream.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap() : FEDataMap(FE_SURFACE_MAP)
{
	m_maxFaceNodes = 0;
	m_format = FMT_MULT;
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(FEDataType dataType) : FEDataMap(FE_SURFACE_MAP, dataType)
{
	m_maxFaceNodes = 0;
	m_format = FMT_MULT;
}

//-----------------------------------------------------------------------------
FESurfaceMap::FESurfaceMap(const FESurfaceMap& map) : FEDataMap(map)
{
	m_maxFaceNodes = map.m_maxFaceNodes;
	m_format = map.m_format;
}

//-----------------------------------------------------------------------------
FESurfaceMap& FESurfaceMap::operator = (const FESurfaceMap& map)
{
	FEDataArray::operator=(map);
	m_name = map.m_name;
	m_maxFaceNodes = map.m_maxFaceNodes;
	return *this;
}

//-----------------------------------------------------------------------------
// return the item list associated with this map
FEItemList* FESurfaceMap::GetItemList()
{
	return const_cast<FEFacetSet*>(m_surf);
}

//-----------------------------------------------------------------------------
int FESurfaceMap::StorageFormat() const
{
	return m_format;
}

//-----------------------------------------------------------------------------
bool FESurfaceMap::Create(const FEFacetSet* ps, double val, Storage_Fmt fmt)
{
	m_surf = ps;
	m_format = fmt;
	if (fmt == FMT_MULT)
	{
		int NF = ps->Faces();
		m_maxFaceNodes = 0;
		for (int i = 0; i < NF; ++i)
		{
			const FEFacetSet::FACET& f = ps->Face(i);

			// TODO: currently, the number of nodes matches the type, but not sure if this will remain the case.
			if (f.ntype > m_maxFaceNodes) m_maxFaceNodes = f.ntype;
		}
		return resize(NF*m_maxFaceNodes, val);
	}
	else if (fmt == FMT_NODE)
	{
		FENodeList nodeList = ps->GetNodeList();
		int NN = nodeList.Size();
		m_maxFaceNodes = 1;
		return resize(NN, val);
	}
	else return false;
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, double v)
{
	int index = n*m_maxFaceNodes;
	for (int i=0; i<m_maxFaceNodes; ++i) set<double>(index+i, v);	
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const vec2d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<vec2d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const Vector3d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<Vector3d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const Matrix3d& v)
{
	int index = n*m_maxFaceNodes;
	for (int i = 0; i<m_maxFaceNodes; ++i) set<Matrix3d>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::setValue(int n, const Matrix3ds& v)
{
	int index = n * m_maxFaceNodes;
	for (int i = 0; i < m_maxFaceNodes; ++i) set<Matrix3ds>(index + i, v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(double v)
{
	set<double>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const vec2d& v)
{
	set<vec2d>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const Vector3d& v)
{
	set<Vector3d>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const Matrix3d& v)
{
	set<Matrix3d>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::fillValue(const Matrix3ds& v)
{
	set<Matrix3ds>(v);
}

//-----------------------------------------------------------------------------
void FESurfaceMap::Serialize(DumpStream& ar)
{
	FEDataMap::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_maxFaceNodes & m_format;
	if (ar.IsSaving())
	{
		FEFacetSet* fs = const_cast<FEFacetSet*>(m_surf);
		ar << fs;
	}
	else
	{
		FEFacetSet* fs = nullptr;
		ar >> fs;
		m_surf = fs;
	}
}

//-----------------------------------------------------------------------------
double FESurfaceMap::value(const FEMaterialPoint& pt)
{
	double v = 0.0;
	switch (m_format)
	{
	case FMT_NODE:
		{
			if (pt.m_elem)
			{
				// get the element this material point is in
				FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
				assert(pe);

				// make sure this element belongs to this domain
				// TODO: Can't check this if map was created through FEFacetSet
			//	assert(pe->GetMeshPartition() == m_dom);

				if (pt.m_index < 0x10000)
				{
					// integration point
					// get shape functions
					double* H = pe->H(pt.m_index);

					int ne = pe->Nodes();
					for (int i = 0; i < ne; ++i)
					{
						double vi = value<double>(pe->m_lnode[i], 0);
						v += vi * H[i];
					}
				}
				else
				{
					// element node
					int n = pt.m_index - 0x10000;
					v = value<double>(pe->m_lnode[n], 0);
				}
				return v;
			}
			else
			{
				// assume material point is a node
				return value<double>(pt.m_index, 0);
			}
		}
		break;
	case FMT_MULT:
		{
			// get the element this material point is in
			FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
			assert(pe);

			// make sure this element belongs to this domain
			// TODO: Can't check this if map was created through FEFacetSet
		//	assert(pe->GetMeshPartition() == m_dom);

			// get its local ID
			int lid = pe->GetLocalID();

			// get shape functions
			if (pt.m_index < 0x10000)
			{
				double* H = pe->H(pt.m_index);

				int ne = pe->Nodes();
				for (int i = 0; i < ne; ++i)
				{
					double vi = value<double>(lid, i);
					v += vi * H[i];
				}
			}
			else
			{
				// element node
				int n = pt.m_index - 0x10000;
				v = value<double>(lid, n);
			}
		}
		break;
	}

	return v;
}

//-----------------------------------------------------------------------------
Vector3d FESurfaceMap::valueVec3d(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
	assert(pe);

	// make sure this element belongs to this domain
	// TODO: Can't check this if map was created through FEFacetSet
	//	assert(pe->GetMeshPartition() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

	// get shape functions
	double* H = pe->H(pt.m_index);

	Vector3d v(0,0,0);
	int ne = pe->Nodes();
	for (int i = 0; i < ne; ++i)
	{
		Vector3d vi = value<Vector3d>(lid, i);
		v += vi*H[i];
	}

	return v;
}

//-----------------------------------------------------------------------------
Matrix3d FESurfaceMap::valueMat3d(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
	assert(pe);

	// make sure this element belongs to this domain
	// TODO: Can't check this if map was created through FEFacetSet
	//	assert(pe->GetMeshPartition() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

	// get shape functions
	double* H = pe->H(pt.m_index);

	Matrix3d v; v.zero();
	int ne = pe->Nodes();
	for (int i = 0; i < ne; ++i)
	{
		Matrix3d vi = value<Matrix3d>(lid, i);
		v += vi*H[i];
	}

	return v;
}


//-----------------------------------------------------------------------------
Matrix3ds FESurfaceMap::valueMat3ds(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FESurfaceElement* pe = dynamic_cast<FESurfaceElement*>(pt.m_elem);
	assert(pe);

	// make sure this element belongs to this domain
	// TODO: Can't check this if map was created through FEFacetSet
	//	assert(pe->GetMeshPartition() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

	// get shape functions
	double* H = pe->H(pt.m_index);

	Matrix3ds v; v.zero();
	int ne = pe->Nodes();
	for (int i = 0; i < ne; ++i)
	{
		Matrix3ds vi = value<Matrix3ds>(lid, i);
		v += vi * H[i];
	}

	return v;
}
