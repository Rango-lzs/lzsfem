#include "materials/FEMaterialPoint.h"
#include "basicio/DumpStream.h"
#include <string.h>

FEMaterialPointData::FEMaterialPointData(FEMaterialPointData* ppt)
{
	m_pPrev = 0;
	m_pNext = ppt;
	if (ppt) ppt->m_pPrev = this;
}

FEMaterialPointData::~FEMaterialPointData()
{ 
	if (m_pNext) delete m_pNext;
	m_pNext = m_pPrev = 0;
}

void FEMaterialPointData::SetPrev(FEMaterialPointData* pt)
{
	m_pPrev = pt;
}

void FEMaterialPointData::SetNext(FEMaterialPointData* pt)
{
	m_pNext = pt;
	if (pt) pt->m_pPrev = this;
}

void FEMaterialPointData::Append(FEMaterialPointData* pt)
{
	if (pt == nullptr) return;
	if (m_pNext) m_pNext->Append(pt);
	else SetNext(pt);
}

void FEMaterialPointData::Init()
{
	if (m_pNext) m_pNext->Init();
}

void FEMaterialPointData::Update(const FETimeInfo& timeInfo)
{
	if (m_pNext) m_pNext->Update(timeInfo);
}

void FEMaterialPointData::Serialize(DumpStream& ar)
{
	if (m_pNext) m_pNext->Serialize(ar);
}

//=================================================================================================
FEMaterialPoint::FEMaterialPoint(FEMaterialPointData* data)
{
	m_data = data;
	m_elem = nullptr;
}

FEMaterialPoint::~FEMaterialPoint()
{

}

void FEMaterialPoint::Init()
{
	if (m_data) m_data->Init();
}

FEMaterialPoint* FEMaterialPoint::Copy()
{
	FEMaterialPoint* mp = new FEMaterialPoint(*this);
	if (m_data) mp->m_data = m_data->Copy();
	return mp;
}

void FEMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	if (m_data) m_data->Update(timeInfo);
}

void FEMaterialPoint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false)
	{
		ar & m_r0 & m_J0 & m_Jt;
	}
	if (m_data) m_data->Serialize(ar);
}

void FEMaterialPoint::Append(FEMaterialPointData* pt)
{
	if (pt == nullptr) return;
	assert(m_data);
	if (m_data) m_data->Append(pt);
}

//=================================================================================================

//-----------------------------------------------------------------------------
FEMaterialPointArray::FEMaterialPointArray(FEMaterialPointData* ppt) : FEMaterialPointData(ppt)
{
	
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::AddMaterialPoint(FEMaterialPoint* pt)
{
	m_mp.push_back(pt);
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Init()
{
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Init();
	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPointData::Update(timeInfo);
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Update(timeInfo);
}
