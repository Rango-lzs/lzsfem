#include "FaceDataRecord.h"
#include "femcore/FEModel.h"
#include "femcore/FESurface.h"
#include "femcore/FEFacetSet.h"

 DEFINE_META_CLASS(FELogFaceData, FELogData, "");

//-----------------------------------------------------------------------------
FELogFaceData::FELogFaceData(FEModel* fem) : FELogData(fem) {}

//-----------------------------------------------------------------------------
FELogFaceData::~FELogFaceData() {}

//-----------------------------------------------------------------------------
FaceDataRecord::FaceDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_FACE) 
{
	m_surface = nullptr;
}

//-----------------------------------------------------------------------------
int FaceDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
void FaceDataRecord::SetData(const char* szexpr)
{
	char szcopy[MAX_STRING] = { 0 };
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
        FELogFaceData* pdata = RANGO_NEW<FELogFaceData>(GetFEModel() ,sz);
		if (pdata) m_Data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
bool FaceDataRecord::Initialize()
{
	return (m_item.empty() == false);
}

//-----------------------------------------------------------------------------
double FaceDataRecord::Evaluate(int item, int ndata)
{
	int nface = item - 1;
	return m_Data[ndata]->value(m_surface->Element(nface));
}

//-----------------------------------------------------------------------------
void FaceDataRecord::SelectAllItems()
{
	assert(false);
}

void FaceDataRecord::SetItemList(FEItemList* itemList, const std::vector<int>& selection)
{
	FEFacetSet* facetSet = dynamic_cast<FEFacetSet*>(itemList); assert(facetSet);
	m_surface = facetSet->GetSurface(); assert(m_surface);
	int n = m_surface->Elements();
	if (selection.empty())
	{
		m_item.resize(n);
		for (int i = 0; i < n; ++i) m_item[i] = i + 1;
	}
	else
	{
		m_item = selection;
	}
}
