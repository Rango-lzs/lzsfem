
#include "FEModelDataRecord.h"

FEModelLogData::FEModelLogData(FEModel* fem) : FELogData(fem) {}
FEModelLogData::~FEModelLogData() {}

FEModelDataRecord::FEModelDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_MODEL) {}
double FEModelDataRecord::Evaluate(int item, int ndata)
{
	assert(item == 0);
	return m_data[ndata]->value();
}

void FEModelDataRecord::ClearData()
{
	for (int i = 0; i < m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

void FEModelDataRecord::SetData(const char* szexpr)
{
	char szcopy[MAX_STRING] = { 0 };
	strcpy(szcopy, szexpr);
	char* sz = szcopy, * ch;
	ClearData();
	strcpy(m_szdata, szexpr);
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
        FEModelLogData* pdata = RANGO_NEW<FEModelLogData>(GetFEModel() , sz);
		if (pdata) m_data.push_back(pdata);
		else throw UnknownDataField(sz);
		sz = ch;
	} while (ch);
}

void FEModelDataRecord::SelectAllItems()
{
	std::vector<int> items;
	items.push_back(0);
	SetItemList(items);
}

int FEModelDataRecord::Size() const
{
	return 1;
}
