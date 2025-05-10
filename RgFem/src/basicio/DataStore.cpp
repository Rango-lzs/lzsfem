#include "DataStore.h"
#include "logger/log.h"
#include "femcore/FEModel.h"
#include "femcore/FEAnalysis/FEAnalysis.h"

//-----------------------------------------------------------------------------
DataStore::DataStore()
{
}

//-----------------------------------------------------------------------------
DataStore::~DataStore()
{
}

//-----------------------------------------------------------------------------
void DataStore::Clear()
{
	for (size_t i=0; i<m_data.size(); ++i) delete m_data[i];
	m_data.clear();
}

//-----------------------------------------------------------------------------

void DataStore::Write()
{
	for (size_t i=0; i<m_data.size(); ++i)
	{
		DataRecord& DR = *m_data[i];
		DR.Write();
	}
}

//-----------------------------------------------------------------------------

void DataStore::AddRecord(DataRecord* prec)
{
	prec->m_nid = (int) m_data.size() + 1;
	m_data.push_back(prec);
}
