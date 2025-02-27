#pragma once
#include "basicio/DataRecord.h"
#include "femcore/fem_export.h"

//-----------------------------------------------------------------------------
class FEM_EXPORT DataStore
{
public:
	DataStore();
	virtual ~DataStore();

	void Clear();

	void Write();

	void AddRecord(DataRecord* prec);

	int Size() { return (int) m_data.size(); }

	DataRecord* GetDataRecord(int i) { return m_data[i]; }

protected:
	std::vector<DataRecord*>	m_data;	//!< the data records
};
