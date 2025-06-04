#pragma once
#include "basicio/DataRecord.h"


//-----------------------------------------------------------------------------
//! Base class for object log data (e.g. rigid bodies)
class FEM_EXPORT FELogObjectData : public FELogData
{
    DECLARE_META_CLASS(FELogObjectData, FELogData);

public:
	FELogObjectData() : FELogData() {}
	virtual ~FELogObjectData(){}
	//virtual double value(FERigidBody& rb) = 0;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT ObjectDataRecord : public DataRecord
{
public:
	ObjectDataRecord(FEModel* pfem);
	double Evaluate(int item, int ndata) override;
	void SetData(const char* sz) override;
	void SelectAllItems() override;
	int Size() const override;

private:
	std::vector<FELogObjectData*>	m_Data;
};
