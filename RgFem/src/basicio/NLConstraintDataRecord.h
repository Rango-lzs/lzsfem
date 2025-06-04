#pragma once
#include "basicio/DataRecord.h"
#include "femcore/FENLConstraint.h"

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints log data (e.g. rigid connectors)
class FEM_EXPORT FELogNLConstraintData : public FELogData
{
    DECLARE_META_CLASS(FELogNLConstraintData, FELogData);

public:
    FELogNLConstraintData() : FELogData() {}
    virtual ~FELogNLConstraintData(){}
    virtual double value(FENLConstraint& rc) = 0;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT NLConstraintDataRecord : public DataRecord
{
public:
	NLConstraintDataRecord();
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
	int Size() const;
    
private:
    std::vector<FELogNLConstraintData*>	m_Data;
};
