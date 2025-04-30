#pragma once
#include "femcore/DataRecord.h"
#include "FENLConstraint.h"

//-----------------------------------------------------------------------------
//! Base class for nonlinear constraints log data (e.g. rigid connectors)
class FEM_EXPORT FELogNLConstraintData : public FELogData
{
    DECLARE_META_CLASS(FELogNLConstraintData, FELogData);

public:
    FELogNLConstraintData(FEModel* fem) : FELogData(fem) {}
    virtual ~FELogNLConstraintData(){}
    virtual double value(FENLConstraint& rc) = 0;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT NLConstraintDataRecord : public DataRecord
{
public:
	NLConstraintDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
	int Size() const;
    
private:
    std::vector<FELogNLConstraintData*>	m_Data;
};
