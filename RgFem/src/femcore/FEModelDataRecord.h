#pragma once
#include "DataRecord.h"

//-----------------------------------------------------------------------------
//! This is the base class for a model data
class FEM_EXPORT FEModelLogData : public FELogData
{
    META_CLASS_DECLARE(FEModelLogData, FELogData);

public:
    FEModelLogData(FEModel* fem);
    virtual ~FEModelLogData();
    virtual double value() = 0;
};

//-----------------------------------------------------------------------------
//! This class records model data
class FEM_EXPORT FEModelDataRecord : public DataRecord
{
public:
    FEModelDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz);
    void ClearData();
    void SelectAllItems();
    int Size() const;

private:
    std::vector<FEModelLogData*> m_data;
};
