#pragma once
#include "basicio/DataRecord.h"

class RgElement;
class RgElementSet;

//-----------------------------------------------------------------------------
//! Base class for element log data
class FEM_EXPORT FELogElemData : public FELogData
{
    DECLARE_META_CLASS(FELogElemData, FELogData);

public:
    FELogElemData();
    virtual ~FELogElemData();
    virtual double value(RgElement& el) = 0;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT ElementDataRecord : public DataRecord
{
    DECLARE_META_CLASS(ElementDataRecord, DataRecord);
    struct ELEMREF
    {
        int ndom;
        int nid;
    };

public:
    ElementDataRecord();
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
    int Size() const;
    void SetElementSet(RgElementSet* pg);

    void SetItemList(FEItemList* itemList, const std::vector<int>& selection) override;

protected:
    void BuildELT();

protected:
    std::vector<ELEMREF> m_ELT;
    int m_offset;
    std::vector<FELogElemData*> m_Data;
};