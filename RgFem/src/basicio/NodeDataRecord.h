#pragma once
#include "femcore/DataRecord.h"

class FENodeSet;
class FENode;

//-----------------------------------------------------------------------------
//! This is the base class for a node data value.
class FEM_EXPORT FELogNodeData : public FELogData
{
    DECLARE_META_CLASS(FELogNodeData, FELogData);

public:
    FELogNodeData(FEModel* fem);
    virtual ~FELogNodeData();
    virtual double value(const FENode& node) = 0;
};

//-----------------------------------------------------------------------------
//! This class records nodal data
//! \todo should I create a different class for each data record? Like for the plot file?
class FEM_EXPORT NodeDataRecord : public DataRecord
{
public:
    NodeDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
    int Size() const;

    void SetItemList(FEItemList* items, const std::vector<int>& selection) override;

private:
    std::vector<FELogNodeData*> m_Data;
};

//-----------------------------------------------------------------------------
// Special class for outputting nodal variables
class FEM_EXPORT FENodeVarData : public FELogNodeData
{
public:
    FENodeVarData(FEModel* pfem, int ndof);
    double value(const FENode& node) override;

private:
    int m_ndof;
    DECLARE_PARAM_LIST();
};
