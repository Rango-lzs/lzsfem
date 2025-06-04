#include "basicio/NodeDataRecord.h"
#include "femcore/FEAnalysis/FEAnalysis.h"
#include "femcore/FEModel.h"
#include "femcore/FEMesh.h"
#include  <assert.h>

DEFINE_META_CLASS(FELogNodeData, FELogData, "");

FELogNodeData::FELogNodeData()
    : FELogData()
{
}

FELogNodeData::~FELogNodeData()
{
}

NodeDataRecord::NodeDataRecord(FEModel* pfem)
    : DataRecord(pfem, FE_DATA_NODE)
{
}

int NodeDataRecord::Size() const
{
    return (int)m_Data.size();
}

void NodeDataRecord::SetData(const char* szexpr)
{
    char szcopy[MAX_STRING] = {0};
    strcpy(szcopy, szexpr);
    char *sz = szcopy, *ch;
    m_Data.clear();
    strcpy(m_szdata, szexpr);
    FEModel* fem = GetFEModel();
    do
    {
        ch = strchr(sz, ';');
        if (ch)
            *ch++ = 0;
        FELogNodeData* pdata = RANGO_NEW<FELogNodeData>(fem, sz);
        if (pdata)
            m_Data.push_back(pdata);
        else
        {
            // see if this refers to a DOF of the model
            int ndof = fem->GetDOFIndex(sz);
            if (ndof >= 0)
            {
                // Add an output for a nodal variable
                pdata = new FENodeVarData(fem, ndof);
                m_Data.push_back(pdata);
            }
            else
                throw UnknownDataField(sz);
        }
        sz = ch;
    }
    while (ch);
}

double NodeDataRecord::Evaluate(int item, int ndata)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int nnode = item - 1;
    assert((nnode >= 0) && (nnode < mesh.Nodes()));
    if ((nnode < 0) || (nnode >= mesh.Nodes()))
        return 0;

    FENode& node = mesh.Node(nnode);
    return m_Data[ndata]->value(node);
}

void NodeDataRecord::SelectAllItems()
{
    int n = GetFEModel()->GetMesh().Nodes();
    m_item.resize(n);
    for (int i = 0; i < n; ++i)
        m_item[i] = i + 1;
}

void NodeDataRecord::SetItemList(FEItemList* items, const std::vector<int>& selection)
{
    // TODO: We don't support using a selection of a node set yet.
    assert(selection.empty());
    FENodeSet* pns = dynamic_cast<FENodeSet*>(items);
    assert(pns);
    int n = pns->Size();
    m_item.resize(n);
    for (int i = 0; i < n; ++i)
        m_item[i] = (*pns)[i] + 1;
}

FENodeVarData::FENodeVarData(FEModel* pfem, int ndof)
    : FELogNodeData()
    , m_ndof(ndof)
{
}

double FENodeVarData::value(const FENode& node)
{
    return node.get(m_ndof);
}
