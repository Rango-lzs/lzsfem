#include "NLConstraintDataRecord.h"
#include "femcore/FEModel.h"

DEFINE_META_CLASS(FELogNLConstraintData, FELogData, "");

//-----------------------------------------------------------------------------
void NLConstraintDataRecord::SetData(const char* szexpr)
{
    /*char szcopy[MAX_STRING] = {0};
    strcpy(szcopy, szexpr);
    char* sz = szcopy, *ch;
    m_Data.clear();
    strcpy(m_szdata, szexpr);
    do
    {
        ch = strchr(sz, ';');
        if (ch) *ch++ = 0;
        FELogNLConstraintData* pdata = fecore_new<FELogNLConstraintData>(sz, GetFEModel());
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    }
    while (ch);*/
}

//-----------------------------------------------------------------------------
NLConstraintDataRecord::NLConstraintDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_NLC) {}

//-----------------------------------------------------------------------------
int NLConstraintDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
double NLConstraintDataRecord::Evaluate(int item, int ndata)
{
    FEModel* fem = GetFEModel();
    int nc = item - 1;
    if ((nc < 0) || (nc >= fem->NonlinearConstraints())) return 0;
    
	FENLConstraint& nlc = *fem->NonlinearConstraint(nc);
	return m_Data[ndata]->value(nlc);
}

//-----------------------------------------------------------------------------
void NLConstraintDataRecord::SelectAllItems()
{
    FEModel* fem = GetFEModel();
    int n = fem->NonlinearConstraints();
	m_item.resize(n);
	for (int i = 0; i<n; ++i) m_item[i] = i + 1;
}
