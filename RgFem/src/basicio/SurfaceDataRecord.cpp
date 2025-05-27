#include "SurfaceDataRecord.h"
#include "femcore/FEModel.h"
#include "femcore/FESurface.h"
#include "femcore/FEMesh.h"

DEFINE_META_CLASS(FELogSurfaceData, FELogData, "");

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SetData(const char* szexpr)
{
    char szcopy[MAX_STRING] = { 0 };
    strcpy(szcopy, szexpr);
    char* sz = szcopy, * ch;
    m_Data.clear();
    strcpy(m_szdata, szexpr);
    do
    {
        ch = strchr(sz, ';');
        if (ch) *ch++ = 0;
        FELogSurfaceData* pdata = RANGO_NEW<FELogSurfaceData>(GetFEModel(),sz);
        if (pdata) m_Data.push_back(pdata);
        else throw UnknownDataField(sz);
        sz = ch;
    } while (ch);
}

//-----------------------------------------------------------------------------
FESurfaceDataRecord::FESurfaceDataRecord(FEModel* pfem) : DataRecord(pfem, FE_DATA_SURFACE) {}

//-----------------------------------------------------------------------------
int FESurfaceDataRecord::Size() const { return (int)m_Data.size(); }

//-----------------------------------------------------------------------------
double FESurfaceDataRecord::Evaluate(int item, int ndata)
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int nd = item - 1;
    if ((nd < 0) || (nd >= mesh.Surfaces())) return 0;

    FESurface& surf = mesh.Surface(nd);
    return m_Data[ndata]->value(surf);
}

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SetSurface(int surfIndex)
{
    m_item.clear();
    m_item.push_back(surfIndex + 1);
}

//-----------------------------------------------------------------------------
void FESurfaceDataRecord::SelectAllItems()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int n = mesh.Surfaces();
    m_item.resize(n);
    for (int i = 0; i < n; ++i) m_item[i] = i + 1;
}
