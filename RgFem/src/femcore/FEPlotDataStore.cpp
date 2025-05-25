#include "FEPlotDataStore.h"
#include "basicio/DumpStream.h"

//-----------------------------------------------------------------------------
FEPlotVariable::FEPlotVariable() {}

//-----------------------------------------------------------------------------
FEPlotVariable::FEPlotVariable(const FEPlotVariable& pv)
{
    m_svar = pv.m_svar;
    m_sdom = pv.m_sdom;
    m_item = pv.m_item;
}

//-----------------------------------------------------------------------------
void FEPlotVariable::operator = (const FEPlotVariable& pv)
{
    m_svar = pv.m_svar;
    m_sdom = pv.m_sdom;
    m_item = pv.m_item;
}

FEPlotVariable::FEPlotVariable(const std::string& var, std::vector<int>& item, const char* szdom)
{
    m_svar = var;
    if (szdom) m_sdom = szdom;
    m_item = item;
}

void FEPlotVariable::Serialize(DumpStream& ar)
{
    ar & m_svar;
    ar & m_sdom;
    ar & m_item;
}

//=======================================================================================
FEPlotDataStore::FEPlotDataStore()
{
    m_plot.clear();
    m_nplot_compression = 0;
}

//-----------------------------------------------------------------------------
FEPlotDataStore::FEPlotDataStore(const FEPlotDataStore& plt)
{
    m_splot_type = plt.m_splot_type;
    m_nplot_compression = plt.m_nplot_compression;
    m_plot = plt.m_plot;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::operator = (const FEPlotDataStore& plt)
{
    m_splot_type = plt.m_splot_type;
    m_nplot_compression = plt.m_nplot_compression;
    m_plot = plt.m_plot;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::AddPlotVariable(const char* szvar, std::vector<int>& item, const char* szdom)
{
    FEPlotVariable var(szvar, item, szdom);
    m_plot.push_back(var);
}

//-----------------------------------------------------------------------------
int FEPlotDataStore::GetPlotCompression() const
{
    return m_nplot_compression;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::SetPlotCompression(int n)
{
    m_nplot_compression = n;
}

//-----------------------------------------------------------------------------
void FEPlotDataStore::SetPlotFileType(const std::string& fileType)
{
    m_splot_type = fileType;
}

void FEPlotDataStore::Serialize(DumpStream& ar)
{
    ar & m_nplot_compression;
    ar & m_splot_type;
    ar & m_plot;
}
