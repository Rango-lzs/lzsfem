#include "femcore/fem_export.h"
#include <vector>
#include <string>

class DumpStream;

class FEM_EXPORT FEPlotVariable
{
public:
	FEPlotVariable();
	FEPlotVariable(const FEPlotVariable& pv);
	void operator = (const FEPlotVariable& pv);

	FEPlotVariable(const std::string& var, std::vector<int>& item, const char* szdom = "");

	void Serialize(DumpStream& ar);

	const std::string& Name() const { return m_svar; }
	const std::string& DomainName() const { return m_sdom; }

public:
	std::string			m_svar;		//!< name of output variable
	std::string			m_sdom;		//!< (optional) name of domain
	std::vector<int>	m_item;		//!< (optional) list of items
};

class FEM_EXPORT FEPlotDataStore
{
public:
	FEPlotDataStore();
	FEPlotDataStore(const FEPlotDataStore&);
	void operator = (const FEPlotDataStore&);

	void AddPlotVariable(const char* szvar, std::vector<int>& item, const char* szdom = "");

	int GetPlotCompression() const;
	void SetPlotCompression(int n);

	void SetPlotFileType(const std::string& fileType);

	void Serialize(DumpStream& ar);

	int PlotVariables() const { return (int)m_plot.size(); }
	FEPlotVariable& GetPlotVariable(int n) { return m_plot[n]; }

private:
	std::string					m_splot_type;
	std::vector<FEPlotVariable>	m_plot;
	int							m_nplot_compression;
};
