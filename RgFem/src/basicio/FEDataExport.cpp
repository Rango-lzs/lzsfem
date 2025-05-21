#include "FEDataExport.h"
#include "datastructure/Vector3d.h"
using namespace std;

//-----------------------------------------------------------------------------
void FEDataExport::Serialize(FEDataStream& ar)
{
	if ((m_type == PLT_VEC3F)&&(m_fmt == FMT_NODE))
	{
		vector<Vector3d>& v = *(static_cast<vector<Vector3d>*>(m_pd));

		int n = (int) v.size();
		for (int i=0; i<n; ++i) ar << v[i];
	}
	else if ((m_type == PLT_FLOAT)&&(m_fmt == FMT_REGION))
	{
		double& d = *(static_cast<double*>(m_pd));
		ar << d;
	}
	else if ((m_type == PLT_FLOAT) && (m_fmt == FMT_NODE))
	{
		vector<double>& v = *(static_cast<vector<double>*>(m_pd));

		int n = (int)v.size();
		for (int i = 0; i<n; ++i) ar << v[i];
	}
	else
	{
		assert(false);
	}
}
