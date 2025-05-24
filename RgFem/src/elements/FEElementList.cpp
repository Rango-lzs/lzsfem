
#include "FEElementList.h"
#include "femcore/FEMesh.h"
#include "femcore/DOmain/FEDomain.h"

FEElement& FEElementList::iterator::operator*()
{ 
	return m_pmesh->Domain(m_ndom).ElementRef(m_nel); 
}

FEElement* FEElementList::iterator::operator->()
{
	return &m_pmesh->Domain(m_ndom).ElementRef(m_nel);
}

FEM_EXPORT FEElementList::iterator::operator FEElement* ()
{
	return &m_pmesh->Domain(m_ndom).ElementRef(m_nel);
}

void FEElementList::iterator::operator ++ ()
{
	if (m_pmesh && (m_ndom >= 0) && (m_nel >= 0))
	{
		m_nel++;
		if (m_nel >= m_pmesh->Domain(m_ndom).Elements())
		{
			m_ndom++;
			m_nel = 0;
			if (m_ndom >= m_pmesh->Domains())
			{
				m_ndom = -1;
				m_nel = -1;
			}
		}
	}
}
