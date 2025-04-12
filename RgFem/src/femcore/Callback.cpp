#include "femcore/Callback.h"

//-----------------------------------------------------------------------------
CallbackHandler::CallbackHandler()
{
	m_event = 0;
}

//-----------------------------------------------------------------------------
CallbackHandler::~CallbackHandler()
{
}

//-----------------------------------------------------------------------------
//! This function adds a callback routine
//
void CallbackHandler::AddCallback(FECORE_CB_FNC pcb, unsigned int nwhen, void *pd, CBInsertPolicy insert)
{
	FECORE_CALLBACK cb;
	cb.m_pcb = pcb;
	cb.m_pd = pd;
	cb.m_nwhen = nwhen;

	if (insert == CBInsertPolicy::CB_ADD_END)
		m_pcb.push_back(cb);
	else
		m_pcb.push_front(cb);
}

//-----------------------------------------------------------------------------
//! Call the callback functions.
//
bool CallbackHandler::DoCallback(FEModel* fem, unsigned int nevent)
{
	m_event = nevent;
	std::list<FECORE_CALLBACK>::iterator it = m_pcb.begin();
	for (int i = 0; i<(int)m_pcb.size(); ++i, ++it)
	{
		// call the callback function
		if (it->m_nwhen & nevent)
		{
			bool bret = (it->m_pcb)(fem, nevent, it->m_pd);
			if (bret == false) return false;
		}
	}
	m_event = 0;

	return true;
}

//! Get the current callback reason (or zero if not inside DoCallback)
unsigned int CallbackHandler::CurrentEvent() const
{
	return m_event;
}
