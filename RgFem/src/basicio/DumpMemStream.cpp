#include "basicio/DumpMemStream.h"
#include <assert.h>
#include <memory.h>

//-----------------------------------------------------------------------------
DumpMemStream::DumpMemStream(FEModel& fem) : DumpStream(fem)
{
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;

	Open(true, true);
}

//-----------------------------------------------------------------------------
void DumpMemStream::clear()
{
	delete [] m_pb;
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;

	// Since we can't read from an empty stream
	// we restore write mode.
	Open(true, true);
}

//-----------------------------------------------------------------------------
void DumpMemStream::Open(bool bsave, bool bshallow)
{
	DumpStream::Open(bsave, bshallow);
	if (m_pb) set_position(0);
}

//-----------------------------------------------------------------------------
bool DumpMemStream::EndOfStream() const
{
	return (bytesSerialized() >= m_nsize);
}

//-----------------------------------------------------------------------------
DumpMemStream::~DumpMemStream()
{
	clear();
}

//-----------------------------------------------------------------------------
void DumpMemStream::set_position(size_t l)
{
	assert((l >= 0) && (l < m_nreserved));
	m_pd = m_pb + l;
}

//-----------------------------------------------------------------------------
void DumpMemStream::grow_buffer(size_t l)
{
	if (l <= 0) return;

	char* pnew = new char[m_nreserved + l];
	if (m_pb)
	{
		memcpy(pnew, m_pb, m_nreserved);
		delete [] m_pb;
	}
	m_pb = pnew;
	m_pd = m_pb + m_nsize;
	m_nreserved += l;
}

//-----------------------------------------------------------------------------
size_t DumpMemStream::write(const void* pd, size_t size, size_t count)
{
	assert(IsSaving());
	size_t nsize = count*size;
	size_t lpos = (size_t)(m_pd - m_pb);
	if (lpos + nsize > m_nreserved) grow_buffer(nsize + 3*m_nreserved/2);
	memcpy(m_pd, pd, nsize);

	m_pd += nsize;
	lpos += nsize;
	if (lpos > m_nsize) m_nsize = lpos;

	return nsize;
}

//-----------------------------------------------------------------------------
size_t DumpMemStream::read(void* pd, size_t size, size_t count)
{
	assert(IsSaving()==false);
	size_t nsize = count*size;
	memcpy(pd, m_pd, nsize);
	m_pd += nsize;
	return nsize;
}
