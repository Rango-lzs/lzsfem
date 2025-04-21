#include "DumpFile.h"

DumpFile::DumpFile(FEModel& fem) : DumpStream(fem)
{
	m_fp = 0;
	m_size = 0;
}

DumpFile::~DumpFile()
{
	Close();
}

bool DumpFile::Open(const char* szfile)
{
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;

	DumpStream::Open(false, false);

	return true;
}

bool DumpFile::Create(const char* szfile)
{
	m_fp = fopen(szfile, "wb");
	if (m_fp == 0) return false;

	DumpStream::Open(true, false);

	return true;
}

bool DumpFile::Append(const char* szfile)
{
	m_fp = fopen(szfile, "a+b");
	if (m_fp == 0) return false;

	DumpStream::Open(true, false);

	return true;
}

void DumpFile::Close()
{
	if (m_fp) fclose(m_fp); 
	m_fp = 0;
}

//! write buffer to archive
size_t DumpFile::write(const void* pd, size_t size, size_t count)
{
	assert(IsSaving());
	size_t elemsWritten = fwrite(pd, size, count, m_fp);
	m_size += size * elemsWritten;
	return size * elemsWritten;
}

//! read buffer from archive
size_t DumpFile::read(void* pd, size_t size, size_t count)
{
	assert(IsLoading());
	size_t elemsRead = fread(pd, size, count, m_fp);
	return size * elemsRead;
}

bool DumpFile::EndOfStream() const
{
	return (feof(m_fp) != 0);
}
