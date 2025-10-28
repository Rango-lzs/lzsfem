#include "elements/FEElementState.h"
#include "materials/FEMaterialPoint.h"


FEElementState::FEElementState()
{
}

FEElementState::~FEElementState()
{
    Clear();
}

FEElementState::FEElementState(const FEElementState& s)
{
    m_data.resize(s.m_data.size());
    for (size_t i = 0; i < m_data.size(); ++i)
    {
        if (s.m_data[i])
            m_data[i] = s.m_data[i]->Copy();
        else
            m_data[i] = 0;
    }
}

FEElementState& FEElementState::operator=(const FEElementState& s)
{
    Clear();
    m_data.resize(s.m_data.size());
    for (size_t i = 0; i < m_data.size(); ++i)
    {
        if (s.m_data[i])
            m_data[i] = s.m_data[i]->Copy();
        else
            m_data[i] = 0;
    }
    return (*this);
}

const std::vector<FEMaterialPoint*>& FEElementState::getMatPoints() const
{
    return m_data;
}

FEMaterialPoint*& FEElementState::operator[](int i)
{
    return m_data[i];
}


void FEElementState::Clear()
{
    for (size_t i = 0; i < m_data.size(); ++i)
        delete m_data[i];
    m_data.clear();
}

//! create
void FEElementState::Create(int n)
{
    m_data.assign(n, static_cast<FEMaterialPoint*>(0));
}