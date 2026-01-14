#include "elements/RgElementState.h"
#include "materials/RgMaterialPoint.h"


RgElementState::RgElementState()
{
}

RgElementState::~RgElementState()
{
    destroy();
}

RgElementState::RgElementState(const RgElementState& s)
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

RgElementState& RgElementState::operator=(const RgElementState& s)
{
    destroy();
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

const std::vector<RgMaterialPoint*>& RgElementState::getMatPoints()
{
    return m_data;
}

int RgElementState::size() const
{
    return m_data.size();
}

RgMaterialPoint*& RgElementState::operator[](int i)
{
    return m_data[i];
}


void RgElementState::destroy()
{
    for (size_t i = 0; i < m_data.size(); ++i)
        delete m_data[i];
    m_data.clear();
}

//! create
void RgElementState::init(int n)
{
    m_data.assign(n, static_cast<RgMaterialPoint*>(0));
}