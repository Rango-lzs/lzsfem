#pragma once

class FEM_EXPORT FEElementState
{
public:
    //! default constructor
    FEElementState()
    {
    }

    //! destructor
    ~FEElementState()
    {
        Clear();
    }

    //! copy constructor
    FEElementState(const FEElementState& s);

    //! assignment operator
    FEElementState& operator=(const FEElementState& s);

    //! clear state data
    void Clear()
    {
        for (size_t i = 0; i < m_data.size(); ++i)
            delete m_data[i];
        m_data.clear();
    }

    //! create
    void Create(int n)
    {
        m_data.assign(n, static_cast<FEMaterialPoint*>(0));
    }

    //! operator for easy access to element data
    FEMaterialPoint*& operator[](int n)
    {
        return m_data[n];
    }

private:
    std::vector<FEMaterialPoint*> m_data;
};