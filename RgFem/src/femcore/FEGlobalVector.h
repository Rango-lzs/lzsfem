#pragma once

#include "femcore/fem_export.h"
#include <vector>

class FEModel;

//-----------------------------------------------------------------------------
//! This class represents a global system array. It provides functions to assemble
//! local (element) vectors into this array
//! TODO: remove FEModel dependency!
class FEM_EXPORT FEGlobalVector
{
public:
    //! constructor
    FEGlobalVector(FEModel& fem, std::vector<double>& R, std::vector<double>& Fr);

    //! destructor
    virtual ~FEGlobalVector();

    //! Assemble the element vector into this global vector
    virtual void Assemble(std::vector<int>& en, std::vector<int>& elm, std::vector<double>& fe, bool bdom = false);

    //! Assemble into this global vector
    virtual void Assemble(std::vector<int>& lm, std::vector<double>& fe);

    //! assemble a nodel value
    virtual void Assemble(int node, int dof, double f);

    //! access operator
    double& operator[](int i)
    {
        return m_R[i];
    }

    //! Get the FE model
    FEModel& GetFEModel()
    {
        return m_fem;
    }

    //! get the size of the vector
    int Size() const
    {
        return (int)m_R.size();
    }

    operator std::vector<double>&()
    {
        return m_R;
    }

protected:
    FEModel& m_fem;             //!< model
    std::vector<double>& m_R;   //!< residual
    std::vector<double>& m_Fr;  //!< nodal reaction forces \todo I want to remove this
};
