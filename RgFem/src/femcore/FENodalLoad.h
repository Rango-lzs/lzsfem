#pragma once
#include "femcore/FEDofList.h"
#include "femcore/FEModelLoad.h"
#include "femcore/FEModelParam.h"
#include "femcore/FENodeDataMap.h"

//-----------------------------------------------------------------------------
class FENodeSet;

//-----------------------------------------------------------------------------
//! Nodal load boundary condition
class FEM_EXPORT FENodalLoad : public FEModelLoad
{
    DECLARE_META_CLASS(FENodalLoad, FEModelLoad);

public:
    //! constructor
    FENodalLoad(FEModel* pfem);

    //! initialization
    bool Init() override;

    //! activation
    void Activate() override;

    //! Get the DOF list
    const FEDofList& GetDOFList() const;

    //! add a node set
    void SetNodeSet(FENodeSet* ns);

    //! get the nodeset
    FENodeSet* GetNodeSet();

    //! serialiation
    void Serialize(DumpStream& ar) override;

public:
    //! Get the DOF list
    //! This must be implemented by derived classes.
    virtual bool SetDofList(FEDofList& dofList) = 0;

    //! Get the nodal value
    //! This must be implemented by derived classes.
    //! The vals array will have the same size as the dof list.
    virtual void GetNodalValues(int inode, std::vector<double>& vals) = 0;

public:
    //! evaluate the contribution to the residual
    virtual void LoadVector(FEGlobalVector& R) override;

    //! evaluate the contribution to the global stiffness matrix
    virtual void StiffnessMatrix(FELinearSystem& LS) override;

private:
    FEDofList m_dofs;
    FENodeSet* m_nodeSet;
    bool m_brelative;
    std::vector<std::vector<double>> m_rval;

    DECLARE_PARAM_LIST();
};

//-----------------------------------------------------------------------------
// Class for prescribing the "load" on a degree of freedom.
class FEM_EXPORT FENodalDOFLoad : public FENodalLoad
{
public:
    FENodalDOFLoad(FEModel* fem);

    //! Set the DOF list
    bool SetDofList(FEDofList& dofList) override;

    //! get/set degree of freedom
    void SetDOF(int ndof)
    {
        m_dof = ndof;
    }
    int GetDOF() const
    {
        return m_dof;
    }

    //! get/set load
    void SetLoad(double s);

    void GetNodalValues(int n, std::vector<double>& val) override;

    double NodeValue(int n);

    void SetDtScale(bool b);

private:
    int m_dof;              //!< degree of freedom index
    FEParamDouble m_scale;  //!< applied load scale factor

    bool m_dtscale;

    DECLARE_PARAM_LIST();
};
