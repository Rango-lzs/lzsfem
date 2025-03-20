/*****************************************************************/ /**
                                                                     * \file   RgNode.h
                                                                     * \brief
                                                                     *
                                                                     * \author 11914
                                                                     * \date   February 2025
                                                                     *********************************************************************/

#pragma once
#include "datastructure/Vector3d.h"
#include "FEM_EXPORT.h"
#include "femcore/DOFS.h"

#include <vector>

class DumpStream;

using GlobalDofId = int;

//-----------------------------------------------------------------------------
//! This class defines a finite element node

//! It stores nodal positions and nodal equations numbers and more.
//!
//! The m_ID array will store the equation number for the corresponding
//! degree of freedom. Its values can be (a) non-negative (0 or higher) which
//! gives the equation number in the linear system of equations, (b) -1 if the
//! dof is fixed, and (c) < -1 if the dof corresponds to a prescribed dof. In
//! that case the corresponding equation number is given by -ID-2.

class FEM_EXPORT RgNode
{
public:
    // Node status flags
    enum Status
    {
        EXCLUDE = 0x01,      // exclude node from analysis
        SHELL = 0x02,        // this node belongs to a shell
        RIGID_CLAMP = 0x04,  // this node should be clamped to a rigid body (only applies to shell nodes)
        HANGING = 0x08       // This is a hanging node
    };

public:
    //! default constructor
    RgNode();

    //! copy constructor
    RgNode(const RgNode& node);

    //! assignment operator
    RgNode& operator=(const RgNode& node);

    //! Set the number of DOFS
    void SetDOFS(int n);

    //! Get the nodal ID
    int GetID() const
    {
        return m_nID;
    }

    //! Set the node ID
    void SetID(int n)
    {
        m_nID = n;
    }

    //! see if status flags are set
    bool HasFlags(unsigned int flags) const
    {
        return ((m_nstate & flags) != 0);
    }

    //! set all the status flags
    void SetAllFlags(unsigned int flags)
    {
        m_nstate = flags;
    }

    //! get the status falgs
    unsigned int Flags() const
    {
        return m_nstate;
    }

    //! Add flags
    void SetFlags(unsigned int flags)
    {
        m_nstate |= flags;
    }

    //! Remove flags
    void UnsetFlags(unsigned int flags)
    {
        m_nstate &= ~flags;
    }

    // Serialize
    void Serialize(DumpStream& ar);

    //! Update nodal values, which copies the current values to the previous array
    void UpdateValues();


public:
    // get/set functions for current value array
    double& get(int n)
    {
        return m_val_t[n];
    }
    double get(int n) const
    {
        return m_val_t[n];
    }
    void set(int n, double v)
    {
        m_val_t[n] = v;
    }
    void add(int n, double v)
    {
        m_val_t[n] += v;
    }
    void sub(int n, double v)
    {
        m_val_t[n] -= v;
    }
    Vector3d get_Vector3d(int i, int j, int k) const
    {
        return Vector3d(m_val_t[i], m_val_t[j], m_val_t[k]);
    }
    void set_Vector3d(int i, int j, int k, const Vector3d& v)
    {
        m_val_t[i] = v.x;
        m_val_t[j] = v.y;
        m_val_t[k] = v.z;
    }

    // get functions for previous value array
    // to set these values, call UpdateValues which copies the current values
    double get_prev(int n) const
    {
        return m_val_p[n];
    }
    Vector3d get_Vector3d_prev(int i, int j, int k) const
    {
        return Vector3d(m_val_p[i], m_val_p[j], m_val_p[k]);
    }

    double get_load(int n) const
    {
        return m_Fr[n];
    }
    Vector3d get_load3(int i, int j, int k) const
    {
        return Vector3d(m_Fr[i], m_Fr[j], m_Fr[k]);
    }

    void set_load(int n, double v)
    {
        m_Fr[n] = v;
    }

public:
    // dof functions
    void set_bc(int ndof, int bcflag)
    {
        m_BC[ndof] = ((m_BC[ndof] & 0xF0) | bcflag);
    }
    void set_active(int ndof)
    {
        m_BC[ndof] |= 0x10;
    }
    void set_inactive(int ndof)
    {
        m_BC[ndof] &= 0x0F;
    }

    int get_bc(int ndof) const
    {
        return (m_BC[ndof] & 0x0F);
    }
    bool is_active(int ndof) const
    {
        return ((m_BC[ndof] & 0xF0) != 0);
    }

    int dofs() const
    {
        return (int)m_dofs.size();
    }

public:
    // return position of shell back-node
    Vector3d s0() const
    {
        return m_r0 - m_d0;
    }
    Vector3d st() const
    {
        return m_rt - m_dt;
    }
    Vector3d sp() const
    {
        return m_rp - m_dp;
    }

private:
    NodeId mId;                    //!< nodal ID

    std::vector<int> m_BC;         //!< boundary condition array, 用于标识自由度的状态，fix ，free， prescribe
    std::vector<double> m_val_t;   //!< current nodal DOF values
    std::vector<double> m_val_p;   //!< previous nodal DOF values
    std::vector<double> m_Fr;      //!< equivalent nodal forces

    std::vector<GlobalDofId> m_dofs;  //!< nodal equation numbers

    // geometry data
    Vector3d m_r0;  //!< initial position , the coordinate~
    Vector3d m_rt;  //!< current position
    Vector3d m_ra;  //!< used by rigid solver

    Vector3d m_at;  //!< nodal acceleration

    Vector3d m_rp;  //!< position of node at previous time step
    Vector3d m_vp;  //!< previous velocity
    Vector3d m_ap;  //!< previous acceleration

    //和 shell的 back node相关
    Vector3d m_d0;  //!< initial director
    Vector3d m_dt;  //!< current director
    Vector3d m_dp;  //!< director at previous time step

    // rigid body data
    unsigned int m_nstate;  //!< node state flags
    int m_rid;              //!< rigid body number
};
