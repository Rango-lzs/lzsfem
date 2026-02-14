/*********************************************************************
 * \file   RgBoundaryCondition.h
 * \brief  Boundary condition classes for FEM
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#pragma once
#include "femcore/FEStepComponent.h"
#include "femcore/FENode.h"
#include "datastructure/Vector3d.h"
#include <string>
#include <vector>

 // Forward declarations
class FEModel;
class FENodeSet;
class FELoadController;

//-----------------------------------------------------------------------------
//! Base class for boundary conditions
class FEM_EXPORT RgBoundaryCondition : public FEStepComponent
{
    DECLARE_META_CLASS(RgBoundaryCondition, FEObjectBase);

public:
    RgBoundaryCondition() {}
    //! Constructor
    RgBoundaryCondition(FEModel* fem);

    //! Destructor
    virtual ~RgBoundaryCondition();

    //! Initialize boundary condition
    virtual bool Init() override;

    //! Update model (called during solution)
    virtual void UpdateModel();

    //! Prepare for time step
    virtual void PrepStep() {}

    //! Serialize
    virtual void Serialize(DumpStream& ar) override;

    //! Set the node set this BC applies to
    void SetNodeSet(FENodeSet* nodeSet);

    //! Get the node set
    FENodeSet* GetNodeSet() { return m_nodeSet; }
    const FENodeSet* GetNodeSet() const { return m_nodeSet; }

protected:
    FEModel* m_fem;           //!< Pointer to FE model
    FENodeSet* m_nodeSet;     //!< Node set this BC applies to
};

//-----------------------------------------------------------------------------
//! Fixed (zero) displacement boundary condition
//! Corresponds to ENCASTRE in Abaqus
class FEM_EXPORT RgFixedBC : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgFixedBC, RgBoundaryCondition);

public:
    enum FIX_DOF_Type {
        DOF_X = 0x01,      //!< Fix X displacement
        DOF_Y = 0x02,      //!< Fix Y displacement
        DOF_Z = 0x04,      //!< Fix Z displacement
        DOF_RX = 0x08,     //!< Fix X rotation
        DOF_RY = 0x10,     //!< Fix Y rotation
        DOF_RZ = 0x20,     //!< Fix Z rotation
        DOF_ALL = 0x3F,    //!< Fix all DOFs (ENCASTRE)
        DOF_TRANSLATIONS = DOF_X | DOF_Y | DOF_Z,
        DOF_ROTATIONS = DOF_RX | DOF_RY | DOF_RZ
    };

public:
    RgFixedBC() {};
    RgFixedBC(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Set which DOFs to fix
    void SetDOFs(int dofMask) { m_dofMask = dofMask; }

    //! Get DOF mask
    int GetDOFs() const { return m_dofMask; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    int m_dofMask;  //!< Bitmask indicating which DOFs to fix
};

//-----------------------------------------------------------------------------
//! Prescribed displacement boundary condition
class FEM_EXPORT RgPrescribedDisplacement : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgPrescribedDisplacement, RgBoundaryCondition);

public:
    RgPrescribedDisplacement(){}
    RgPrescribedDisplacement(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update model
    void UpdateModel() override;

    //! Prepare for time step
    void PrepStep() override;

    //! Set DOF index (0=x, 1=y, 2=z for translations)
    void SetDOF(int dof) { m_dof = dof; }

    //! Get DOF index
    int GetDOF() const { return m_dof; }

    //! Set prescribed value (scale factor)
    void SetScale(double scale) { m_scale = scale; }

    //! Get prescribed value
    double GetScale() const { return m_scale; }

    //! Set load controller
    void SetLoadController(FELoadController* lc) { m_loadController = lc; }

    //! Get load controller
    FELoadController* GetLoadController() { return m_loadController; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    int m_dof;                          //!< DOF index to prescribe
    double m_scale;                     //!< Scale factor
    FELoadController* m_loadController; //!< Optional load controller for time variation
};

//-----------------------------------------------------------------------------
//! Prescribed rotation boundary condition (for shells/beams)
class FEM_EXPORT RgPrescribedRotation : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgPrescribedRotation, RgBoundaryCondition);

public:
    RgPrescribedRotation() {}
    RgPrescribedRotation(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update model
    void UpdateModel() override;

    //! Set rotation DOF (3=rx, 4=ry, 5=rz)
    void SetDOF(int dof) { m_dof = dof; }

    //! Get rotation DOF
    int GetDOF() const { return m_dof; }

    //! Set prescribed rotation value
    void SetScale(double scale) { m_scale = scale; }

    //! Get prescribed rotation
    double GetScale() const { return m_scale; }

    //! Set load controller
    void SetLoadController(FELoadController* lc) { m_loadController = lc; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    int m_dof;                          //!< Rotation DOF index (3-5)
    double m_scale;                     //!< Rotation value (radians)
    FELoadController* m_loadController; //!< Optional load controller
};

//-----------------------------------------------------------------------------
//! Prescribed velocity boundary condition (for explicit dynamics)
class FEM_EXPORT RgPrescribedVelocity : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgPrescribedVelocity, RgBoundaryCondition);

public:
    RgPrescribedVelocity() {}
    RgPrescribedVelocity(FEModel* fem);

    bool Init() override;
    void Activate() override;
    void UpdateModel() override;

    void SetDOF(int dof) { m_dof = dof; }
    int GetDOF() const { return m_dof; }

    void SetScale(double scale) { m_scale = scale; }
    double GetScale() const { return m_scale; }

    void SetLoadController(FELoadController* lc) { m_loadController = lc; }

    void Serialize(DumpStream& ar) override;

private:
    int m_dof;
    double m_scale;
    FELoadController* m_loadController;
};

//-----------------------------------------------------------------------------
//! General prescribed DOF boundary condition
//! Can handle displacement, rotation, temperature, etc.
class FEM_EXPORT RgPrescribedDOF : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgPrescribedDOF, RgBoundaryCondition);

public:
    RgPrescribedDOF(){}
    RgPrescribedDOF(FEModel* fem);

    bool Init() override;
    void Activate() override;
    void UpdateModel() override;
    void PrepStep() override;

    //! Set DOF by index
    void SetDOF(int dof) { m_dof = dof; }
    int GetDOF() const { return m_dof; }

    //! Set DOF by name
    void SetDOFName(const std::string& dofName) { m_dofName = dofName; }
    std::string GetDOFName() const { return m_dofName; }

    //! Set prescribed value
    void SetValue(double value) { m_value = value; }
    double GetValue() const { return m_value; }

    //! Set load controller
    void SetLoadController(FELoadController* lc) { m_loadController = lc; }
    FELoadController* GetLoadController() { return m_loadController; }

    void Serialize(DumpStream& ar) override;

private:
    int m_dof;                          //!< DOF index
    std::string m_dofName;              //!< DOF name (alternative to index)
    double m_value;                     //!< Prescribed value
    FELoadController* m_loadController; //!< Load controller
};

//-----------------------------------------------------------------------------
//! Symmetry boundary condition
class FEM_EXPORT RgSymmetryBC : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgSymmetryBC, RgBoundaryCondition);

public:
    RgSymmetryBC();

    bool Init() override;
    void Activate() override;

    //! Set symmetry plane normal (automatically determines which DOF to constrain)
    void SetNormal(const Vector3d& normal) { m_normal = normal; }
    Vector3d GetNormal() const { return m_normal; }

    void Serialize(DumpStream& ar) override;

private:
    Vector3d m_normal;  //!< Symmetry plane normal
    int m_dofMask;      //!< DOFs to constrain (computed from normal)
};

//-----------------------------------------------------------------------------
//! Anti-symmetry boundary condition
class FEM_EXPORT RgAntiSymmetryBC : public RgBoundaryCondition
{
    DECLARE_META_CLASS(RgAntiSymmetryBC, RgBoundaryCondition);

public:
    RgAntiSymmetryBC();
    RgAntiSymmetryBC(FEModel* fem);

    bool Init() override;
    void Activate() override;

    void SetNormal(const Vector3d& normal) { m_normal = normal; }
    Vector3d GetNormal() const { return m_normal; }

    void Serialize(DumpStream& ar) override;

private:
    Vector3d m_normal;
    int m_dofMask;
};

//-----------------------------------------------------------------------------
//! Boundary condition factory
class FEM_EXPORT RgBCFactory
{
public:
    //! Create BC from type string
    static RgBoundaryCondition* Create(FEModel* fem, const std::string& type);

    //! Create fixed BC
    static RgFixedBC* CreateFixed(FEModel* fem, FENodeSet* nodeSet, int dofMask);

    //! Create prescribed displacement BC
    static RgPrescribedDisplacement* CreatePrescribed(
        FEModel* fem, FENodeSet* nodeSet, int dof, double value);

    //! Create encastre BC (all DOFs fixed)
    static RgFixedBC* CreateEncastre(FEModel* fem, FENodeSet* nodeSet);

    //! Create symmetry BC
    static RgSymmetryBC* CreateSymmetry(FEModel* fem, FENodeSet* nodeSet,
        const Vector3d& normal);
};