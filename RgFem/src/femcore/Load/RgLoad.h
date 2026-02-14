/*********************************************************************
 * \file   RgLoad.h
 * \brief  Load classes for FEM
 *
 * \author Leizs
 * \date   January 2025
 *********************************************************************/

#pragma once
#include "femcore/FENode.h"
#include "datastructure/Vector3d.h"
#include <string>
#include <vector>
#include "femcore/FEStepComponent.h"

// Forward declarations
class FEModel;
class FENodeSet;
class FEFacetSet;
class FESurface;
class RgLoadController;
class RgMaterialPoint;

//-----------------------------------------------------------------------------
//! Base class for loads
class FEM_EXPORT RgLoad : public FEStepComponent
{
    DECLARE_META_CLASS(RgLoad, FEStepComponent);

public:
    RgLoad() {}
    RgLoad(FEModel* fem);

    //! Destructor
    virtual ~RgLoad();

    //! Initialize load
    virtual bool Init() override;

    //! Activate load
    virtual void Activate() override;

    //! Deactivate load
    virtual void Deactivate() override;

    //! Update load (evaluate at current time)
    virtual void Update();

    //! Get current load magnitude
    virtual double GetMagnitude() const { return m_magnitude; }

    //! Serialize
    virtual void Serialize(DumpStream& ar) override;

    //! Set load magnitude/scale
    void SetMagnitude(double mag) { m_magnitude = mag; }

    //! Set load controller
    void SetLoadController(RgLoadController* lc) { m_loadController = lc; }

    //! Get load controller
    RgLoadController* GetLoadController() { return m_loadController; }

protected:
    FEModel* m_fem;                     //!< Pointer to FE model
    double m_magnitude;                  //!< Load magnitude/scale factor
    RgLoadController* m_loadController;  //!< Load controller for time variation
};

//-----------------------------------------------------------------------------
//! Nodal force load
class FEM_EXPORT RgNodalLoad : public RgLoad
{
    DECLARE_META_CLASS(RgNodalLoad, RgLoad);

public:
    RgNodalLoad() {}
    RgNodalLoad(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update
    void Update() override;

    //! Set node set
    void SetNodeSet(FENodeSet* nodeSet) { m_nodeSet = nodeSet; }

    //! Get node set
    FENodeSet* GetNodeSet() { return m_nodeSet; }

    //! Set force vector
    void SetForce(const Vector3d& force) { m_force = force; }

    //! Get force vector
    Vector3d GetForce() const { return m_force; }

    //! Set DOF (0-2 for x,y,z forces; 3-5 for moments)
    void SetDOF(int dof) { m_dof = dof; }

    //! Get DOF
    int GetDOF() const { return m_dof; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    FENodeSet* m_nodeSet;  //!< Node set to apply load
    Vector3d m_force;      //!< Force vector
    int m_dof;             //!< DOF index (-1 for vector, 0-5 for single DOF)
};

//-----------------------------------------------------------------------------
//! Distributed surface load (pressure, traction)
class FEM_EXPORT RgSurfaceLoad : public RgLoad
{
    DECLARE_META_CLASS(RgSurfaceLoad, RgLoad);

public:
    enum LoadType {
        PRESSURE,        //!< Normal pressure
        TRACTION,        //!< Tangential traction
        GENERAL_TRACTION //!< General traction with components
    };

public:
    RgSurfaceLoad() {}
    RgSurfaceLoad(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update
    void Update() override;

    //! Set surface/facet set
    void SetSurface(FESurface* surface) { m_surface = surface; }
    void SetFacetSet(FEFacetSet* facetSet) { m_facetSet = facetSet; }

    //! Get surface
    FESurface* GetSurface() { return m_surface; }
    FEFacetSet* GetFacetSet() { return m_facetSet; }

    //! Set load type
    void SetLoadType(LoadType type) { m_loadType = type; }

    //! Get load type
    LoadType GetLoadType() const { return m_loadType; }

    //! Set traction vector (for GENERAL_TRACTION)
    void SetTraction(const Vector3d& traction) { m_traction = traction; }

    //! Get traction vector
    Vector3d GetTraction() const { return m_traction; }

    //! Set whether traction follows surface (deformed normal)
    void SetFollower(bool follower) { m_bFollower = follower; }

    //! Is follower load
    bool IsFollower() const { return m_bFollower; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    FESurface* m_surface;      //!< Surface to apply load
    FEFacetSet* m_facetSet;    //!< Facet set (alternative to surface)
    LoadType m_loadType;        //!< Type of surface load
    Vector3d m_traction;        //!< Traction vector
    bool m_bFollower;           //!< Follower load flag
};

//-----------------------------------------------------------------------------
//! Body force load (gravity, centrifugal, etc.)
class FEM_EXPORT RgBodyLoad : public RgLoad
{
    DECLARE_META_CLASS(RgBodyLoad, RgLoad);

public:
    enum BodyLoadType {
        GRAVITY,         //!< Gravitational load
        CENTRIFUGAL,     //!< Centrifugal load
        CONSTANT         //!< Constant body force
    };

public:
    RgBodyLoad() {}
    RgBodyLoad(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update
    void Update() override;

    //! Set body load type
    void SetBodyLoadType(BodyLoadType type) { m_bodyLoadType = type; }

    //! Get body load type
    BodyLoadType GetBodyLoadType() const { return m_bodyLoadType; }

    //! Set force vector (for GRAVITY or CONSTANT)
    void SetForce(const Vector3d& force) { m_force = force; }

    //! Get force vector
    Vector3d GetForce() const { return m_force; }

    //! Set rotation axis (for CENTRIFUGAL)
    void SetAxis(const Vector3d& axis) { m_axis = axis; }

    //! Get rotation axis
    Vector3d GetAxis() const { return m_axis; }

    //! Set rotation origin (for CENTRIFUGAL)
    void SetOrigin(const Vector3d& origin) { m_origin = origin; }

    //! Get rotation origin
    Vector3d GetOrigin() const { return m_origin; }

    //! Set angular velocity (for CENTRIFUGAL)
    void SetAngularVelocity(double omega) { m_omega = omega; }

    //! Get angular velocity
    double GetAngularVelocity() const { return m_omega; }

    //! Evaluate body force at a material point
    Vector3d EvaluateForce(const RgMaterialPoint& mp) const;

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    BodyLoadType m_bodyLoadType;  //!< Type of body load
    Vector3d m_force;              //!< Force vector (for GRAVITY/CONSTANT)
    Vector3d m_axis;               //!< Rotation axis (for CENTRIFUGAL)
    Vector3d m_origin;             //!< Rotation origin (for CENTRIFUGAL)
    double m_omega;                //!< Angular velocity (for CENTRIFUGAL)
};

//-----------------------------------------------------------------------------
//! Concentrated moment load
class FEM_EXPORT RgMomentLoad : public RgLoad
{
    DECLARE_META_CLASS(RgMomentLoad, RgLoad);

public:
    RgMomentLoad() {}
    RgMomentLoad(FEModel* fem);

    //! Initialize
    bool Init() override;

    //! Activate
    void Activate() override;

    //! Update
    void Update() override;

    //! Set node set
    void SetNodeSet(FENodeSet* nodeSet) { m_nodeSet = nodeSet; }

    //! Get node set
    FENodeSet* GetNodeSet() { return m_nodeSet; }

    //! Set moment vector
    void SetMoment(const Vector3d& moment) { m_moment = moment; }

    //! Get moment vector
    Vector3d GetMoment() const { return m_moment; }

    //! Serialize
    void Serialize(DumpStream& ar) override;

private:
    FENodeSet* m_nodeSet;  //!< Node set to apply moment
    Vector3d m_moment;     //!< Moment vector
};

//-----------------------------------------------------------------------------
//! Load factory
class FEM_EXPORT RgLoadFactory
{
public:
    //! Create load from type string
    static RgLoad* Create(FEModel* fem, const std::string& type);

    //! Create nodal force load
    static RgNodalLoad* CreateNodalLoad(FEModel* fem, FENodeSet* nodeSet, 
                                        const Vector3d& force);

    //! Create nodal force load (single DOF)
    static RgNodalLoad* CreateNodalLoad(FEModel* fem, FENodeSet* nodeSet, 
                                        int dof, double magnitude);

    //! Create pressure load
    static RgSurfaceLoad* CreatePressure(FEModel* fem, FESurface* surface, 
                                         double pressure);

    //! Create pressure load on facet set
    static RgSurfaceLoad* CreatePressure(FEModel* fem, FEFacetSet* facetSet, 
                                         double pressure);

    //! Create traction load
    static RgSurfaceLoad* CreateTraction(FEModel* fem, FESurface* surface, 
                                         const Vector3d& traction);

    //! Create gravity load
    static RgBodyLoad* CreateGravity(FEModel* fem, const Vector3d& g);

    //! Create centrifugal load
    static RgBodyLoad* CreateCentrifugal(FEModel* fem, const Vector3d& axis, 
                                         const Vector3d& origin, double omega);

    //! Create moment load
    static RgMomentLoad* CreateMoment(FEModel* fem, FENodeSet* nodeSet, 
                                      const Vector3d& moment);
};
