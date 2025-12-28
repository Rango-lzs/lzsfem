#pragma once

#include "femcore/fecore_api.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/tens4d.h"
#include "FEMaterial.h"
#include "femcore/FETimeInfo.h"
#include "RgMaterial.h"
#include "RgMaterialPoint.h"

#include <vector>


//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.
//!

class FEM_EXPORT RgMaterialPointData
{
public:
    RgMaterialPointData(RgMaterialPointData* pParant = nullptr);
    virtual ~RgMaterialPointData();

public:
    //! The init function is used to intialize data
    virtual void init();

    //! The Update function is used to update material point data
    //! Note that this gets called at the start of the time step during PreSolveUpdate
    virtual void update(const FETimeInfo& timeInfo);

    //! copy material point data (for running restarts) \todo Is this still used?
    virtual RgMaterialPointData* copy()
    {
        return 0;
    }

    // serialization
    virtual void serialize(DumpStream& ar);

public:
    //! Get the child material point data
    const std::vector<RgMaterialPointData*>& child()
    {
        return m_child;
    }

    //! Get the parent material point data
    RgMaterialPointData* parent()
    {
        return m_parent;
    }

    // assign the parent pointer
    void setParent(RgMaterialPointData* pt);

    //! add a child material point
    void addChild(RgMaterialPointData* pt);

public:
    //! Extract data (\todo Is it safe for a plugin to use this function?)
    template <class T>
    T* ExtractData();

    template <class T>
    const T* ExtractData() const;

protected:
    RgMaterialPointData* m_parent;
    std::vector<RgMaterialPointData*> m_child; //这种应该不需要

    friend class RgMaterialPoint;
};

//-----------------------------------------------------------------------------
template <class T>
inline T* RgMaterialPointData::ExtractData()
{
    // first see if this is the correct type
    T* p = dynamic_cast<T*>(this);
    if (p)
        return p;

    // check all the child classes
    for (auto child : m_child)
    {
        p = dynamic_cast<T*>(child);
        if (p)
            return p;
    }

    // search up to parent
    RgMaterialPointData* parent = m_parent;
    while (parent)
    {
        p = dynamic_cast<T*>(parent);
        if (p)
            return p;
        parent = parent->m_parent;
    }

    // Everything has failed. Material point data can not be found
    return nullptr;
}

//-----------------------------------------------------------------------------
template <class T>
inline const T* RgMaterialPointData::ExtractData() const
{
    // first see if this is the correct type
    const T* p = dynamic_cast<const T*>(this);
    if (p)
        return p;

    // check all the child classes
    for (auto child : m_child)
    {
        p = dynamic_cast<const T*>(child);
        if (p)
            return p;
    }

    // search up to parent
    const RgMaterialPointData* parent = m_parent;
    while (parent)
    {
        p = dynamic_cast<const T*>(parent);
        if (p)
            return p;
        parent = parent->m_parent;
    }

    // Everything has failed. Material point data can not be found
    return nullptr;
}


namespace RgFem
{
    namespace LargeDef
    {
        /// Data structure for storing kinematic and stress state at a material point
        class LargeDefMaterialPointData : public RgMaterialPointData
        {
        public:
            /// Kinematic variables
            Matrix3d F;      ///< Deformation gradient (current trial configuration)
            Matrix3d Fprev;  ///< Deformation gradient at last converged step
            double J;        ///< Jacobian determinant (volume ratio) J = det(F)
            Matrix3d C;      ///< Right Cauchy-Green tensor C = F^T * F
            Matrix3d E;      ///< Green-Lagrange strain (material measure)
            Matrix3d e;      ///< Almansi or small strain (spatial measure)

            /// Stress measures
            Matrix3d S;      ///< Second Piola-Kirchhoff stress (material measure)
            Matrix3d sigma;  ///< Cauchy stress (spatial)

            /// Internal variables for history-dependent materials
            std::vector<double> ivar_committed;  ///< Committed internal variables
            std::vector<double> ivar_trial;      ///< Trial internal variables

        public:
            LargeDefMaterialPointData(RgMaterialPointData* ppt = nullptr);
            virtual ~LargeDefMaterialPointData() = default;

            /// Update kinematic variables from deformation gradient
            void updateKinematicsFromF(const Matrix3d& F_new);

            /// Push forward stress from material to spatial configuration
            void pushForwardStress();

            /// Commit state variables
            void commit();

            /// Revert state variables to last committed state
            void revert();

            /// Initialize material point data
            void init() override;

            /// Update material point data
            void update(const FETimeInfo& timeInfo) override;

            /// Serialize material point data
            void serialize(DumpStream& ar) override;
        };
    }  // namespace LargeDef


    namespace SmallDef
    {
        /// Data structure for storing kinematic and stress state at a material point
        class SmallDefMaterialPointData : public RgMaterialPointData
        {
        public:
            /// Kinematic variables
            Matrix3d gradU;        ///< Displacement gradient tensor
            Matrix3ds strain;       ///< Infinitesimal strain tensor (Cauchy strain) ε = 1/2(∇u + ∇u^T)
            Matrix3ds strain_prev;  ///< Strain tensor at last converged step

            /// Stress measures
            Matrix3ds stress;       ///< Cauchy stress tensor (spatial measure)
            Matrix3ds stress_prev;  ///< Stress tensor at last converged step

            /// Internal variables for history-dependent materials
            std::vector<double> ivar_committed;  ///< Committed internal variables
            std::vector<double> ivar_trial;      ///< Trial internal variables

            /// Material tangent stiffness
            tens4d C;            ///< Fourth-order elasticity tensor (material tangent stiffness)

        public:
            SmallDefMaterialPointData(RgMaterialPointData* ppt = nullptr);
            virtual ~SmallDefMaterialPointData() = default;

            /// Update kinematic variables from displacement gradient
            void updateKinematicsFromGradU(const Matrix3d& gradU_new);

            /// Commit state variables
            void commit();

            /// Revert state variables to last committed state
            void revert();

            /// Initialize material point data
            void init() override;

            /// Update material point data
            void update(const FETimeInfo& timeInfo) override;

            /// Serialize material point data
            void serialize(DumpStream& ar) override;
        };
    }  // namespace SmallDef

}  // namespace RgFem