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
    RgMaterialPointData(RgMaterialPointData* ppt = 0);
    virtual ~RgMaterialPointData();

public:
    //! The init function is used to intialize data
    virtual void Init();

    //! The Update function is used to update material point data
    //! Note that this gets called at the start of the time step during PreSolveUpdate
    virtual void Update(const FETimeInfo& timeInfo);

    //! copy material point data (for running restarts) \todo Is this still used?
    virtual RgMaterialPointData* Copy()
    {
        return 0;
    }

    // serialization
    virtual void Serialize(DumpStream& ar);

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
    std::vector<RgMaterialPointData*> m_child;

    friend class RgMaterialPoint;
};

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
        };
    }  // namespace SmallDef

}  // namespace RgFem