#pragma once

#include "RgMaterial.h"
#include "FEMaterial.h"
#include "FEMaterialPoint.h"
#include "datastructure/Matrix3d.h"
#include "datastructure/tens4d.h"
#include "femcore/FETimeInfo.h"
#include <vector>

namespace RgFem {
namespace Framework {

/// Data structure for storing kinematic and stress state at a material point
class MaterialPointData : public FEMaterialPointData {
public:
    /// Kinematic variables
    Matrix3d F;          ///< Deformation gradient (current trial configuration)
    Matrix3d Fprev;      ///< Deformation gradient at last converged step
    double J;            ///< Jacobian determinant (volume ratio) J = det(F)
    Matrix3d C;          ///< Right Cauchy-Green tensor C = F^T * F
    Matrix3d E;          ///< Green-Lagrange strain (material measure)
    Matrix3d e;          ///< Almansi or small strain (spatial measure)
    
    /// Stress measures
    Matrix3d S;          ///< Second Piola-Kirchhoff stress (material measure)
    Matrix3d sigma;      ///< Cauchy stress (spatial)

    /// Internal variables for history-dependent materials
    std::vector<double> ivar_committed;  ///< Committed internal variables
    std::vector<double> ivar_trial;      ///< Trial internal variables

public:
    MaterialPointData(FEMaterialPointData* ppt = nullptr);
    virtual ~MaterialPointData() = default;

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

} // namespace Framework
} // namespace RgFem