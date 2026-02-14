#pragma once

#include "materials/RgMaterialPointData.h"
#include "datastructure/Matrix.h"
#include "datastructure/Matrix3d.h"
#include "basicio/DumpStream.h"

namespace SmallDef {

/// Material point data for linear elastic materials
/// This class extends the general small deformation material point data
/// with specific fields and methods for linear elastic behavior
class RgLinearElasticMatData : public SmallDefMaterialPointData
{
public:
    /// Constructor
    RgLinearElasticMatData(RgMaterialPointData* ppt = nullptr);

    /// Destructor
    virtual ~RgLinearElasticMatData() = default;

    /// Initialize material point data
    void init() override;

    /// Update material point data
    void update(const FETimeInfo& timeInfo) override;

    /// Copy material point data
    RgMaterialPointData* copy() override;

    /// Serialize material point data
    void serialize(DumpStream& ar) override;

public:
    // For linear elastic materials, we can reuse the basic small deformation
    // fields which include:
    // - strain: infinitesimal strain tensor
    // - stress: Cauchy stress tensor
    // These are inherited from SmallDef::SmallDefMaterialPointData
};

} // namespace SmallDef
