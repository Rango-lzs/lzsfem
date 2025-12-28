#include "materials/SmallDeformation/RgLinearElasticMatData.h"

namespace RgFem {
namespace SmallDef {

RgLinearElasticMatData::RgLinearElasticMatData(RgMaterialPointData* ppt)
    : SmallDefMaterialPointData(ppt)
{
    // For linear elastic materials, initialization is handled by parent class
    // The base SmallDef::SmallDefMaterialPointData already initializes strain, stress, etc.
}

void RgLinearElasticMatData::init()
{
    // Call parent initialization which handles the basic fields
    SmallDefMaterialPointData::init();
    
    // For linear elastic materials, no additional initialization is required
    // as the base SmallDef::SmallDefMaterialPointData handles all necessary fields
}

void RgLinearElasticMatData::update(const FETimeInfo& timeInfo)
{
    // Call parent update
    SmallDefMaterialPointData::update(timeInfo);
    
    // For linear elastic materials, no special update behavior is needed
    // The base class handles updating of history-dependent variables
}

RgMaterialPointData* RgLinearElasticMatData::copy()
{
    // Create a new instance of this class
    RgLinearElasticMatData* newData = new RgLinearElasticMatData(*this);
    
    // The copy implementation depends on the parent's Copy implementation
    // The base RgMaterialPointData class already has a Copy method that returns nullptr by default
    return newData;
}

void RgLinearElasticMatData::serialize(DumpStream& ar)
{
    // Serialize parent data
    SmallDefMaterialPointData::serialize(ar);
    
    // For linear elastic materials, no additional serialization is needed
    // as all required data is handled by the base SmallDef::SmallDefMaterialPointData
}

} // namespace SmallDef
} // namespace RgFem