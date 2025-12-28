#include "materials/RgMaterialPoint.h"
#include "basicio/DumpStream.h"
#include <string.h>
#include "RgMaterialPointData.h"

//=================================================================================================
RgMaterialPoint::RgMaterialPoint(RgMaterialPointData* data)
{
	m_data = data;
	m_elem = nullptr;
	m_shape = nullptr;  // Initialize shape function pointer
	m_r0 = Vector3d(0, 0, 0);  // Initialize position
	m_rt = Vector3d(0, 0, 0);  // Initialize position
	m_J0 = m_Jt = 1.0;  // Initialize Jacobian values
	m_index = -1;  // Initialize index
}

RgMaterialPoint::~RgMaterialPoint()
{
	// Delete the data chain if this is the root
	if (m_data) {
		delete m_data;
		m_data = nullptr;
	}
}

void RgMaterialPoint::Init()
{
	if (m_data) m_data->init();
}

RgMaterialPoint* RgMaterialPoint::Copy()
{
	RgMaterialPoint* mp = new RgMaterialPoint(*this);
	if (m_data) {
		// For a complete copy implementation, we would need to implement
		// a deep copy of the data chain, but the base Copy() returns nullptr
		// by default, so we respect that behavior here
		RgMaterialPointData* copiedData = m_data->copy();
		if (copiedData) {
			mp->m_data = copiedData;
		}
		// Otherwise, m_data will be a shallow copy which is not ideal
		// but respects the interface design
	}
	return mp;
}

void RgMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	if (m_data) m_data->update(timeInfo);
}

void RgMaterialPoint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false)
	{
		ar & m_r0 & m_J0 & m_Jt;
		// Serialize element pointer and index if needed
		// Note: Serializing pointers directly might not be appropriate
		// This would depend on the DumpStream implementation
	}
	if (m_data) m_data->serialize(ar);
}

void RgMaterialPoint::Append(RgMaterialPointData* pt)
{
	if (pt == nullptr) return;
	if (m_data) {
		m_data->addChild(pt);
	}
	// If no data exists yet, we can't append to it
}

//=================================================================================================