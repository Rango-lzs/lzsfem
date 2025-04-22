#pragma once
#include "FEDataArray.h"

class FEMaterialPoint;
class FEItemList;

//-----------------------------------------------------------------------------
// Base class for all data maps. A data map needs to be able to evaluate data across a domain
// TODO: This is a work in progress. 
// This class was added to create a base for FESurfaceMap and FEDomainMap so that both could be used in 
// FEMappedValue. 
class FEM_EXPORT FEDataMap : public FEDataArray
{
public:
	FEDataMap(FEDataMapType mapType, FEDataType dataType = FE_INVALID_TYPE);
	FEDataMap(const FEDataMap& map);

	//! set the name
	void SetName(const std::string& name);

	//! get the name
	const std::string& GetName() const;

public:
	// This function needs to be overridden by derived classes
	virtual double value(const FEMaterialPoint& mp) = 0;
	virtual Vector3d valueVec3d(const FEMaterialPoint& mp) = 0;
	virtual Matrix3d valueMat3d(const FEMaterialPoint& mp) = 0;
	virtual Matrix3ds valueMat3ds(const FEMaterialPoint& mp) = 0;

	// return the item list associated with this map
	virtual FEItemList* GetItemList() = 0;

public:
	void Serialize(DumpStream& ar) override;

protected:
	std::string	m_name;					// name of data map TODO: Move to base class?
};
