#pragma once
#include "basicio/DataRecord.h"

class FESurface;
class FESurfaceElement;

//-----------------------------------------------------------------------------
//! This is the base class for a face data value.
class FEM_EXPORT FELogFaceData : public FELogData
{
    DECLARE_META_CLASS(FELogFaceData, FELogData);

public:
	FELogFaceData(FEModel* fem);
	virtual ~FELogFaceData();
	virtual double value(FESurfaceElement& el) = 0;
};

//-----------------------------------------------------------------------------
//! This class records surface data
class FEM_EXPORT FaceDataRecord : public DataRecord
{
public:
	FaceDataRecord(FEModel* pfem);
	double Evaluate(int item, int ndata) override;
	bool Initialize() override;
	void SetData(const char* sz) override;
	void SelectAllItems() override;
	int Size() const override;

	void SetItemList(FEItemList* itemList, const std::vector<int>& selection) override;

private:
	FESurface*	m_surface;
	std::vector<FELogFaceData*>	m_Data;
};
