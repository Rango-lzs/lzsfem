#include "basicio/DataRecord.h"

class FESurface;

//-----------------------------------------------------------------------------
//! Base class for surface log data
class FEM_EXPORT FELogSurfaceData : public FELogData
{
    DECLARE_META_CLASS(FELogSurfaceData, FELogData);

public:
    FELogSurfaceData() : FELogData() {}
    virtual ~FELogSurfaceData() {}
    virtual double value(FESurface& surface) = 0;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FESurfaceDataRecord : public DataRecord
{
public:
    FESurfaceDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SetSurface(int surfIndex);
    void SelectAllItems();
    int Size() const;

private:
    std::vector<FELogSurfaceData*>	m_Data;
};
