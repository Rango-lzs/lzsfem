#pragma once
#include "basicio/DataRecord.h"
#include "ElementDataRecord.h"

class RgDomain;

//-----------------------------------------------------------------------------
//! Base class for domain log data
class FEM_EXPORT FELogDomainData : public FELogData
{
    DECLARE_META_CLASS(FELogDomainData, FELogData);

public:
    FELogDomainData() : FELogData() {}
    virtual ~FELogDomainData() {}
    virtual double value(RgDomain& rc) = 0;

    virtual bool SetParameters(std::vector<std::string>& params) { return false; }
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FEDomainDataRecord : public DataRecord
{
public:
    FEDomainDataRecord(FEModel* pfem);
    double Evaluate(int item, int ndata);
    void SetData(const char* sz) override;
    void SelectAllItems();
    void SetDomain(int domainIndex);
    int Size() const;

private:
    std::vector<FELogDomainData*>	m_Data;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELogAvgDomainData : public FELogDomainData
{
public:
    FELogAvgDomainData();
    ~FELogAvgDomainData();
    double value(RgDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELogPctDomainData : public FELogDomainData
{
public:
    FELogPctDomainData();
    ~FELogPctDomainData();
    double value(RgDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    double          m_pct;
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELogIntegralDomainData : public FELogDomainData
{
public:
	FELogIntegralDomainData();
	~FELogIntegralDomainData();
	double value(RgDomain& rc) override;

	bool SetParameters(std::vector<std::string>& params);

private:
	FELogElemData* m_elemData;
};
