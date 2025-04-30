#pragma once
#include "femcore/DataRecord.h"
#include "ElementDataRecord.h"

class FEDomain;

//-----------------------------------------------------------------------------
//! Base class for domain log data
class FEM_EXPORT FELogDomainData : public FELogData
{
    DECLARE_META_CLASS(FELogDomainData, FELogData);

public:
    FELogDomainData(FEModel* fem) : FELogData(fem) {}
    virtual ~FELogDomainData() {}
    virtual double value(FEDomain& rc) = 0;

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
    FELogAvgDomainData(FEModel* pfem);
    ~FELogAvgDomainData();
    double value(FEDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELogPctDomainData : public FELogDomainData
{
public:
    FELogPctDomainData(FEModel* pfem);
    ~FELogPctDomainData();
    double value(FEDomain& rc) override;

    bool SetParameters(std::vector<std::string>& params);

private:
    double          m_pct;
    FELogElemData* m_elemData;
};

//-----------------------------------------------------------------------------
class FEM_EXPORT FELogIntegralDomainData : public FELogDomainData
{
public:
	FELogIntegralDomainData(FEModel* pfem);
	~FELogIntegralDomainData();
	double value(FEDomain& rc) override;

	bool SetParameters(std::vector<std::string>& params);

private:
	FELogElemData* m_elemData;
};
