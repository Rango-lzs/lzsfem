#pragma once
#include "femcore/FEObjectBase.h"
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <vector>

//-----------------------------------------------------------------------------
// forward declaration
class FEModel;
class DumpStream;
class FEItemList;

//-----------------------------------------------------------------------------
enum FEDataRecordType
{
    FE_DATA_NODE = 0x01,
    FE_DATA_FACE,
    FE_DATA_ELEM,
    FE_DATA_RB,
    FE_DATA_NLC,
    FE_DATA_SURFACE,
    FE_DATA_DOMAIN,
    FE_DATA_MODEL
};

//-----------------------------------------------------------------------------
// Exception thrown when parsing fails
class FEM_EXPORT UnknownDataField : public std::runtime_error
{
public:
    UnknownDataField(const char* sz);
};

//-----------------------------------------------------------------------------

class FEM_EXPORT DataRecord : public FEObjectBase
{
    DECLARE_META_CLASS(DataRecord, FEObjectBase);

public:
    enum
    {
        MAX_DELIM = 16,
        MAX_STRING = 1024
    };

public:
    DataRecord(FEModel* pfem, int ntype);
    virtual ~DataRecord();

    bool SetFileName(const char* szfile);

    bool Write();

    void SetItemList(const std::vector<int>& items);
    virtual void SetItemList(FEItemList* itemList, const std::vector<int>& selection);

    void SetName(const char* sz);
    void SetDelim(const char* sz);
    void SetFormat(const char* sz);
    void SetComments(bool b)
    {
        m_bcomm = b;
    }

public:
    virtual bool Initialize();
    virtual double Evaluate(int item, int ndata) = 0;
    virtual void SelectAllItems() = 0;
    virtual void Serialize(DumpStream& ar);
    virtual void SetData(const char* sz) = 0;
    virtual int Size() const = 0;

private:
    std::string printToString(int i);
    std::string printToFormatString(int i);

public:
    int m_nid;                //!< ID of data record
    std::vector<int> m_item;  //!< item list
    int m_type;               //!< type of data record

protected:
    bool m_bcomm;               //!< export comments or not
    char m_szname[MAX_STRING];  //!< name of expression
    char m_szdelim[MAX_DELIM];  //!< data delimitor
    char m_szdata[MAX_STRING];  //!< data expression
    char m_szfmt[MAX_STRING];   //!< max format string

protected:
    char m_szfile[MAX_STRING];  //!< file name of data record
    FILE* m_fp;
};

//=========================================================================
// Super class for log data classes.
class FEM_EXPORT FELogData : public FEObjectBase
{
    DECLARE_META_CLASS(FELogData, FEObjectBase);

public:
    FELogData(FEModel* fem);
};
