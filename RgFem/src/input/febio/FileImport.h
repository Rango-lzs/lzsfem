#pragma once
#include "datastructure/Matrix3d.h"
#include "datastructure/Vector3d.h"
#include "femcore/fem_export.h"
#include "input/febio/febioxml_api.h"
#include "input/febio/FEModelBuilder.h"
#include "input/XML/XMLReader.h"

#include <map>
#include <stdio.h>

//-----------------------------------------------------------------------------
// Forward declarations
class FEModel;
class FEFileImport;

//-----------------------------------------------------------------------------
// Base class for FEBio import exceptions
// Derived classes should set the error string in their constructor
class FEM_EXPORT FEFileException
{
public:
    enum
    {
        MAX_ERR_STRING = 1024
    };

    FEFileException();

    FEFileException(const char* sz, ...);

public:
    // retrieve the error string
    const char* GetErrorString() const
    {
        return m_szerr;
    }

protected:
    // set the error string (used by derived classes)
    void SetErrorString(const char* sz, ...);

protected:
    char m_szerr[MAX_ERR_STRING];
};

//-----------------------------------------------------------------------------
// Class for handling unrecognized tags.
// Use FEFileSection::SetInvalidTagHandler to set the handler for unrecognized parameters.
class FEInvalidTagHandler
{
public:
    FEInvalidTagHandler()
    {
    }
    virtual ~FEInvalidTagHandler()
    {
    }

    virtual bool ProcessTag(XMLTag& tag)
    {
        return false;
    }
};

//-----------------------------------------------------------------------------
// The FEObsoleteParamHandler class tries to map an unrecognized tag to
// a model parameter. Obsolete parameters are added with the AddParam function.
class FEObsoleteParamHandler : public FEInvalidTagHandler
{
    struct FEObsoleteParam
    {
        const char* oldName = nullptr;
        const char* newName = nullptr;
        int paramType = FE_PARAM_INVALID;
        bool readIn = false;
        union
        {
            bool bVal;
            int iVal;
            double gVal;
        };
    };

public:
    FEObsoleteParamHandler(XMLTag& tag, FEObjectBase* pc);

    // Add an obsolete parameter.
    // The oldname is a relative path w.r.t. the XMLTag passed in the constructor.
    // the newName is a ParamString, relative to the pc parameter.
    // To mark a parameter as ignored, set newName to nullptr, and paramType to FE_PARAM_INVALID.
    void AddParam(const char* oldName, const char* newName, FEParamType paramType);

    // This will try to find an entry in the m_param list. If a match is not find, this function returns false.
    bool ProcessTag(XMLTag& tag) override;

    // This function will try to map the obsolete parameters to model parameters.
    // Or it will print a warning if the obsolete parameter will be ignored.
    virtual void MapParameters();

private:
    std::string m_root;
    FEObjectBase* m_pc;
    std::vector<FEObsoleteParam> m_param;
};

//-----------------------------------------------------------------------------
// Base class for XML sections parsers
class FEBIOXML_API FEFileSection
{
public:
    FEFileSection(FEFileImport* pim)
    {
        m_pim = pim;
    }
    virtual ~FEFileSection()
    {
    }

    virtual void Parse(XMLTag& tag) = 0;

    FEFileImport* GetFileReader()
    {
        return m_pim;
    }

    FEModel* GetFEModel();

    FEModelBuilder* GetBuilder();

    // Set the handler for unrecognized tags
    void SetInvalidTagHandler(FEInvalidTagHandler* ith);

public:
    //! read a nodal ID
    //! This assumes the node ID is defined via the "id" attribute
    int ReadNodeID(XMLTag& tag);

public:
    bool ReadParameter(XMLTag& tag, FEParameterList& pl, const char* szparam = 0, FEObjectBase* pc = 0,
                       bool parseAttributes = true);
    bool ReadParameter(XMLTag& tag, FEObjectBase* pc, const char* szparam = 0, bool parseAttributes = true);
    void ReadParameterList(XMLTag& tag, FEParameterList& pl);
    void ReadParameterList(XMLTag& tag, FEObjectBase* pc);
    void ReadAttributes(XMLTag& tag, FEObjectBase* pc);

public:
    void value(XMLTag& tag, int& n);
    void value(XMLTag& tag, double& g);
    void value(XMLTag& tag, bool& b);
    void value(XMLTag& tag, Vector3d& v);
    void value(XMLTag& tag, Matrix3d& m);
    void value(XMLTag& tag, Matrix3ds& m);
    void value(XMLTag& tag, tens3drs& m);
    void value(XMLTag& tag, char* szstr);
    int value(XMLTag& tag, int* pi, int n);
    int value(XMLTag& tag, double* pf, int n);
    void value(XMLTag& tag, std::string& v);
    void value(XMLTag& tag, std::vector<int>& v);
    void value(XMLTag& tag, std::vector<double>& v);

protected:
    bool parseEnumParam(FEParam* pp, const char* val);

private:
    FEFileImport* m_pim;
    FEInvalidTagHandler* m_ith = nullptr;
};

//-----------------------------------------------------------------------------
// class that manages file section parsers
class FEBIOXML_API FEFileSectionMap : public std::map<std::string, FEFileSection*>
{
public:
    ~FEFileSectionMap();
    void Clear();

    void Parse(XMLTag& tag);
};

//-----------------------------------------------------------------------------
//! Base class for file import classes.
//! FEBio import files are XML formatted files, where each major section (children of root) is represented
//! by an FEFileSection.
//! This class also offers a simple error reporting mechanism and manages the FILE* pointer.
//! This class also manages "xml parameters". This is a feature of FEBio files that allow users to use parameters
//! as values for xml tag. A parameter is defined by a name-value pair and referenced in the input file using the
//! $(parameter_name) syntax.

class FEM_EXPORT FEFileImport
{
public:
    //! constructor
    FEFileImport();

    //! destructor
    virtual ~FEFileImport();

    //! get the error message
    void GetErrorMessage(char* szerr);

    //! Get the current FE model that is being processed
    FEModel* GetFEModel();

    //! Get the model builder
    FEModelBuilder* GetBuilder();

    // return the file path
    const char* GetFilePath();

    // set file version
    void SetFileVerion(int nversion);

    // get file version
    int GetFileVersion() const;

    // throw exception if an unknown attribute is found
    void SetStopOnUnknownAttribute(bool b);
    bool StopOnUnknownAttribute() const;

protected:
    //! open a file
    bool Open(const char* szfile, const char* szmode);

    //! close the file
    void Close();

    //! helper function for reporting errors
    bool errf(const char* szerr, ...);

    //! parse the file
    bool ParseFile(XMLTag& tag);

protected:
    FILE* m_fp;          //!< file pointer
    char m_szfile[256];  //!< file name
    char m_szerr[256];   //!< error message
    char m_szpath[512];  //!< file path

protected:
    FEFileSectionMap m_map;
    FEModelBuilder* m_builder;
    bool m_stopOnUnknownAttribute;

private:
    int m_nversion;  // version of file
};
