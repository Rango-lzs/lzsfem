#include "FEBioGlobalsSection.h"
#include "femcore/FEGlobalData.h"
#include "femcore/FEModel.h"
#include "femcore/RTTI/MetaClass.h"

//-----------------------------------------------------------------------------
//!  This function reads the global variables from the xml file
//!
void FEBioGlobalsSection::Parse(XMLTag& tag)
{
    ++tag;
    do
    {
        if (tag == "Constants")
            ParseConstants(tag);
        else if (tag == "Solutes")
            ParseGlobalData(tag);
        else if (tag == "SolidBoundMolecules")
            ParseGlobalData(tag);
        else if (tag == "Variables")
            ParseVariables(tag);
        else
            throw XMLReader::InvalidTag(tag);
        ++tag;
    }
    while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseConstants(XMLTag& tag)
{
    FEModel& fem = *GetFEModel();
    ++tag;
    std::string s;
    double v;
    do
    {
        s = std::string(tag.Name());
        tag.value(v);
        fem.SetGlobalConstant(s, v);
        ++tag;
    }
    while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseGlobalData(XMLTag& tag)
{
    FEModel& fem = *GetFEModel();

    // read the global solute data
    ++tag;
    do
    {
        // create new global data
        FEGlobalData* pgd = static_cast<FEGlobalData*>(FEGlobalData::staic_meta()->create());
        //fecore_new<FEGlobalData>(tag.Name(), &fem);
        if (pgd == 0)
            throw XMLReader::InvalidTag(tag);
        fem.AddGlobalData(pgd);

        // TODO: We have to call the Init member here because solute data
        //       allocates the concentration dofs and they have to be allocated before
        //       materials are read in. I'd like to move this to FEModel::Init but not sure
        //       yet how.
        pgd->Init();

        // read solute properties
        ReadParameterList(tag, pgd);

        ++tag;
    }
    while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseVariables(XMLTag& tag)
{
    FEModel& fem = *GetFEModel();
    fem.GetParameterList();
    if (tag.isleaf())
        return;
    ++tag;
    do
    {
        if (tag == "var")
        {
            const char* szname = tag.AttributeValue("name");
            double v = 0;
            tag.value(v);
            fem.AddGlobalVariable(szname, v);
        }
        else
            throw XMLReader::InvalidTag(tag);
        ++tag;
    }
    while (!tag.isend());
}
