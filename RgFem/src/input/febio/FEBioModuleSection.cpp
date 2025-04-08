#include "FEBioModuleSection.h"
#include "femcore/FEModel.h"

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//! It is currently only used to allocate the FESolver.
void FEBioModuleSection::Parse(XMLTag& tag)
{
    int nversion = GetFileReader()->GetFileVersion();

    // For version 2.5 and up this tag can only be read in once and we
    // determine that by inspecting the m_szmod member.
    if (nversion >= 0x0205)
    {
        std::string moduleName = GetBuilder()->GetModuleName();
        if (moduleName.empty() == false)
            throw XMLReader::InvalidTag(tag);
    }

    // get the type attribute
    const char* szt = tag.AttributeValue("type");

    // some special case
    if (strcmp(szt, "explicit-solid") == 0)
    {
        szt = "solid";
        GetBuilder()->SetDefaultSolver("explicit-solid");
    }
    else if (strcmp(szt, "CG-solid") == 0)
    {
        szt = "solid";
        GetBuilder()->SetDefaultSolver("CG-solid");
    }

    GetBuilder()->SetActiveModule(szt);

    if (tag.isempty())
        return;

    ++tag;
    do
    {
        if (tag == "units")
        {
            const char* szunits = tag.szvalue();
            GetFEModel()->SetUnits(szunits);
        }

        ++tag;
    }
    while (!tag.isend());
}
