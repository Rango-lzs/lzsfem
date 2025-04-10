#include "FEBioControlSection4.h"

#include "femcore/FEAnalysis/FEAnalysis.h"
#include "femcore/FEModel.h"

//-----------------------------------------------------------------------------
FEBioControlSection4::FEBioControlSection4(FEFileImport* pim)
    : FEFileSection(pim)
{
}

//-----------------------------------------------------------------------------
void FEBioControlSection4::Parse(XMLTag& tag)
{
    // get the step (don't allocate solver)
    FEAnalysis* pstep = GetBuilder()->GetStep(false);
    if (pstep == 0)
    {
        throw XMLReader::InvalidTag(tag);
    }

    // read the step parameters
    ReadParameterList(tag, pstep);
}
