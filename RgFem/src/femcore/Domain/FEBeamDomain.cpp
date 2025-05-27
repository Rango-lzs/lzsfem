#include "FEBeamDomain.h"
#include "femcore/FEMesh.h"

DEFINE_META_CLASS(FEBeamDomain, FEDomain, "");

//-----------------------------------------------------------------------------
FEBeamDomain::FEBeamDomain(FEModel* fem) : FEDomain(FE_DOMAIN_BEAM, fem)
{

}
