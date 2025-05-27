#include "femcore/FELoadController.h"
#include "basicio/DumpStream.h"

DEFINE_META_CLASS(FELoadController, FEModelComponent, "");

FELoadController::FELoadController(FEModel* fem) : FEModelComponent(fem)
{
}

void FELoadController::Evaluate(double time)
{
	m_value = GetValue(time);
}

void FELoadController::Serialize(DumpStream& ar)
{
	FEObjectBase::Serialize(ar);
	ar & m_value;
}
