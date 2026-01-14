#include "elements/RgElement/RgStructureElement.h"
#include "basicio/DumpStream.h"
#include "materials/FEMaterialPoint.h"
#include <math.h>



double* RgStructureElement::H(int order, int n)
{
	if (order == -1) return m_pTraits->m_H[n];
	else return m_pTraits->m_Hp[order][n];
}

int RgStructureElement::ShapeFunctions(int order) const
{
	return (order == -1 ? NodeSize() : m_pTraits->ShapeFunctions(order));
}

int RgStructureElement::GaussPointSize() const {
    return m_pTraits->m_nint;
}

