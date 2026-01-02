#include "RgElasticAssembler.h"
#include "femcore/FEModel.h"
#include "femcore/FEGlobalVector.h"
#include "femcore/FELinearSystem.h"
#include "elements/RgElement/RgSolidElement.h"
#include "femcore/FENode.h"
#include "femcore/tools.h"
#include "femcore/Matrix/FEGlobalMatrix.h"


//-----------------------------------------------------------------------------
RgElasticAssembler::RgElasticAssembler(FEModel* pfem) : RgAssembler(pfem), m_domain(nullptr)
{
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::InternalForces(FEGlobalVector& R)
{
    // Clear the global vector
    //R.Zero();
    
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pElem->dofs());
            ElementInternalForce(*pElem, fe);
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            
            // add element force vector to global force vector
            R.Assemble(pElem->getNodeIds(), lm, fe);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pElem->dofs());
            ElementBodyForce(bf, *pElem, fe);
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            
            // add element force vector to global force vector
            R.Assemble(pElem->getNodeIds(), lm, fe);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
    // Clear the global vector
    //R.Zero();
    
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pElem->dofs());
            
            // Note: This would require an additional virtual method in derived classes
            // to actually compute the inertial forces
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            
            // add element force vector to global force vector
            R.Assemble(pElem->getNodeIds(), lm, fe);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::StiffnessMatrix(FELinearSystem& LS)
{
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // element stiffness matrix
            FEElementMatrix ke;
            ElementStiffness(pElem->getId(), ke);
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global stiffness matrix
            LS.Assemble(ke);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
    // For most cases, body force stiffness is zero or negligible
    // but we provide the framework for it
    
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // element stiffness matrix
            FEElementMatrix ke;
            ke.resize(pElem->dofs(), pElem->dofs());
            // Normally body force stiffness would be computed here
            ke.zero();
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global stiffness matrix
            LS.Assemble(ke);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::MassMatrix(FELinearSystem& LS, double scale)
{
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](RgElement& el) {
        // Cast to RgSolidElement
        RgSolidElement* pElem = dynamic_cast<RgSolidElement*>(&el);
        if (pElem) {
            // element mass matrix
            FEElementMatrix ke;
            ElementMassMatrix(pElem->getId(), ke, scale);
            
            // get the element's LM vector
            UnpackLM(*pElem, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global mass matrix
            LS.Assemble(ke);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::UnpackLM(RgSolidElement& el, std::vector<int>& lm)
{
    // Default implementation - create empty LM vector
    // Derived classes should override this to provide actual DOF mapping
    lm.clear();
    
    // If we have a domain, try to use its UnpackLM method
    if (m_domain) {
        m_domain->UnpackLM(el, lm);
    }
}