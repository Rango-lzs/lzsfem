#include "RgElasticAssembler.h"
#include "FEModel.h"
#include "FEGlobalVector.h"
#include "FELinearSystem.h"
#include "FESolidElement.h"
#include "FENode.h"
#include "tools.h"

//-----------------------------------------------------------------------------
RgElasticAssembler::RgElasticAssembler(FEModel* pfem) : RgAssembler(pfem), m_domain(nullptr)
{
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::InternalForces(FEGlobalVector& R)
{
    // Clear the global vector
    R.Zero();
    
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pse->Dofs());
            ElementInternalForce(*pse, fe);
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
            
            // add element force vector to global force vector
            R.Assemble(pse->m_node, lm, fe);
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
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pse->Dofs());
            ElementBodyForce(bf, *pse, fe);
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
            
            // add element force vector to global force vector
            R.Assemble(pse->m_node, lm, fe);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
    // Clear the global vector
    R.Zero();
    
    // Check if we have a domain
    if (!m_domain) return;
    
    // Vector to store element degrees of freedom
    std::vector<int> lm;
    
    // Loop over all elements using domain's ForEachElement method
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // get the element force vector
            std::vector<double> fe;
            fe.resize(pse->Dofs());
            
            // Note: This would require an additional virtual method in derived classes
            // to actually compute the inertial forces
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
            
            // add element force vector to global force vector
            R.Assemble(pse->m_node, lm, fe);
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
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // element stiffness matrix
            Matrix ke;
            ElementStiffness(pse->GetID(), ke);
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
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
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // element stiffness matrix
            Matrix ke;
            ke.resize(pse->Dofs(), pse->Dofs());
            // Normally body force stiffness would be computed here
            ke.zero();
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
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
    m_domain->ForEachElement([&](FEElement& el) {
        // Cast to FESolidElement
        FESolidElement* pse = dynamic_cast<FESolidElement*>(&el);
        if (pse) {
            // element mass matrix
            Matrix ke;
            ElementMassMatrix(pse->GetID(), ke, scale);
            
            // get the element's LM vector
            UnpackLM(*pse, lm);
            ke.SetIndices(lm);
            
            // assemble element matrix in global mass matrix
            LS.Assemble(ke);
        }
    });
}

//-----------------------------------------------------------------------------
void RgElasticAssembler::UnpackLM(FESolidElement& el, std::vector<int>& lm)
{
    // Default implementation - create empty LM vector
    // Derived classes should override this to provide actual DOF mapping
    lm.clear();
    
    // If we have a domain, try to use its UnpackLM method
    if (m_domain) {
        m_domain->UnpackLM(el, lm);
    }
}