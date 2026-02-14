#include "RgSolidDomain.h"
#include "femcore/FEMesh.h"
#include "RgAssembler.h"
#include "elements/RgElement/RgHex8Element.h"
#include "elements/RgElement/RgTet4Element.h"
#include "elements/RgElement/RgTet10Element.h"
#include "logger/log.h"

DEFINE_META_CLASS(RgSolidDomain, RgDomain, "");

//-----------------------------------------------------------------------------
RgSolidDomain::RgSolidDomain(FEModel* fem)
    : RgDomain(fem)
{
    m_assembler = nullptr;
    m_dof_set = RgDofSet::Solid3D();
}

//-----------------------------------------------------------------------------
RgSolidDomain::RgSolidDomain()
    : RgDomain(nullptr)
{
    m_assembler = nullptr;
    m_dof_set = RgDofSet::Solid3D();
}

//-----------------------------------------------------------------------------
RgSolidDomain::~RgSolidDomain()
{
    if (m_assembler)
        delete m_assembler;

    // Clean up elements
    for (size_t i = 0; i < m_elems.size(); ++i)
    {
        delete m_elems[i];
    }
    m_elems.clear();
}

//-----------------------------------------------------------------------------
bool RgSolidDomain::Create(int nsize, ElementType type)
{
    // Allocate elements based on element specification
    m_elems.resize(nsize);

    for (int i = 0; i < nsize; ++i)
    {
        // Create appropriate element type based on espec
        // This would depend on your element types (HEX8, TET4, etc.)
        m_elems[i] = CreateSolidElement(type);
        if (m_elems[i] == nullptr)
        {
            feLogError("Failed to create solid element %d", i);
            return false;
        }

        // Set element ID
        m_elems[i]->setId(i + 1);

        // Set domain reference
        m_elems[i]->setDomain(this);
    }

    return true;
}

//-----------------------------------------------------------------------------
RgSolidElement* RgSolidDomain::CreateSolidElement(ElementType type)
{
    // Factory method to create appropriate element type
    // You'll need to implement this based on your element types
    RgSolidElement* elem = nullptr;

    switch (type)
    {
        case FE_HEX8G8:
             elem = new RgHex8Element();
            break;
        case FE_TET4G1:
             elem = new RgTet4Element();
            break;
        case FE_TET10G4:
             elem = new RgTet10Element();
            break;
            // Add other element types as needed
        default:
            return nullptr;
    }

    return elem;
}

//-----------------------------------------------------------------------------
int RgSolidDomain::Elements() const
{
    return (int)m_elems.size();
}

//-----------------------------------------------------------------------------
RgElement& RgSolidDomain::ElementRef(int n)
{
    assert(n >= 0 && n < (int)m_elems.size());
    return *m_elems[n];
}

//-----------------------------------------------------------------------------
const RgElement& RgSolidDomain::ElementRef(int n) const
{
    assert(n >= 0 && n < (int)m_elems.size());
    return *m_elems[n];
}


//-----------------------------------------------------------------------------
void RgSolidDomain::ForEachElement(std::function<void(RgElement& el)> f)
{
    int NE = Elements();
    for (int i = 0; i < NE; ++i)
    {
        f(ElementRef(i));
    }
}

//-----------------------------------------------------------------------------
void RgSolidDomain::SetMaterial(RgMaterial* pmat)
{
    RgDomain::SetMaterial(pmat);

    // Set material for all elements
    /*for (int i = 0; i < Elements(); ++i)
    {
        ElementRef(i).SetMaterial(pmat);
    }*/
}

//-----------------------------------------------------------------------------
bool RgSolidDomain::Init()
{
    // Call base class init
    if (RgDomain::Init() == false)
        return false;

    // Initialize all elements
    for (int i = 0; i < Elements(); ++i)
    {
        RgElement& el = ElementRef(i);
        /* if (el.Init() == false)
         {
             feLogError("Failed to initialize element %d in domain %s", i, GetName().c_str());
             return false;
         }*/
    }

    // Initialize the assembler if needed
    if (m_assembler)
    {
        /*if (m_assembler->Init() == false)
            return false;*/
    }

    return true;
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Reset()
{
    RgDomain::Reset();

    // Reset all elements
    for (int i = 0; i < Elements(); ++i)
    {
        //ElementRef(i).Reset();
    }
}

//-----------------------------------------------------------------------------
void RgSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    RgDomain::PreSolveUpdate(timeInfo);

    // Update all elements
    for (int i = 0; i < Elements(); ++i)
    {
        //ElementRef(i).PreSolveUpdate(timeInfo);
    }
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Activate()
{
    RgDomain::Activate();

    // Activate all domain nodes
    for (int i = 0; i < Nodes(); ++i)
    {
        FENode& node = Node(i);
        // Activate displacement DOFs (typically DOF 0, 1, 2 for x, y, z)
        node.set_active(0);
        node.set_active(1);
        node.set_active(2);
    }
}

//-----------------------------------------------------------------------------
RgMaterial* RgSolidDomain::GetMaterial()
{
    return m_pMat;
}

//-----------------------------------------------------------------------------
void RgSolidDomain::UnpackLM(RgElement& el, std::vector<int>& lm)
{
    int neln = el.NodeSize();
    int ndof = 3;  // For solid mechanics, typically 3 DOFs per node (x, y, z displacement)

    lm.resize(neln * ndof);

    for (int i = 0; i < neln; ++i)
    {
        int nid = el.getNodeId(i);
        FENode& node = m_pMesh->Node(nid);

        for (int j = 0; j < ndof; ++j)
        {
            lm[i * ndof + j] = node.getDofs()[j];
        }
    }
}

//-----------------------------------------------------------------------------
void RgSolidDomain::Update(const FETimeInfo& tp)
{
    // Update element stresses and state variables
    for (int i = 0; i < Elements(); ++i)
    {
       /* RgSolidElement& el = ElementRef(i);
        el.Update(tp);*/
    }
}

//-----------------------------------------------------------------------------
void RgSolidDomain::InternalForces(FEGlobalVector& R)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InternalForces(R);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::BodyForce(FEGlobalVector& R, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForce(R, bf);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->InertialForces(R, F);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->StiffnessMatrix(LS);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->MassMatrix(LS, scale);
		return;
	}
}

//-----------------------------------------------------------------------------
void RgSolidDomain::BodyForceStiffness(FELinearSystem& LS, FEBodyForce& bf)
{
	// use the assembler if available
	if (m_assembler) 
	{
		m_assembler->BodyForceStiffness(LS, bf);
		return;
	}
}