#pragma once
#include "elements/RgElement.h"
#include "FEBModel.h"
#include "femcore/FEModel.h"
#include "femcore/FEStepComponent.h"

#include <string>

class FESolver;
class FEPointFunction;
class FEMeshDataGenerator;
class FEDomainMap;
class FENodalLoad;
class FEEdgeLoad;
class FESurfaceLoad;
class FEBodyLoad;
class FEDomain;
class FESurface;
class FEFacetSet;

// This is a helper class for building the FEModel from file input.
class FEM_EXPORT FEModelBuilder
{
public:
    struct ELEMENT
    {
        int nid;
        std::vector<int> nodes;
    };

    struct NodeSetPair
    {
        char szname[256];
        FENodeSet* set1;
        FENodeSet* set2;
    };

    struct NodeSetSet
    {
        NodeSetSet()
        {
            count = 0;
        }
        enum
        {
            MAX_SETS = 32
        };
        char szname[256];
        FENodeSet* set[MAX_SETS];
        int count;

        void add(FENodeSet* ps)
        {
            set[count++] = ps;
        }
    };

    struct MappedParameter
    {
        FEParam* pp;
        FEParamObject* pc;
        const char* szname;
        int index;
    };

    struct MapLCToFunction
    {
        int lc;
        double scale;
        FEPointFunction* pf;
    };

public:
    //! constructor
    FEModelBuilder(FEModel& fem);
    virtual ~FEModelBuilder();

    //! set the active module
    void SetActiveModule(const std::string& moduleName);

    //! Get the module name
    std::string GetModuleName() const;

    // create a new analysis step
    FEAnalysis* CreateNewStep(bool allocSolver = true);

    // create a material
    FEMaterial* CreateMaterial(const char* sztype);

    // get the current step (will create a new one if no step was defined yet)
    FEAnalysis* GetStep(bool allocSolver = true);

    // add component to current step
    void AddComponent(FEStepComponent* mc);

    // reset some data for reading next step
    void NextStep();

    //! Create a domain
    virtual FEDomain* CreateDomain(FE_Element_Spec espec, FEMaterial* mat);

    //! Get the mesh
    FEMesh& GetMesh();

    //! get the FE model
    FEModel& GetFEModel();

public:
    bool BuildSurface(FESurface& s, FEFacetSet& f, bool bnodal = false);

    // bool BuildEdge(FEEdge& s, FESegmentSet& f);

    FE_Element_Spec ElementSpec(const char* sz);

    // Call this to initialize default variables when reading older files.
    void SetDefaultVariables();

public:
    virtual void AddMaterial(FEMaterial* pmat);

    void AddBC(FEBoundaryCondition* pbc);
    void AddNodalLoad(FENodalLoad* pfc);
    //void AddEdgeLoad(FEEdgeLoad* pel);
    void AddSurfaceLoad(FESurfaceLoad* psl);
    void AddInitialCondition(FEInitialCondition* pic);
    void AddContactInterface(FESurfacePairConstraint* pci);
    void AddModelLoad(FEModelLoad* pml);
    void AddNonlinearConstraint(FENLConstraint* pnc);

    // TODO: Try to remove these
    virtual void AddRigidComponent(FEStepComponent* prc);

public:
    void AddNodeSetPair(NodeSetPair& p)
    {
        m_nsetPair.push_back(p);
    }
    NodeSetPair* FindNodeSetPair(const char* szname);

    void AddNodeSetSet(NodeSetSet& p)
    {
        m_nsetSet.push_back(p);
    }
    NodeSetSet* FindNodeSetSet(const char* szname);

    FENodeSet* FindNodeSet(const std::string& setName);

public:
    void MapLoadCurveToFunction(FEPointFunction* pf, int lc, double scale = 1.0);

protected:
    FESolver* BuildSolver(FEModel& fem);

public:
    // Build the node ID table
    void BuildNodeList();

    // find a node index from its ID
    int FindNodeFromID(int nid);

    // convert an array of nodal ID to nodal indices
    void GlobalToLocalID(int* l, int n, std::vector<int>& m);

public:
    void AddMappedParameter(FEParam* p, FEObjectBase* parent, const char* szmap, int index = 0);

    void AddMeshDataGenerator(FEMeshDataGenerator* gen, FEDomainMap* map, FEParamDouble* pp);

    // This will associate all mapped parameters to their assigned maps.
    void ApplyParameterMaps();

    void ApplyLoadcurvesToFunctions();

    bool GenerateMeshDataMaps();

    // finish the build process
    bool Finish();

    FEBModel& GetFEBModel();

    void SetDefaultSolver(const std::string& s)
    {
        m_defaultSolver = s;
    }

private:
    FEModel& m_fem;       //!< model that is being constructed
    FEAnalysis* m_pStep;  //!< pointer to current analysis step
    int m_nsteps;         //!< nr of step sections read

    FEBModel m_feb;

    std::string m_defaultSolver;  //!< default solver

public:
    int m_maxid;              //!< max element ID

    bool m_b3field_hex;       //!< three-field element flag for hex (and wedge elements)
    bool m_b3field_tet;       //!< three-field element flag for quadratic tets
    bool m_b3field_shell;     //!< three-field element flag for shells
    bool m_b3field_quad;      //!< three-field element flag for quad shells
    bool m_b3field_tri;       //!< three-field element flag for tri shells
    bool m_but4;              //!< use UT4 formulation flag
    int m_default_shell;      //!< shell formulation
    bool m_shell_norm_nodal;  //!< shell normal flag (nodal or face)
    double m_ut4_alpha;       //!< UT4 integration alpha value
    bool m_ut4_bdev;          //!< UT4 integration deviatoric formulation flag
    double m_udghex_hg;       //!< hourglass parameter for UDGhex integration
    ElementType m_nhex8;      //!< hex integration rule
    ElementType m_ntet4;      //!< tet4 integration rule
    ElementType m_ntet10;     //!< tet10 integration rule
    ElementType m_ntet15;     //!< tet15 integration rule
    ElementType m_ntet20;     //!< tet20 integration rule
    ElementType m_ntri3;      //!< tri3 integration rule
    ElementType m_ntri6;      //!< tri6 integration rule
    ElementType m_ntri7;      //!< tri7 integration rule
    ElementType m_ntri10;     //!< tri10 integration rule
    ElementType m_nquad4;     //!< quad4 integration rule
    ElementType m_nquad8;     //!< quad8 integration rule
    ElementType m_nquad9;     //!< quad9 integration rule

protected:
    std::vector<NodeSetPair> m_nsetPair;
    std::vector<NodeSetSet> m_nsetSet;
    std::vector<MappedParameter> m_mappedParams;
    std::vector<MapLCToFunction> m_lc2fnc;

protected:
    int m_node_off;                //!< node offset (i.e. lowest node ID)
    std::vector<int> m_node_list;  //!< map node ID's to their nodes.
};
