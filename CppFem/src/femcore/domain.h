/*****************************************************************//**
 * \file   domain.h
 * \brief  Description of the fem model
 * 
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef domain_h
#define domain_h

/**
 * Class and object Domain. Domain contains mesh description, or if program runs in parallel then it contains
 * description of domain associated to particular processor or thread of execution. Generally, it contain and
 * manages lists of Dof managers, elements, boundary conditions, cross sections and materials - these describe
 * the geometry of problem, its constitutive properties and applied boundary conditions. Services for accessing
 * these objects are provided. Domain is attribute of engineering model - which represent type of analysis
 * which should be performed.
 *
 * Domain also provides services for reading its description from
 * input stream and instantiating corresponding components accordingly. The basic Domain task are following
 * - Reading its description from input and creating corresponding objects.
 * - Provides services for accessing its particular components.
 * - Checking yourself.
 */

class Domain
{
public:
    /// Element list.
    std :: vector< std :: unique_ptr< Element > > elementList;
    /// Dof manager list.
    std :: vector< std :: unique_ptr< DofManager > > dofManagerList;
    /// Cross section list.
    std :: vector< std :: unique_ptr< CrossSection > > crossSectionList;    
private:
    /// Material list.
    std :: vector< std :: unique_ptr< Material > > materialList;
    /// Boundary condition list.
    std :: vector< std :: unique_ptr< GeneralBoundaryCondition > > bcList;
    /// Initial condition list.
    std :: vector< std :: unique_ptr< InitialCondition > > icList;
    /// Load time function list.
    std :: vector< std :: unique_ptr< Function > > functionList;
    /// Set list.
    std :: vector< std :: unique_ptr< Set > > setList;
    /// Nonlocal barrier list.
    std :: vector< std :: unique_ptr< NonlocalBarrier > > nonlocalBarrierList;

    /// Default dofs for a node (depends on the domain type).
    IntArray defaultNodeDofIDArry;

    /**
     * Domain type. Determined by input data. It determines the problem type (like plane stress or plane strain mode).
     * According to this mode the default number of Dofs per node (or side) and their physical meaning are determined.
     * These default settings can be redefined by particular node or side. See related documentation for details.
     * @see Node
     * @see ElementSide
     */
    domainType dType;

    /**
     * Associated Engineering model. An abstraction for type of analysis which will be prformed.
     */
    EngngModel *engineeringModel;

    /**
     * Domain connectivity table. Table is build upon request.
     * Provides connectivity information of current domain.
     */
    std :: unique_ptr< ConnectivityTable > connectivityTable;
    /**
     * Spatial Localizer. It is build upon request.
     * Provides the spatial localization services.
     */
    std :: unique_ptr< SpatialLocalizer > spatialLocalizer;
    /// Output manager, allowing to filter the produced output.
    std :: unique_ptr< OutputManager > outputManager;
    /// Domain number.
    int number;
    /// Domain serial (version) number. Used for domain version identification during Adaptive computations.
    int serialNumber;
    /// Number of spatial dimensions
    int nsd;
    bool axisymm;
    /// nodal recovery object associated to receiver.
    std :: unique_ptr< NodalRecoveryModel > smoother; ///@todo I don't see why this has to be stored, and there is only one? /Mikael

    std :: string mDomainType;
    /**
     * For nonlocal models of integral type
     * it is necessary, mainly due to resulting efficiency, to compute variable(s)
     * which are nonlocally averaged in advance, before average process begins.
     * The loop over all  integration points is typically made to compute these variables.
     * To prevent doing this multiple times at the same solution state,
     * the modification time mark is kept.
     * This state counter could not be kept in static global variable,
     * because in case of multiple domains stateCounter should be kept independently for each domain.
     */
    StateCounterType nonlocalUpdateStateCounter;
    /// XFEM Manager
    std :: unique_ptr< XfemManager > xfemManager;

    /// Fracture Manager
    std :: unique_ptr< FractureManager > fracManager;

    /// Contact Manager
    std :: unique_ptr< ContactManager > contactManager;

    /// BC tracker (keeps track of BCs applied wia sets to components)
    BCTracker bcTracker;
    
    /**
     * Map from an element's global number (label) to its place
     * in the element array. Added by ES 140326.
     */
    // elementGlobal2LocalMap
    std::unordered_map< int, int > elementGlobal2LocalMap;

    /**
     * Map from a dofmans's global number (label) to its place
     * in the dofman array.
     */
    std::unordered_map< int, int > dofmanGlobal2LocalMap;

    /**
     * Map from material number to elements that have the
     * given material number. Added by ES 140718.
     */
    std::unordered_map< int, IntArray> materialNum2ElMap;

    /// Topology description
    std :: unique_ptr< TopologyDescription > topology;

public:
    /// Keeps track of next free dof ID (for special Lagrange multipliers, XFEM and such)
    int freeDofID;
private:

#ifdef __PARALLEL_MODE
    /**
     * Transaction manager. The purpose of this class is to
     * make the domain modification (in terms of adding and deleting components) versatile.
     */
    std :: unique_ptr< DomainTransactionManager > transactionManager;
    /// Global dof manager map (index is global of man number).
    std :: map< int, DofManager * >dmanMap;
    /// dmanMap init flag.
    bool dmanMapInitialized;
    /// Global element map (index is global of man number).
    std :: map< int, Element * >elementMap;
    /// dmanMap init flag.
    bool elementMapInitialized;


    /**@name Load Balancing data structures */
    //@{
    /// List of received elements.
    std :: list< Element * >recvElemList; ///@todo This seems like dead, old unused code. Remove it? / Mikael
    //@}
#endif

public:
    /**
     * Constructor. Creates empty n-th domain belonging to given engineering model.
     * @param n Domain number.
     * @param serNum Serial number
     * @param e Engineering model domain will belong to.
     * @see giveSerialNumber
     */
    Domain(int n, int serNum, EngngModel * e);

    Domain(const Domain& src) = delete;
    Domain &operator = (const Domain &src) = delete;

    /// Destructor.
    ~Domain();

    /// Returns domain number.
    int giveNumber() { return this->number; }
    /// Returns domain number.
    void setNumber(int nn) { this->number = nn; }
    /// Returns domain serial (version) number.
    int giveSerialNumber() { return this->serialNumber; }

    // management of the mesh components
    /**
     * Service for accessing particular domain fe element.
     * Generates error if no such element is defined.
     * @param n Pointer to n-th element is returned.
     */
    Element *giveElement(int n);
    std :: vector< std :: unique_ptr< Element > > &giveElements() { return this->elementList; }
    /**
     * Service for accessing particular domain fe element.
     * Generates error if no such element is defined.
     * @param n Pointer to the element with id n
     */
    Element *giveGlobalElement(int n);
};
#endif // domain_h
