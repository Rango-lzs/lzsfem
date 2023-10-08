/*****************************************************************//**
 * \file   domain.h
 * \brief  Description of the fem model
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef domain_h
#define domain_h

#include "domaintype.h"
#include "intarray.h"

#include <unordered_map>
#include <memory>

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

#define ParameterListDeclare(ParameterListDefiner) \
class Parameter\
{\
public:\
	ParameterListDefiner\
};\
static Parameter param; \

#define ParameterDefiner(ParamName) std::string ParamName = #ParamName

#define DoaminParamListDefiner \
ParameterDefiner(name); \
ParameterDefiner(type); \

ParameterListDeclare(DoaminParamListDefiner)

///@name Input fields for domains
#define _IFT_Domain_type "domain" ///< This is trouble, will not work with dynamic input record
#define _IFT_Domain_ndofman "ndofman"
#define _IFT_Domain_nelem "nelem"
#define _IFT_Domain_nmat "nmat"
#define _IFT_Domain_ncrosssect "ncrosssect"
#define _IFT_Domain_nbc "nbc"
#define _IFT_Domain_nic "nic"
#define _IFT_Domain_nfunct "nltf"
#define _IFT_Domain_nset "nset"
#define _IFT_Domain_nbarrier "nbarrier"
#define _IFT_Domain_topology "topology"
#define _IFT_Domain_nxfemman "nxfemman" /// [in,optional] Specifies if there is an xfem-manager.
#define _IFT_Domain_ncontactman "ncontactman" /// [in,optional] Specifies if there is a contact manager.
#define _IFT_Domain_numberOfSpatialDimensions "nsd" ///< [in,optional] Specifies how many spatial dimensions the domain has.
#define _IFT_Domain_nfracman "nfracman" /// [in,optional] Specifies if there is a fracture manager.
#define _IFT_Domain_axisymmetric "axisymm" /// [optional] Specifies if the problem is axisymmetric.

namespace fem 
{
	class Element;
	class DofManager;
	class CrossSection;
	class Material;
	class GeneralBoundaryCondition;
	class InitialCondition;
	class Function;
	class Set;

	class EngngModel;
	class ConnectivityTable;
	class OutputManager;
	class TopologyDescription;

	class DataReader;

	class Domain
	{
		//ÓÃºê¶¨Òå
		class Parameter
		{
		public:
			const std::string domain_type = "domain";
			const std::string ndofman = "ndofman";
		};

		static Parameter param;

	private:
		/// Element list.
		std::vector< std::unique_ptr< Element > > m_elementList;
		/// Dof manager list.
		std::vector< std::unique_ptr< DofManager > > m_dofManagerList;
		/// Cross section list.
		std::vector< std::unique_ptr< CrossSection > > m_crossSectionList;
		/// Material list.
		std::vector< std::unique_ptr< Material > > m_materialList;
		/// Boundary condition list.
		std::vector< std::unique_ptr< GeneralBoundaryCondition > > m_bcList;
		/// Initial condition list.
		std::vector< std::unique_ptr< InitialCondition > > m_icList;
		/// Load time function list.
		std::vector< std::unique_ptr< Function > > m_functionList;
		/// Set list.
		std::vector< std::unique_ptr<Set> > m_setList;
				
		/**
		 * Associated Engineering model. An abstraction for type of analysis which will be prformed.
		 */
		EngngModel* m_engineeringModel;

		/**
		 * Domain connectivity table. Table is build upon request.
		 * Provides connectivity information of current domain.
		 */
		//std::unique_ptr< ConnectivityTable > connectivityTable;

		/// Output manager, allowing to filter the produced output.
		std::unique_ptr< OutputManager > outputManager;
		/// Domain number.
		int number;
		/// Domain serial (version) number. Used for domain version identification during Adaptive computations.
		int serialNumber;
		/// Number of spatial dimensions
		int nsd;
		bool axisymm;

		DomainType m_dType;
		std::string m_strdType;

		/**
		 * Map from an element's global number (label) to its place
		 * in the element array.
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

	public:
		/// Keeps track of next free dof ID (for special Lagrange multipliers, XFEM and such)
		int freeDofID;
	public:
		/**
		 * Constructor. Creates empty n-th domain belonging to given engineering model.
		 * @param n Domain number.
		 * @param serNum Serial number
		 * @param e Engineering model domain will belong to.
		 * @see giveSerialNumber
		 */
		Domain(int n, int serNum, EngngModel* e);

		Domain(const Domain& src) = delete;
		Domain& operator = (const Domain& src) = delete;

		/// Destructor.
		~Domain();

		/// Returns domain number.
		int giveNumber() { return this->number; }
		/// Returns domain number.
		void setNumber(int nn) { this->number = nn; }
		/// Returns domain serial (version) number.
		int giveSerialNumber() { return this->serialNumber; }

		EngngModel* giveEngngModel() { return m_engineeringModel; }

		void BuildDofManPlaceInArrayMap();
		void BuildElementPlaceInArrayMap();

		// management of the mesh components
		/**
		 * Service for accessing particular domain fe element.
		 * Generates error if no such element is defined.
		 * @param n Pointer to n-th element is returned.
		 */
		Element* giveElement(int n);
		std::vector< std::unique_ptr< Element > >& giveElements() { return this->m_elementList; }
		/**
		 * Service for accessing particular domain fe element.
		 * Generates error if no such element is defined.
		 * @param n Pointer to the element with id n
		 */
		Element* giveGlobalElement(int n);

		//void loadFile(DataReader& dr);

		int instanciateYourself(DataReader& dr);

		void resolveDomainDofsDefaults(const char* typeName);

	private:
	};
} //end namespace fem
#endif // domain_h
