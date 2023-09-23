/*****************************************************************//**
 * \file   domain.h
 * \brief  Description of the fem model
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef domain_h
#define domain_h

#include <unordered_map>

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

namespace fem 
{
	class Element;
	class DofManager;
	class CrossSection;
	class Material;
	class BoundaryCondition;
	class InitialCondition;
	class Function;
	class FemObjectSet;

	class EngngModel;
	class ConnectivityTable;
	class OutputManager;
	class TopologyDescription;

	class DataReader;

	class Domain
	{
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
		std::vector< std::unique_ptr< BoundaryCondition > > m_bcList;
		/// Initial condition list.
		std::vector< std::unique_ptr< InitialCondition > > m_icList;
		/// Load time function list.
		std::vector< std::unique_ptr< Function > > m_functionList;
		/// Set list.
		std::vector< std::unique_ptr< FemObjectSet > > m_setList;
				
		/**
		 * Associated Engineering model. An abstraction for type of analysis which will be prformed.
		 */
		EngngModel* m_engineeringModel;

		/**
		 * Domain connectivity table. Table is build upon request.
		 * Provides connectivity information of current domain.
		 */
		std::unique_ptr< ConnectivityTable > connectivityTable;

		/// Output manager, allowing to filter the produced output.
		std::unique_ptr< OutputManager > outputManager;
		/// Domain number.
		int number;
		/// Domain serial (version) number. Used for domain version identification during Adaptive computations.
		int serialNumber;
		/// Number of spatial dimensions
		int nsd;
		bool axisymm;

		std::string mDomainType;

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

		/// Topology description
		std::unique_ptr< TopologyDescription > topology;

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

		void loadFile(DataReader& dr);

	private:
	};
} //end namespace fem
#endif // domain_h
