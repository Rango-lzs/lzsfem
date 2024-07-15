
#ifndef femcmpnn_h
#define femcmpnn_h

#include "interfacetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "fem_export.h"

#include <string>

namespace fem
{
	class DataStream;
	class Domain;
	class Interface;
	class TimeStep;
	class InputRecord;
	class DynamicInputRecord;
	class oofegGraphicContext;
	class EntityRenumberingFunctor;
	class FloatArray;
	class IntArray;
	class FloatMatrix;

	/**
	 * The top abstract class of all classes constituting the finite element mesh.
	 * Defines the attributes and methods common to all components of mesh:
	 * elements, nodes, time steps, materials, loads and load-time functions.
	 * This class defines the two attributes common to all component classes ;
	 * 'number' is primarily used for reading data in the data file. 'domain' is
	 * used for communicating with other components (e.g., for an element to obtain its material),
	 * for accessing the linear system and the data file.
	 * @see error handles error reporting.
	 * @see checkConsistency to ensure, whether internal data structures are consistent.
	 */

	/*
	* 提供一些通用基础能力
	*/
	class FEM_EXPORT FEMComponent
	{
	protected:
		/// Component number
		int number;
		/// Link to domain object, useful for communicating with other FEM components
		Domain* domain;

	public:
		/**
		 * Regular constructor, creates component with given number and belonging to given domain.
		 * @param n Component number in particular domain. For instance, can represent
		 * node number in particular domain.
		 * @param d Domain to which component belongs to.
		 */
		FEMComponent(int n, Domain* d) : number(n), domain(d) { }
		/// Virtual destructor.
		virtual ~FEMComponent() = default;
	};
} // end namespace fem
#endif // femcmpnn_h
