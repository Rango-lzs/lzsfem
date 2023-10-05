
#ifndef nodalload_h
#define nodalload_h

#include "load.h"
#include "bcgeomtype.h"
#include "valuemodetype.h"

namespace fem
{
	///@name Input fields for nodal loads
	//@{
#define _IFT_NodalLoad_Name "nodalload"
#define _IFT_NodalLoad_cstype "cstype"
//@}

	class TimeStep;

	/**
	 * Class implementing a concentrated load (force, moment,...) that acts
	 * directly on a dof manager (node or element side, if it has associated DOFs).
	 * This load could not be applied on an element.
	 * A nodal load is usually attribute of one or more nodes or element sides.
	 *
	 * The component array, which size should be same as number of DOFs in particular
	 * node/side is read from input.
	 *
	 * The attribute componentArray contains, for example for the case of a
	 * plane beam structure, 2 forces and 1 moment on a right place. (6 dof per
	 * node is assumed)
	 *
	 * @note{Load is not restricted to structural problems. For example, in thermal
	 * analysis, a nodal load would be a concentrated heat source.}
	 */
	class FEM_EXPORT NodalLoad : public Load
	{
	protected:
		/**
		 * Load coordinate system.
		 * It is actually used only when local coordinate system in node is defined and load is specified in global
		 * coordinate system
		 */
		CoordSystType coordSystemType;

	public:
		/**
		 * Constructor. Creates nodal load object with given number, belonging to given domain.
		 * @param n Load  number.
		 * @param d Domain to which new object will belongs.
		 */
		NodalLoad(int n, Domain* d) : Load(n, d) { }

		const char* giveInputRecordName() const override { return _IFT_NodalLoad_Name; }
		void computeValueAt(FloatArray& answer, TimeStep* tStep, const FloatArray& coords, ValueModeType mode) override
		{
			computeComponentArrayAt(answer, tStep, mode);
		}
		CoordSystType giveCoordSystMode() override { return coordSystemType; }

		void initializeFrom(InputRecord& ir) override;
		void giveInputRecord(DynamicInputRecord& input) override;
		const char* giveClassName() const override { return "NodalLoad"; }
		bcGeomType giveBCGeoType() const override { return NodalLoadBGT; }

		void saveContext(DataStream& stream, ContextMode mode) override;
		void restoreContext(DataStream& stream, ContextMode mode) override;
	};
} // end namespace oofem
#endif // nodalload_h
