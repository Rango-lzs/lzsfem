/*****************************************************************//**
 * \file   primaryfield.h
 * \brief  
 * 
 * \author Leizs
 * \date   October 2023
 *********************************************************************/

#ifndef primaryfield_h
#define primaryfield_h

#include "field.h"
#include "floatarray.h"
#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "timestep.h"

#include <vector>

namespace fem
{
	class PrimaryField;
	class Dof;
	class BoundaryCondition;
	class InitialCondition;
	class UnknownNumberingScheme;


	/**
	 * Abstract class representing field of primary variables (those, which are unknown and are typically associated to nodes).
	 * In the current design the primary field is understood as simple database, that allows to keep track
	 * of the history of a solution vector representing primary field. The history is kept as sequence of
	 * solution vectors. The history depth kept can be selected. PrimaryField class basically provides access to
	 * time-dependent vectors of the field of unknowns. It adds the possibility to
	 * further interpolate the field values using element interpolation functions.
	 * The prescribed values of the field are not maintained, since they can be obtained directly from
	 * corresponding DOFs of associated domain.
	 *
	 * As the PrimaryField stores the state directly in solution vectors that are usually directly
	 * updated by EngngModel, it may contain a mix of different fields (this is especially true for
	 * strongly coupled problems). Then masked primary field can be used to select only certain DOFs
	 * (based on DofID) from its master PrimaryField.
	 *
	 * @note This primary field will always assume default numbering schemes.
	 */
	class OOFEM_EXPORT PrimaryField : public Field
	{
	protected:
		int actualStepNumber;
		int actualStepIndx;
		int nHistVectors;
		std::vector< FloatArray >solutionVectors;
		std::vector< FloatArray >prescribedVectors;
		std::vector< TimeStep >solStepList;
		EngngModel* emodel;
		int domainIndx;

	public:
		/**
		 * Constructor. Creates a field of given type associated to given domain.
		 * Not using pointer to domain, because this will prevent the use of PrimaryField as an
		 * EngngModel attribute. This is because the domain does not exists when
		 * PrimaryField is created (this is when EngngModel is created).
		 * @param a Engineering model which field belongs to.
		 * @param idomain Index of domain for field.
		 * @param ft Type of stored field.
		 * @param nHist Number of old time steps to store.
		 */
		PrimaryField(EngngModel* a, int idomain, FieldType ft, int nHist);
		virtual ~PrimaryField();

		/**
		 * Copy unknowns from previous solution or DOF's dictionary to the solution vector
		 * @param mode what the unknown describes (increment, total value etc.).
		 * @param tStep Time of interest.
		 * @param answer Resulting vector.
		 */
		virtual void initialize(ValueModeType mode, TimeStep* tStep, FloatArray& answer, const UnknownNumberingScheme& s);

		// These functions are hardcoded to assume the default numbering scheme (as is the rest of primaryfield)
		void storeDofManager(TimeStep* tStep, DofManager& dman);
		void storeInDofDictionaries(TimeStep* tStep);
		void readDofManager(TimeStep* tStep, DofManager& dman);
		void readFromDofDictionaries(TimeStep* tStep);

		/**
		 * Applies the default initial values values for all DOFs (0) in given domain.
		 * @param domain Domain number
		 */
		virtual void applyDefaultInitialCondition();
		/**
		 * Applies initial condition to all DOFs.
		 * @param ic Initial condition for DOFs
		 */
		void applyInitialCondition(InitialCondition& ic);
		/**
		 * Applies all boundary conditions to all prescribed DOFs.
		 * @param tStep Current time step.
		 */
		virtual void applyBoundaryCondition(TimeStep* tStep);
		/**
		 * Applies the boundary condition to all prescribed DOFs in given domain.
		 * @param bc Boundary condition.
		 * @param domain tStep Time step for when bc applies.
		 */
		void applyBoundaryCondition(BoundaryCondition& bc, TimeStep* tStep);

		/**
		 * @param dof Pointer to DOF.
		 * @param mode What the unknown describes (increment, total value etc.).
		 * @param tStep Time step of interest.
		 * @return Value of interest at given DOF.
		 */
		virtual double giveUnknownValue(Dof* dof, ValueModeType mode, TimeStep* tStep);

		int evaluateAt(FloatArray& answer, const FloatArray& coords, ValueModeType mode, TimeStep* tStep) override;
		int evaluateAt(FloatArray& answer, DofManager* dman, ValueModeType mode, TimeStep* tStep) override;

		/**
		 * Evaluates the field at given DOF manager, allows to select specific
		 * dofs using mask
		 * @param answer Evaluated field at dman.
		 * @param dman DOF manager of interest.
		 * @param mode Mode of evaluated unknowns.
		 * @param tStep Time step of interest.
		 * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
		 * @return Error code (0=ok, 1=point not found in domain).
		 */
		virtual int __evaluateAt(FloatArray& answer, DofManager* dman,
			ValueModeType mode, TimeStep* tStep, IntArray* dofId);
		/**
		 * Evaluates the field at given point, allows to select specific
		 * dofs using mask.
		 * @param answer Evaluated field at coords.
		 * @param coords Coordinates of the point of interest.
		 * @param mode Mode of evaluated unknowns.
		 * @param tStep Time step of interest.
		 * @param dofId Dof mask, id set to NULL, all Dofs evaluated.
		 * @return Error code (0=ok, 1=point not found in domain)
		 */
		virtual int __evaluateAt(FloatArray& answer, const FloatArray& coords,
			ValueModeType mode, TimeStep* tStep, IntArray* dofId);
		/**
		 * @param tStep Time step to take solution for.
		 * @return Solution vector for requested time step.
		 */
		virtual FloatArray* giveSolutionVector(TimeStep* tStep);

		/**
		 * Project vectorToStore back to DOF's dictionary
		 * @param vectorToStore Vector with the size of number of equations.
		 * @param mode Mode of the unknown (increment, total value etc.)
		 * @param tStep Time step unknowns belong to.
		 */
		virtual void update(ValueModeType mode, TimeStep* tStep, const FloatArray& vectorToStore, const UnknownNumberingScheme& s);

		/**
		 * Brings up a new solution vector for given time step.
		 * @param tStep Time step for new solution vector.
		 */
		virtual void advanceSolution(TimeStep* tStep);

		void saveContext(DataStream& stream) override;
		void restoreContext(DataStream& stream) override;

		const char* giveClassName() const override { return "PrimaryField"; }

		int giveActualStepNumber() { return actualStepNumber; }
	protected:
		int resolveIndx(TimeStep* tStep, int shift);
		FloatArray* giveSolutionVector(int);
		FloatArray* givePrescribedVector(int);
	};
} // end namespace fem
#endif // primaryfield_h
