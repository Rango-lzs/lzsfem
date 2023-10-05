#ifndef field_h
#define field_h

#include "valuemodetype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "enumitem.h"
#include <string>
#include <memory>

namespace fem
{
	///@todo FieldType and UnknownType basically determine the same thing. Should be possible to stick to one. Combinations of fields should be possible with logical bitfields.
#define FieldType_DEF \
    ENUM_ITEM_WITH_VALUE(FT_Unknown, 0) \
    ENUM_ITEM_WITH_VALUE(FT_Velocity, 1) \
    ENUM_ITEM_WITH_VALUE(FT_Displacements, 2) \
    ENUM_ITEM_WITH_VALUE(FT_VelocityPressure, 3) \
    ENUM_ITEM_WITH_VALUE(FT_Pressure, 4) \
    ENUM_ITEM_WITH_VALUE(FT_Temperature, 5) \
    ENUM_ITEM_WITH_VALUE(FT_HumidityConcentration, 6) \
    ENUM_ITEM_WITH_VALUE(FT_TransportProblemUnknowns, 7) \
    ENUM_ITEM_WITH_VALUE(FT_TemperatureAmbient, 8)

/// Physical type of field.
	enum FieldType {
		FieldType_DEF
	};
#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

	class TimeStep;
	class FloatArray;
	class DofManager;
	class DataStream;

	class Field;
	typedef std::shared_ptr<Field> FieldPtr;

	/**
	 * Abstract class representing field. Field represent the spatial distribution of certain variable.
	 * Field is able to evaluate its value at any point of interest. The field is usually associated to
	 * specific domain.
	 */
	class FEM_EXPORT Field
	{
	protected:
		FieldType type;

	public:
		/**
		 * Constructor. Creates a field of given type associated to given domain.
		 */
		Field(FieldType b) : type(b) { }
		virtual ~Field() { }
		/**
		 * Evaluates the field at given point.
		 * @param coords Coordinates of the point of interest
		 * @param answer Field evaluated at coordinate.
		 * @param tStep Time step to evaluate for.
		 * @param mode Mode of value (total, velocity,...).
		 * @return Zero if ok, otherwise nonzero.
		 */
		virtual int evaluateAt(FloatArray& answer, const FloatArray& coords,
			ValueModeType mode, TimeStep* tStep) = 0;

		/**
		 * Evaluates the field at given DofManager. This potentially can be resolved quickly, as
		 * receiver data may be described using values at dofManagers. Here an additional issue
		 * exists: one needs to make sure, that passed dman is from the same domain, so that its
		 * number can be used to perform suggested quick evaluation.
		 *
		 * If this is not the case (the field is described differently),
		 * the response can be evaluated using dofman coordinates in a standard way.
		 * @param[out] answer Evaluated field for dman.
		 * @param dman Reference to dofManager.
		 * @param mode Mode of value (total, velocity,...).
		 * @param tStep Time step to evaluate for.
		 * @return Zero if ok, nonzero Error code (0-ok, 1-failed)
		 */
		virtual int evaluateAt(FloatArray& answer, DofManager* dman,
			ValueModeType mode, TimeStep* tStep) = 0;

		/// Returns the type of receiver
		FieldType giveType() { return type; }

		/// Sets the type of receiver
		void setType(FieldType b) { type = b; }

		/**
		 * Stores receiver state to output stream.
		 * @param stream Output stream.
		 * @exception Throws an ContextIOERR exception if error encountered.
		 */
		virtual void saveContext(DataStream& stream) = 0;
		/**
		 * Restores the receiver state previously written in stream.
		 * @param stream Input stream.
		 * @exception Throws an ContextIOERR exception if error encountered.
		 */
		virtual void restoreContext(DataStream& stream) = 0;


		/// Returns string for prepending output (used by error reporting macros).
		std::string errorInfo(const char* func) const;

		/// @return Class name of the receiver.
		virtual const char* giveClassName() const = 0;
	};
} // end namespace fem
#endif // field_h
