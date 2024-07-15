/*****************************************************************//**
 * \file   dynamic_input_record.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef dynamicinputrecord_h
#define dynamicinputrecord_h

#include "inputrecord.h"

#include <map>
#include <set>

namespace fem {
	class FEMComponent;
	class DynamicInputRecord;

	/// Helper function for creating a dynamic input record for a node
	FEM_EXPORT std::unique_ptr<DynamicInputRecord> CreateNodeIR(int i, InputFieldType nodeType, FloatArray coord);

	/// Helper function for creating elements (with optional cross-section number).
	FEM_EXPORT std::unique_ptr<DynamicInputRecord> CreateElementIR(int i, InputFieldType elementType, IntArray nodes, int cs = 0);

	/**
	 * Class representing the a dynamic Input Record.
	 * The input record is represented as a list of fields.
	 * This is intended for internal usage, where new input records and such are created dynamically.
	 * @author Mikael Öhman
	 */
	class FEM_EXPORT DynamicInputRecord : public InputRecord
	{
	protected:
		std::string recordKeyword;
		int recordNumber;

		std::set< std::string >emptyRecord; ///< Fields without values
		std::map< std::string, int >intRecord;
		std::map< std::string, double >doubleRecord;
		std::map< std::string, bool >boolRecord;
		std::map< std::string, std::string >stringRecord;
		std::map< std::string, FloatArray >floatArrayRecord;
		std::map< std::string, IntArray >intArrayRecord;
		std::map< std::string, FloatMatrix >matrixRecord;
		std::map< std::string, std::vector< std::string > >stringListRecord;
		std::map< std::string, Dictionary >dictionaryRecord;
		std::map< std::string, std::list< Range > >rangeRecord;
		std::map< std::string, ScalarFunction >scalarFunctionRecord;

	public:
		/// Creates an empty input record.
		DynamicInputRecord(std::string answer = "", int value = 0);
		/// Extracts input record from given component
		DynamicInputRecord(FEMComponent& femc); ///@todo Make FEMComponent const
		/// Copy constructor.
		DynamicInputRecord(const DynamicInputRecord&);
		/// Destructor.
		virtual ~DynamicInputRecord();
		/// Assignment operator.
		DynamicInputRecord& operator = (const DynamicInputRecord&);

		std::unique_ptr<InputRecord> clone() const override { return std::make_unique<DynamicInputRecord>(*this); }
		void finish(bool wrn = true) override;

		std::string giveRecordAsString() const override;

		void giveRecordKeywordField(std::string& answer, int& value) override;
		void giveRecordKeywordField(std::string& answer) override;
		void giveField(int& answer, InputFieldType id) override;
		void giveField(double& answer, InputFieldType id) override;
		void giveField(bool& answer, InputFieldType id) override;
		void giveField(std::string& answer, InputFieldType id) override;
		void giveField(FloatArray& answer, InputFieldType id) override;
		void giveField(IntArray& answer, InputFieldType id) override;
		void giveField(FloatMatrix& answer, InputFieldType id) override;
		void giveField(std::vector< std::string >& answer, InputFieldType id) override;
		void giveField(Dictionary& answer, InputFieldType id) override;
		void giveField(std::list< Range >& answer, InputFieldType id) override;
		void giveField(ScalarFunction& function, InputFieldType id) override;

		bool hasField(InputFieldType id) override;
		void printYourself() override;

		// Setters, unique for the dynamic input record
		void setRecordKeywordField(std::string keyword, int number);
		void setRecordKeywordNumber(int number);
		void setField(int item, InputFieldType id);
		void setField(double item, InputFieldType id);
		void setField(bool item, InputFieldType id);
		void setField(std::string item, InputFieldType id);
		void setField(FloatArray item, InputFieldType id);
		void setField(IntArray item, InputFieldType id);
		void setField(FloatMatrix item, InputFieldType id);
		void setField(std::vector< std::string > item, InputFieldType id);
		void setField(const Dictionary& item, InputFieldType id);
		void setField(const std::list< Range >& item, InputFieldType id);
		void setField(const ScalarFunction& function, InputFieldType id);
		/// Sets an empty field with given id.
		void setField(InputFieldType id);
		/// Removes given field from record.
		void unsetField(InputFieldType id);
	};
} // end namespace fem
#endif // dynamicinputrecord_h
