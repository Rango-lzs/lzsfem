/*****************************************************************//**
 * \file   text_data_reader.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef FEMTXTDataReader_h
#define FEMTXTDataReader_h

#include "datareader.h"
#include "text_input_record.h"

#include <fstream>

namespace fem {
	/**
	 * Class representing the implementation of plain text date reader.
	 * It reads a sequence of input records from data file
	 * and creates the corresponding input records.
	 * There is no check for record type requested, it is assumed that records are
	 * written in correct order, which determined by the coded sequence of
	 * component initialization and described in input manual.
	 */
	class FEM_EXPORT FEMTXTDataReader : public DataReader
	{
	protected:
		std::string dataSourceName;
		std::list< FEMTXTInputRecord > recordList;

		/// Keeps track of the current position in the list
		std::list< FEMTXTInputRecord > ::iterator it;

	public:
		/// Constructor.
		FEMTXTDataReader(std::string inputfilename);
		FEMTXTDataReader(const FEMTXTDataReader& x);
		virtual ~FEMTXTDataReader();

		InputRecord& giveInputRecord(InputRecordType, int recordId) override;
		bool peakNext(const std::string& keyword) override;
		void finish() override;
		std::string giveReferenceName() const override { return dataSourceName; }

	protected:
		/**
		 * Reads one line from inputStream
		 * Parts within quotations have case preserved.
		 */
		bool giveLineFromInput(std::ifstream& stream, int& lineNum, std::string& line);
		/// Reads one line from stream.
		bool giveRawLineFromInput(std::ifstream& stream, int& lineNum, std::string& line);
	};
} // end namespace fem
#endif // FEMTXTDataReader_h
