/*****************************************************************//**
 * \file   datareader.h
 * \brief
 *
 * \author Leizs
 * \date   September 2023
 *********************************************************************/

#ifndef datareader_h
#define datareader_h

#include "inputrecord.h"
#include <memory>

namespace fem {
	/**
	 * Abstract class representing the input data resource.
	 * 数据源可能是自定义格式的文件，其他求解器格式文件，动态输入文件
	 */
	class FEM_EXPORT DataReader
	{
	protected:
		/// Output file name (first line in OOFEM input files).
		std::string m_outputFileName;
		/// Description line (second line in OOFEM input files).
		std::string m_description;

		std::vector<std::unique_ptr<InputRecord>> m_records;

	public:
		/// Determines the type of input record.
		enum InputRecordType {
			IR_domainRec, IR_outManRec, IR_domainCompRec, IR_geometryRec, IR_gbpmRec,
			IR_emodelRec, IR_mstepRec, IR_expModuleRec, IR_dofmanRec, IR_elemRec,
			IR_crosssectRec, IR_matRec, IR_nlocBarRec, IR_bcRec, IR_icRec, IR_funcRec, IR_setRec,
			IR_xfemManRec, IR_enrichFuncRec, IR_geoRec, IR_enrichItemRec,
			IR_enrichFrontRec, IR_propagationLawRec, IR_crackNucleationRec, IR_fracManRec, IR_failCritRec,
			IR_contactManRec, IR_contactDefRec
		};

		DataReader() { }
		virtual ~DataReader() { }

		/**
		 * Returns input record corresponding to given InputRecordType value and its record_id.
		 * The returned InputRecord reference is valid only until the next call.
		 * @param irType Determines type of record to be returned.
		 * @param recordId Determines the record  number corresponding to component number.
		 */
		virtual InputRecord& giveInputRecord(InputRecordType irType, int recordId) = 0;

		//Gives the output file name
		const std::string& giveOutputFileName() { return this->m_outputFileName; }

		//Gives the problem description
		const std::string& giveDescription() { return this->m_description; }
	};
} // end namespace fem
#endif // datareader_h
