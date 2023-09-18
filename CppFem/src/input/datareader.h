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

namespace rango {
/**
 * Class representing the abstraction for input data source.
 * Its role is to provide input records for particular components.
 * The input records are identified by record type and component number.
 * The order of input records is in fact determined by the coded sequence of
 * component initialization. The input record identification facilitates the
 * implementation of database readers with direct or random access.
 */
class FEM_EXPORT DataReader
{
protected:
    /// Output file name (first line in OOFEM input files).
    std :: string outputFileName;
    /// Description line (second line in OOFEM input files).
    std :: string description;

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
    virtual InputRecord &giveInputRecord(InputRecordType irType, int recordId) = 0;

    /**
     * Peak in advance into the record list.
     * @return True if next keyword is a set.
     */
    virtual bool peakNext(const std :: string &keyword) { return false; }

    /**
     * Allows to detach all data connections.
     */
    virtual void finish() = 0;

    /// Gives the reference file name (e.g. file name)
    virtual std :: string giveReferenceName() const = 0;
    /// Gives the output file name
    std :: string giveOutputFileName() { return this->outputFileName; }
    /// Gives the problem description
    std :: string giveDescription() { return this->description; }
};
} // end namespace fem
#endif // datareader_h
