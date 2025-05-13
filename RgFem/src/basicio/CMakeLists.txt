set(HEADERS
    basicio/DataRecord.h
    basicio/DataStore.h
    basicio/DomainDataRecord.h
    basicio/DumpFile.h
    basicio/DumpMemStream.h
    basicio/DumpStream.h
    basicio/ElementDataRecord.h
    basicio/FEModelDataRecord.h
    basicio/FEPlotData.h
    basicio/FEPlotDataStore.h
    basicio/FaceDataRecord.h
    basicio/NLConstraintDataRecord.h
    basicio/NodeDataRecord.h
    basicio/ObjectDataRecord.h
    basicio/SurfaceDataRecord.h
)

set(SOURCES
    basicio/DataRecord.cpp
    basicio/DataStore.cpp
    basicio/DomainDataRecord.cpp
    basicio/DumpFile.cpp
    basicio/DumpMemStream.cpp
    basicio/DumpStream.cpp
    basicio/ElementDataRecord.cpp
    basicio/FEModelDataRecord.cpp
    basicio/FEPlotData.cpp
    basicio/FEPlotDataStore.cpp
    basicio/FaceDataRecord.cpp
    basicio/NLConstraintDataRecord.cpp
    basicio/NodeDataRecord.cpp
    basicio/ObjectDataRecord.cpp
    basicio/SurfaceDataRecord.cpp
)

set(BASICIO_SOURCES_
    ${HEADERS}
    ${SOURCES}
)

set(PROJECT_SOURCES 
    ${PROJECT_SOURCES} 
    ${BASICIO_SOURCES_} 
    PARENT_SCOPE
)

set(BASICIO_SOURCES ${BASICIO_SOURCES_} PARENT_SCOPE)
