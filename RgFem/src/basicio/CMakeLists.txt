set(HEADERS
    DataRecord.h
    DataStore.h
    DomainDataRecord.h
    DumpFile.h
    DumpMemStream.h
    DumpStream.h
    ElementDataRecord.h
    FEModelDataRecord.h
    FEPlotData.h
    FEPlotDataStore.h
    FaceDataRecord.h
    NLConstraintDataRecord.h
    NodeDataRecord.h
    ObjectDataRecord.h
    SurfaceDataRecord.h
)

set(SOURCES
    DataRecord.cpp
    DataStore.cpp
    DomainDataRecord.cpp
    DumpFile.cpp
    DumpMemStream.cpp
    DumpStream.cpp
    ElementDataRecord.cpp
    FEModelDataRecord.cpp
    FEPlotData.cpp
    FEPlotDataStore.cpp
    FaceDataRecord.cpp
    NLConstraintDataRecord.cpp
    NodeDataRecord.cpp
    ObjectDataRecord.cpp
    SurfaceDataRecord.cpp
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
