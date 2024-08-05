include_directories(../graph ../tools)

# add_library(pivoter pivoter.h pivoter.cpp)
# add_library(pivoter pivoter.h pivoterKCListDouble.cpp)
# add_library(pivoter pivoter.h pivoterSG.cpp)
# add_library(pivoter pivoter.h pivoterSmallPivot.cpp)
add_library(pivoter pivoter.h pivoterSG.cpp)
target_link_libraries(pivoter Graph linearSet types)

add_library(kpivoter kpivoter.h kpivoterSGClique.cpp)
target_link_libraries(kpivoter Graph linearSet types)

add_library(kCListClique kCListClique.h kCListClique.cpp)
target_link_libraries(kCListClique Graph linearSet types)

add_library(kccPivot kccPivot.h kccPivot.cpp MCEGO.cpp)
target_link_libraries(kccPivot Graph linearSet types)

add_library(kccPlusPivoter kccPlusPivoter.h kccPlusPivoter.cpp)
target_link_libraries(kccPlusPivoter Graph linearSet types)
