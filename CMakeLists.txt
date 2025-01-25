cmake_minimum_required(VERSION 3.6)

project(multiPivot)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native -mavx -O3")
SET(EXECUTABLE_OUTPUT_PATH "bin/")

# set(exName1 "run") #maximal k plex
# set(exName2 "d2kRun")
# set(exName3 "saveAsBin")
# set(exName3_2 "saveAsTxt")
set(exName4 "pc")
set(exName5 "plv2")
# set(exName6 "msdc")
set(exName7 "pd")
# set(exName8 "pivoterRun")
# set(exName9 "kpivoterRun")
# set(exName10 "kcc") #kCListClique
set(exName11 "sdcc")
# set(exName12 "pikcc")
# set(exName13 "kpp")

# set(exName14 "cLocal")
set(exName15 "sELocal")
set(exName16 "pELocal")
set(exName17 "sVLocal")
set(exName18 "pVLocal")
#local counting per edge


add_subdirectory(src/graph)
add_subdirectory(src/plex)
add_subdirectory(src/sdc)
# add_subdirectory(src/clique)
add_subdirectory(src/local)
add_subdirectory(src/tools)

# add_executable(${exName1} src/run/run.cpp)
# add_executable(${exName2} src/run/d2kRun.cpp)
# add_executable(${exName3} src/run/saveAsBin.cpp)
# add_executable(${exName3_2} src/run/saveAsTxt.cpp)
add_executable(${exName4} src/run/runPlexCounting.cpp)
add_executable(${exName5} src/run/runKPList.cpp)
# add_executable(${exName6} src/run/msdcRunner.cpp)
add_executable(${exName7} src/run/kDListRunner.cpp)
# add_executable(${exName8} src/run/pivoterRunner.cpp)
# add_executable(${exName9} src/run/kpivoterRunner.cpp)
# add_executable(${exName10} src/run/kclistCliqueRunner.cpp)
add_executable(${exName11} src/run/runSdcCounting.cpp)
# add_executable(${exName12} src/run/kccPivotRunner.cpp)
# add_executable(${exName13} src/run/kppRunner.cpp)
# add_executable(${exName14} src/run/runCliqueLocal.cpp)
add_executable(${exName15} src/run/runSdcELocal.cpp)
add_executable(${exName16} src/run/runPlexELocal.cpp)
add_executable(${exName17} src/run/runSdcVLocal.cpp)
add_executable(${exName18} src/run/runPlexVLocal.cpp)



# target_link_libraries(${exName1} getArgs)
# target_link_libraries(${exName1} Graph)
# target_link_libraries(${exName1} plexEnumerator)

# target_link_libraries(${exName2} getArgs)
# target_link_libraries(${exName2} Graph)
# target_link_libraries(${exName2} d2k)

# target_link_libraries(${exName3} getArgs)
# target_link_libraries(${exName3} Graph)
# target_link_libraries(${exName3_2} getArgs)
# target_link_libraries(${exName3_2} Graph)

target_link_libraries(${exName4} getArgs)
target_link_libraries(${exName4} Graph)
target_link_libraries(${exName4} plexCounting)

target_link_libraries(${exName5} getArgs)
target_link_libraries(${exName5} Graph)
target_link_libraries(${exName5} kPList)

# target_link_libraries(${exName6} getArgs)
# target_link_libraries(${exName6} Graph)
# target_link_libraries(${exName6} maximalSDC)

target_link_libraries(${exName7} getArgs)
target_link_libraries(${exName7} Graph)
target_link_libraries(${exName7} kDList)

# target_link_libraries(${exName8} getArgs)
# target_link_libraries(${exName8} Graph)
# target_link_libraries(${exName8} pivoter)

# target_link_libraries(${exName9} getArgs)
# target_link_libraries(${exName9} Graph)
# target_link_libraries(${exName9} kpivoter)

# target_link_libraries(${exName10} getArgs)
# target_link_libraries(${exName10} Graph)
# target_link_libraries(${exName10} kCListClique)

target_link_libraries(${exName11} getArgs)
target_link_libraries(${exName11} Graph)
target_link_libraries(${exName11} sdcCountingDP)

# target_link_libraries(${exName12} getArgs)
# target_link_libraries(${exName12} Graph)
# target_link_libraries(${exName12} kccPivot)

# target_link_libraries(${exName13} getArgs)
# target_link_libraries(${exName13} Graph)
# target_link_libraries(${exName13} kccPlusPivoter)

# target_link_libraries(${exName14} getArgs)
# target_link_libraries(${exName14} Graph)
# target_link_libraries(${exName14} cliqueLocal)

target_link_libraries(${exName15} getArgs)
target_link_libraries(${exName15} Graph)
target_link_libraries(${exName15} sdcELocal)

target_link_libraries(${exName16} getArgs)
target_link_libraries(${exName16} Graph)
target_link_libraries(${exName16} plexELocal)

target_link_libraries(${exName17} getArgs)
target_link_libraries(${exName17} Graph)
target_link_libraries(${exName17} sdcVLocal)

target_link_libraries(${exName18} getArgs)
target_link_libraries(${exName18} Graph)
target_link_libraries(${exName18} plexVLocal)

