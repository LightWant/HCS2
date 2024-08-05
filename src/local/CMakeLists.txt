include_directories(../graph ../tools)

add_library(cliqueLocal cliqueLocal.h cliqueLocal.cpp)
target_link_libraries(cliqueLocal Graph linearSet types)

add_library(sdcELocal sdcLocal.h sdcELocal.cpp)
target_link_libraries(sdcELocal Graph linearSet types)

add_library(plexELocal plexLocal.h plexELocal.cpp)
target_link_libraries(plexELocal Graph linearSet types)

add_library(sdcVLocal sdcLocal.h sdcVLocal.cpp)
target_link_libraries(sdcVLocal Graph linearSet types)

add_library(plexVLocal plexLocal.h plexVLocal.cpp)
target_link_libraries(plexVLocal Graph linearSet types)
