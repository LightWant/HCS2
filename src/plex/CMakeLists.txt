include_directories(../graph ../tools)

add_library(plexEnumerator plexEnumerator.h plexEnumerator4.cpp candidate.hpp)
target_link_libraries(plexEnumerator Graph linearSet types)

add_library(d2k d2k.h d2kng.cpp)
target_link_libraries(d2k Graph linearSet types)

add_library(plexCounting plexCounting.h plexCountingV10.cpp)
target_link_libraries(plexCounting Graph linearSet types)

add_library(kPList kPList.h kPList.cpp)
target_link_libraries(kPList Graph types linearSet)