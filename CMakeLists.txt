cmake_minimum_required(VERSION 3.3)
project(clique-partitioning)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

include_directories(include)
file(GLOB headers
    include/*.hxx
)
enable_testing()

find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})

# unit tests
add_executable(test-clique-partitioning-bnc src/unit-tests/clique-partitioning-bnc.cxx ${headers})
target_link_libraries(test-clique-partitioning-bnc ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
add_test(test-clique-partitioning-bnc test-clique-partitioning-bnc)

add_executable(test-property-maps src/unit-tests/property-maps.cxx ${headers})
add_test(test-property-maps test-property-maps)

add_executable(test-load src/unit-tests/load.cxx ${headers})
add_test(test-load test-load)

add_executable(test-callback src/unit-tests/callback.cxx ${headers})
add_test(test-callback test-callback)

# inequalities
add_executable(test-inequality src/unit-tests/inequalities/inequality.cxx ${headers})
add_test(test-inequality test-inequality)

add_executable(test-odd-wheel src/unit-tests/inequalities/odd-wheel.cxx ${headers})
add_test(test-odd-wheel test-odd-wheel)

add_executable(test-two-partition src/unit-tests/inequalities/two-partition.cxx ${headers})
add_test(test-two-partition test-two-partition)

add_executable(test-hypermetric src/unit-tests/inequalities/hypermetric.cxx ${headers})
add_test(test-hypermetric test-hypermetric)

add_executable(test-separator src/unit-tests/inequalities/separator.cxx ${headers})
add_test(test-separator test-separator)

add_executable(test-complete src/unit-tests/inequalities/complete.cxx ${headers})
add_test(test-complete test-complete)

add_executable(test-chvatal-gomory src/unit-tests/inequalities/chvatal-gomory.cxx ${headers})
target_link_libraries(test-chvatal-gomory ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
add_test(test-chvatal-gomory test-chvatal-gomory)

add_executable(test-half-chorded-odd-cycle src/unit-tests/inequalities/half-chorded-odd-cycle.cxx ${headers})
add_test(test-half-chorded-odd-cycle test-half-chorded-odd-cycle)

add_executable(test-clique-web src/unit-tests/inequalities/clique-web.cxx ${headers})
target_link_libraries(test-clique-web ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
add_test(test-clique-web test-clique-web)

# main
add_executable(main src/main.cxx ${headers})
target_link_libraries(main ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
target_compile_features(main PUBLIC cxx_std_17)
