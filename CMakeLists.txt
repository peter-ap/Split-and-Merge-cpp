cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(CornerSearch)


find_package(Eigen3 3.3 REQUIRED NO_MODULE)

FIND_PATH(BOOST_DIR "boost")

find_package(PythonLibs REQUIRED)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIRS}
    ${BOOST_DIR}
    usr/include/python3.6
    ${PYTHON_INCLUDE_DIRS}
  )

add_library(SaMlib SHARED
  src/CornerSearch.cpp
  src/dataReader.cpp
  src/visualizing.cpp
)
target_link_libraries(SaMlib ${PYTHON_LIBRARIES})

add_executable (CornerSearch main.cpp)
target_link_libraries(CornerSearch ${PYTHON_LIBRARIES} SaMlib)

# Testing, you could use testing frameworks such as CTest in the future
# then you can run "make test" for example
add_executable (testDataReader test/testDataReader.cpp)
target_link_libraries(testDataReader SaMlib)
