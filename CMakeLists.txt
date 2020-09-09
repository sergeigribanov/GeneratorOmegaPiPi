cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(GeneratorOmegaPiPi)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(GSL REQUIRED)
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})
find_package(MEUtils REQUIRED)
find_package(MElCalc REQUIRED)
find_package(RandomStar REQUIRED)
find_package(EventOmegaPiPi REQUIRED)
include_directories(${MEUtils_INCLUDE_DIRS})
include_directories(${MElCalc_INCLUDE_DIRS})
include_directories(${RandomStar_INCLUDE_DIRS})
include_directories(${EventOmegaPiPi_INCLUDE_DIRS})

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB SOURCES_NOMAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(REMOVE_ITEM SOURCES_NOMAIN main.cpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/gen-omegapipi-setup.sh.in gen-omegapipi-setup.sh)
add_library(GeneratorOmegaPiPi SHARED ${SOURCES_NOMAIN})
target_link_libraries(GeneratorOmegaPiPi
  PRIVATE
  GSL::gsl
  GSL::gslcblas
  ${ROOT_LIBRARIES}
  ${MEUtils_LIBRARIES}
  ${MElCalc_LIBRARIES}
  ${RandomStar_LIBRARIES}
  ${EventOmegaPiPi_LIBRARIES})
target_include_directories(GeneratorOmegaPiPi
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
  ${GSL_INCLUDE_DIRS})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(GeneratorOmegaPiPi
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
add_executable(gen-omegapipi-test ${SOURCES_NOMAIN})
target_link_libraries(gen-omegapipi-test
  GSL::gsl
  GSL::gslcblas
  ${ROOT_LIBRARIES}
  ${MEUtils_LIBRARIES}
  ${MElCalc_LIBRARIES}
  ${RandomStar_LIBRARIES}
  ${EventOmegaPiPi_LIBRARIES})
target_include_directories(gen-omegapipi-test
  PRIVATE
  ${GSL_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/include)
install(TARGETS GeneratorOmegaPiPi
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)
install(TARGETS gen-omegapipi-test DESTINATION bin)
install(FILES ${PROJECT_BINARY_DIR}/gen-omegapipi-setup.sh DESTINATION env)
