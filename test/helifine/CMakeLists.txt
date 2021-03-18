cmake_minimum_required(VERSION 3.5.1)

project(out C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/TailorNew")

find_package(MPI REQUIRED)
find_package(Threads REQUIRED)
find_package(METIS REQUIRED)
find_package(MPFR REQUIRED)
find_package(Boost REQUIRED COMPONENTS mpi serialization program_options log thread system filesystem)

set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-deprecated")

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

include_directories(~/TailorNew/include ${Boost_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${METIS_INCLUDE_DIRS} /truba/home/osibliyev/amgcl/include)
link_directories(~/TailorNew/build)
add_executable(out main.cpp)
target_link_libraries(out tailor ${MPI_LIBRARIES} ${Boost_LIBRARIES} ${METIS_LIBRARIES} ${GMP_LIBRARIES} ${MPFR_LIBRARIES} Threads::Threads)
#target_link_libraries(out tailor ${MPI_LIBRARIES} ${Boost_LIBRARIES} ${METIS_LIBRARIES} Threads::Threads)

if(MPI_COMPILE_FLAGS)
    set_target_properties(out PROPERTIES
        COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(out PROPERTIES
        LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
