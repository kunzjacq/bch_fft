cmake_minimum_required(VERSION 3.1)
project(fft_bch)

set(OPT_FLAGS "-O3 -march=native -mtune=native")
set(DEBUG_FLAGS "-g -D_DEBUG")
set(WARN_FLAGS "-Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${OPT_FLAGS} -DNDEBUG")
set(CMAKE_C_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 ${WARN_FLAGS}")

include_directories("${PROJECT_SOURCE_DIR}/include")

#gf2x lib and include should be in in 'external' subdir.
include_directories("${PROJECT_SOURCE_DIR}/external/include")
link_directories("${PROJECT_SOURCE_DIR}/external/lib")

aux_source_directory(libbase BASE_LIB_SRC)
aux_source_directory(libbch BCH_LIB_SRC)
aux_source_directory(libfft FFT_LIB_SRC)

aux_source_directory(bch BCH_DEMO_SRC)
aux_source_directory(fft FFT_DEMO_SRC)

add_library(base ${BASE_LIB_SRC})
add_library(bch ${BCH_LIB_SRC})
add_library(fft ${FFT_LIB_SRC})

add_executable(fftdemo ${FFT_DEMO_SRC})
target_link_libraries(fftdemo fft base gf2x)

add_executable(bchdemo ${BCH_DEMO_SRC})
target_link_libraries(bchdemo bch fft base gf2x)

