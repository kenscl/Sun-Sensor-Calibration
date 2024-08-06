cmake_minimum_required(VERSION 3.5)

project(sun-sensor-simulator CXX)
add_compile_options(-Wall)
set(CMAKE_CXX_STANDARD 11)

set(MATH
    math/matlib.h
    math/matlib.cpp
    math/matrix.h
    math/matrix.cpp
    math/quaternion.h
    math/quaternion.cpp
)

set(POLFIT
    polynomial_fit/polynomial_fit.h
    polynomial_fit/polynomial_fit.cpp
)

set(LUT
    lut/lut.h
    lut/lut.cpp
    lut/c-spline.h
    lut/c-spline.cpp
    lut/b-spline.h
    lut/b-spline.cpp
)
set(SINFIT
    sin_fit/sin_fit.h
    sin_fit/sin_fit.cpp
    sin_fit/extended_sin_fit.h
    sin_fit/extended_sin_fit.cpp
    sin_fit/extended_sin_fit_alt.h
    sin_fit/extended_sin_fit_alt.cpp
)

set(SOURCES 
    sun_sensor_data_parser.h
    sun_sensor.h
    sun_sensor.cpp
    main.cpp
)

set(TEST
    test/test.cpp
    test/test.h
)

set(ALIGNMENT
    alignment/main.cpp
    alignment/sat_model.h
    alignment/sat_model.cpp
    alignment/angle_corrector.h
    alignment/angle_corrector.cpp
    alignment/ekf.h
    alignment/ekf.cpp
    alignment/vasco.h
    alignment/vasco.cpp
    alignment/vasco.cpp
    alignment/data_reader.h
    alignment/data_reader.cpp
)

if(NOT DEFINED ENV{GITHUB_ACTIONS})
    set(MATPLOTLIB_CPP_DIR "matplotlib-cpp/")

    find_package(Python3 COMPONENTS Development NumPy REQUIRED)

    include_directories(${MATPLOTLIB_CPP_DIR})
    include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})

    add_executable(${PROJECT_NAME} ${MATH} ${SOURCES} ${POLFIT} ${SINFIT} ${LUT} ${TEST})

    target_include_directories(${PROJECT_NAME} PRIVATE ${MATPLOTLIB_CPP_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${Python3_LIBRARIES})

    add_executable(alignment ${ALIGNMENT} ${MATH} ${LUT} ${POLFIT} ${})
    target_include_directories(alignment PRIVATE ${MATPLOTLIB_CPP_DIR})
    target_link_libraries(alignment  PRIVATE ${Python3_LIBRARIES})

endif()

