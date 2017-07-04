project(ExtendedKF)

cmake_minimum_required (VERSION 2.8)

add_definitions(-std=c++11)
add_compile_options(-Werror -Wall)

set(sources src/tools.cpp src/FusionEKF.cpp src/kalman_filter.cpp src/tools.h src/FusionEKF.h src/kalman_filter.h)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 

include_directories(/usr/local/include)
include_directories(/usr/local/opt/openssl/include)
link_directories(/usr/local/lib)
link_directories(/usr/local/opt/openssl/lib)
link_directories(/usr/local/Cellar/libuv/1.11.0/lib)

endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin") 


add_executable(ExtendedKF src/main.cpp ${sources})
target_link_libraries(ExtendedKF z ssl uv uWS)

#add_executable(ExtendedKF_test src/main2.cpp ${sources})
