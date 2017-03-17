project(RigidBodyTest)

set(SOURCES
  RigidBodyTest.cpp
)

cuda_add_executable(RigidBodyTest ${SOURCES})
target_link_libraries(RigidBodyTest ${CUDA_LIBRARIES})
install(TARGETS RigidBodyTest DESTINATION bin)
