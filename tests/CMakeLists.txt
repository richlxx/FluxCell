# FluxCell Tests
add_executable(fluxcell_tests
    test_physics.cpp
    test_discretization.cpp
    test_spm.cpp
    test_spme.cpp
    test_dfn.cpp
    test_ida_solver.cpp
    test_io.cpp
)

target_link_libraries(fluxcell_tests PRIVATE
    fluxcell_core
    GTest::gtest_main
)

# IO 测试需要读取源码树中的 params/*.yml
target_compile_definitions(fluxcell_tests PRIVATE
    FLUXCELL_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

include(GoogleTest)
gtest_discover_tests(fluxcell_tests)

set_target_properties(fluxcell_tests PROPERTIES FOLDER "tests")
