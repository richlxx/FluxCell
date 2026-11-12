# FluxCell Core Library
# 头文件显式列出：VS2022 中可见、可搜索（Ctrl+T）、可跳转（ClInclude）
add_library(fluxcell_core STATIC
    include/fluxcell/core/types.hpp
    include/fluxcell/core/model/cell_model.hpp
    include/fluxcell/core/model/spm_model.hpp
    include/fluxcell/core/model/spme_model.hpp
    include/fluxcell/core/model/dfn_model.hpp
    include/fluxcell/core/physics/physics.hpp
    include/fluxcell/core/discretization/discretization.hpp
    include/fluxcell/core/simulation/simulation.hpp
    include/fluxcell/core/solver/dae_problem.hpp
    include/fluxcell/core/io/io.hpp

    src/spm_model.cpp
    src/spme_model.cpp
    src/dfn_model.cpp
    src/simulation.cpp
    src/physics.cpp
    src/discretization.cpp
    src/io.cpp
)

if(FLUXCELL_WITH_SUNDIALS)
    # SUNDIALS IDA 自适应求解器（仅启用时编译进库）
    target_sources(fluxcell_core PRIVATE
        include/fluxcell/core/solver/ida_solver.hpp
        src/ida_solver.cpp
    )
endif()

target_include_directories(fluxcell_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fluxcell_core PUBLIC
    Eigen3::Eigen
    yaml-cpp::yaml-cpp
    nlohmann_json::nlohmann_json
    hdf5::hdf5-shared
)

if(FLUXCELL_WITH_SUNDIALS)
    target_link_libraries(fluxcell_core PUBLIC ${FLUXCELL_SUNDIALS_LIBS})
    target_compile_definitions(fluxcell_core PUBLIC FLUXCELL_WITH_SUNDIALS=1)
    if(FLUXCELL_WITH_SUNDIALS_KLU)
        target_compile_definitions(fluxcell_core PUBLIC FLUXCELL_WITH_SUNDIALS_KLU=1)
    endif()
endif()

# C++20
target_compile_features(fluxcell_core PUBLIC cxx_std_20)

set_target_properties(fluxcell_core PROPERTIES FOLDER "fluxcell")
