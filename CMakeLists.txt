cmake_minimum_required(VERSION 3.25)
project(FluxCell VERSION 0.1.0 LANGUAGES CXX)

# ── C++ Standard ─────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Output Directories ───────────────────────────────────────
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ── Options ──────────────────────────────────────────────────
option(FLUXCELL_BUILD_TESTS "Build unit tests" ON)
option(FLUXCELL_BUILD_GUI   "Build GUI application" ON)

# ── Compiler Options ─────────────────────────────────────────
if(MSVC)
    # /utf-8      源码字符集（中文注释）
    # /W3         警告级别
    # /MP         多进程并行编译，加速 VS2022 全量构建
    # /permissive- 严格标准符合性
    add_compile_options(/utf-8 /W3 /MP /permissive-)
endif()

# ── VS2022 Solution Layout ───────────────────────────────────
# 解决方案资源管理器按文件夹分组
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
# F5 直接启动 GUI 应用（默认启动项目）
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT fluxcell_app)

# ── Dependencies ─────────────────────────────────────────────
include(cmake/FetchDependencies.cmake)

# ── Sub-directories ──────────────────────────────────────────
add_subdirectory(core)
add_subdirectory(tools)

if(FLUXCELL_BUILD_GUI)
    add_subdirectory(gui)
    add_subdirectory(app)
endif()

if(FLUXCELL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ── Install / Packaging ─────────────────────────────────────────────
# 生成可分发的安装目录（可执行文件 + 运行时 DLL + 默认参数文件）：
#   cmake --install build --config Release --prefix dist/FluxCell
# 或一键 ZIP 便携包：
#   cpack -C Release -G ZIP（输出于 build/ 目录）
# 说明：
#   - 静态库（fluxcell_core/gui）与测试不参与安装，仅安装运行时产物
#   - vcpkg x64-windows 依赖为 shared DLL（hdf5/sundials/glfw3/yaml-cpp 等），
#     $<TARGET_RUNTIME_DLLS:...> 自动解析传递闭包（含 zlib 等间接依赖）
#   - MSVC /MD 动态运行库需目标机安装 VC Redist，或随 bin 附 vcruntime140 等
if(FLUXCELL_BUILD_GUI)
    install(TARGETS fluxcell_app RUNTIME DESTINATION bin)
    install(FILES $<TARGET_RUNTIME_DLLS:fluxcell_app> DESTINATION bin)
endif()
install(TARGETS spm_cli RUNTIME DESTINATION bin)
if(NOT FLUXCELL_BUILD_GUI)
    # CLI-only 场景：DLL 依赖由 spm_cli 的链接闭包解析（不含 GLFW）
    install(FILES $<TARGET_RUNTIME_DLLS:spm_cli> DESTINATION bin)
endif()
# 默认参数文件置于安装根目录 params/（GUI/CLI 以相对路径加载）
install(DIRECTORY params/ DESTINATION params)

# ── CPack（ZIP 便携包元数据）───────────────────────────────────────
set(CPACK_PACKAGE_NAME "FluxCell")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_VENDOR "FluxCell")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "FluxCell 锂离子电池电化学仿真软件（SPM/SPMe/DFN，含 SUNDIALS IDA 自适应求解）")
set(CPACK_GENERATOR "ZIP")
set(CPACK_SYSTEM_NAME "win64")
set(CPACK_PACKAGE_FILE_NAME
    "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_SYSTEM_NAME}")
include(CPack)

# ── VS 元目标归组（ALL_BUILD / ZERO_CHECK / INSTALL / RUN_TESTS）──
foreach(meta_target ALL_BUILD ZERO_CHECK INSTALL RUN_TESTS)
    if(TARGET ${meta_target})
        set_target_properties(${meta_target} PROPERTIES FOLDER "CMake")
    endif()
endforeach()
