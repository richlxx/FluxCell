include(FetchContent)

# ── Dear ImGui (from source) ─────────────────────────────────
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        v1.91.8-docking
    GIT_SHALLOW    TRUE
)

# ── ImPlot (from source) ─────────────────────────────────────
FetchContent_Declare(
    implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        v0.16
    GIT_SHALLOW    TRUE
)

# ── ImGuiFileDialog (from source) ────────────────────────────
# SOURCE_SUBDIR 指向不存在的目录：仅获取源码，不执行上游 CMakeLists.txt
# （上游 CMake 未感知 imgui include 路径，我们自行编译 ImGuiFileDialog.cpp）
FetchContent_Declare(
    imguifiledialog
    GIT_REPOSITORY https://github.com/aiekick/ImGuiFileDialog.git
    GIT_TAG        v0.6.7
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  upstream-cmake-disabled
)

# ── Google Test ──────────────────────────────────────────────
# 测试专用（静态构建），禁用上游 install 规则，避免 gtest 头/库
# 混入项目安装产物（include/gtest、lib/gtest 等）
set(INSTALL_GTEST OFF CACHE INTERNAL "Disable googletest install rules")
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG        v1.15.2
    GIT_SHALLOW    TRUE
)

# ── Fetch all dependencies ───────────────────────────────────
FetchContent_MakeAvailable(imgui implot imguifiledialog googletest)

# vcpkg: GLFW3 (target: glfw)
find_package(glfw3 CONFIG REQUIRED)

# ══════════════════════════════════════════════════════════════
#  Build ImGui as a static library
# ══════════════════════════════════════════════════════════════
add_library(imgui_lib STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui_lib PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui_lib PUBLIC glfw)
find_package(OpenGL REQUIRED)
target_link_libraries(imgui_lib PUBLIC OpenGL::GL)
set_target_properties(imgui_lib PROPERTIES FOLDER "third_party")

# ══════════════════════════════════════════════════════════════
#  Build ImPlot as a static library
# ══════════════════════════════════════════════════════════════
add_library(implot_lib STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
    ${implot_SOURCE_DIR}/implot_demo.cpp
)
target_include_directories(implot_lib PUBLIC ${implot_SOURCE_DIR})
target_link_libraries(implot_lib PUBLIC imgui_lib)
set_target_properties(implot_lib PROPERTIES FOLDER "third_party")

# ══════════════════════════════════════════════════════════════
#  ImGuiFileDialog (header-only usage, compiled as static lib)
# ══════════════════════════════════════════════════════════════
add_library(imgui_file_dialog STATIC
    ${imguifiledialog_SOURCE_DIR}/ImGuiFileDialog.cpp
)
target_include_directories(imgui_file_dialog PUBLIC ${imguifiledialog_SOURCE_DIR})
target_link_libraries(imgui_file_dialog PUBLIC imgui_lib)
set_target_properties(imgui_file_dialog PROPERTIES FOLDER "third_party")

# ══════════════════════════════════════════════════════════════
#  Find vcpkg packages
# ══════════════════════════════════════════════════════════════
find_package(Eigen3 CONFIG REQUIRED)
find_package(yaml-cpp CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(HDF5 CONFIG REQUIRED)

# vcpkg: SUNDIALS（IDA 自适应求解器；shared 库，SUNDIALS 7.x API）
option(FLUXCELL_WITH_SUNDIALS "Enable SUNDIALS IDA adaptive solver" ON)
if(FLUXCELL_WITH_SUNDIALS)
    find_package(SUNDIALS CONFIG REQUIRED)
    set(FLUXCELL_SUNDIALS_LIBS
        SUNDIALS::ida_shared
        SUNDIALS::nvecserial_shared
        SUNDIALS::sunmatrixdense_shared
        SUNDIALS::sunlinsoldense_shared
        SUNDIALS::sunlinsolspgmr_shared
    )
    # sundials[klu]（项目 overlay port）：大稀疏 DAE 系统用 KLU 稀疏直解
    # 替代无预条件 GMRES（厚膜长循环崩溃根治，详见 tests/test_ida_solver.cpp）
    if(TARGET SUNDIALS::sunlinsolklu_shared)
        list(APPEND FLUXCELL_SUNDIALS_LIBS
            SUNDIALS::sunmatrixsparse_shared
            SUNDIALS::sunlinsolklu_shared)
        set(FLUXCELL_WITH_SUNDIALS_KLU TRUE)
        message(STATUS "SUNDIALS KLU sparse solver: enabled")
    else()
        message(STATUS "SUNDIALS KLU sparse solver: NOT available (falling back to GMRES)")
    endif()
endif()

# googletest 由 FetchContent 创建的目标归入 third_party 文件夹
foreach(gtest_target gtest gtest_main gmock gmock_main)
    if(TARGET ${gtest_target})
        set_target_properties(${gtest_target} PROPERTIES FOLDER "third_party/googletest")
    endif()
endforeach()
