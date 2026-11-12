# FluxCell GUI Library
# 头文件显式列出：VS2022 中可见、可搜索（Ctrl+T）、可跳转（ClInclude）
add_library(fluxcell_gui STATIC
    include/fluxcell/gui/application.hpp
    include/fluxcell/gui/main_window.hpp
    include/fluxcell/gui/plot_panel.hpp
    include/fluxcell/gui/simulation_panel.hpp

    src/application.cpp
    src/main_window.cpp
    src/plot_panel.cpp
    src/simulation_panel.cpp
)

target_include_directories(fluxcell_gui PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(fluxcell_gui PUBLIC
    fluxcell_core
    imgui_lib
    implot_lib
    imgui_file_dialog
)

target_compile_features(fluxcell_gui PUBLIC cxx_std_20)

set_target_properties(fluxcell_gui PROPERTIES FOLDER "fluxcell")
