# FluxCell Application（app.rc：应用图标资源，fluxcell.ico 由根目录 Simulator.png 生成）
add_executable(fluxcell_app main.cpp app.rc)

target_link_libraries(fluxcell_app PRIVATE
    fluxcell_core
    fluxcell_gui
)

# GUI 子系统（不弹控制台窗口）；标准 main() 入口经链接器显式指定，免去 WinMain 适配
target_link_options(fluxcell_app PRIVATE "/ENTRY:mainCRTStartup")

set_target_properties(fluxcell_app PROPERTIES
    OUTPUT_NAME "FluxCell"
    WIN32_EXECUTABLE TRUE
    FOLDER "fluxcell"
    # F5 调试时工作目录 = 项目根，params/*.yml 相对路径可直接加载
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
