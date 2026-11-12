# FluxCell CLI Tools
add_executable(spm_cli
    spm_cli.cpp
)

target_link_libraries(spm_cli PRIVATE
    fluxcell_core
)

set_target_properties(spm_cli PROPERTIES
    FOLDER "tools"
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
