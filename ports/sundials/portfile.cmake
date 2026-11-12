vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO LLNL/sundials
    REF "v${VERSION}"
    SHA512 e026a7b8d3de8a4a53de9c3b4bdd3cdb1fc61feb2b9a31e20017474b36964216ddd8a84822d2111f5ca63f129b7bbad278f534968cef09f973b92d3f24c675c2
    HEAD_REF master
)

string(COMPARE EQUAL "${VCPKG_LIBRARY_LINKAGE}" "static" SUN_BUILD_STATIC)
string(COMPARE EQUAL "${VCPKG_LIBRARY_LINKAGE}" "dynamic" SUN_BUILD_SHARED)

if(VCPKG_TARGET_IS_ANDROID)
    set(POSIX_TIMERS "-DSUNDIALS_POSIX_TIMERS=TRUE")
endif()

# FluxCell overlay：无条件启用 KLU 稀疏线性求解器（suitesparse-klu 为主
# 依赖）。上游 portfile 不开 SUNDIALS_ENABLE_KLU，且本 vcpkg 版本的
# portfile 上下文拿不到 FEATURES 变量，故不走 feature 开关。
# 同时注入自定义 FindKLU：vcpkg suitesparse-klu 的 KLUConfig 仅导出
# imported target，不提供 SUNDIALS 查找逻辑所需的
# KLU_LIBRARIES / KLU_INCLUDE_DIR 变量
vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        ${POSIX_TIMERS}
        -D_BUILD_EXAMPLES=OFF
        -DSUNDIALS_TEST_UNITTESTS=OFF
        -DBUILD_STATIC_LIBS=${SUN_BUILD_STATIC}
        -DBUILD_SHARED_LIBS=${SUN_BUILD_SHARED}
        -DSUNDIALS_ENABLE_KLU=ON
        # 注：不注入自定义 FindKLU——SUNDIALS 自带 cmake/tpl/FindKLU.cmake
        # （位于自身 CMAKE_MODULE_PATH），走 config 模式命中 vcpkg
        # suitesparse-klu 的 KLUConfig，创建 SUNDIALS::KLU alias；
        # 自定义模块反而会遮蔽官方模块并丢失 target 创建逻辑
)

vcpkg_cmake_install()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
file(REMOVE "${CURRENT_PACKAGES_DIR}/LICENSE")
file(REMOVE "${CURRENT_PACKAGES_DIR}/debug/LICENSE")

vcpkg_copy_pdbs()
vcpkg_cmake_config_fixup(CONFIG_PATH "lib/cmake/${PORT}")
