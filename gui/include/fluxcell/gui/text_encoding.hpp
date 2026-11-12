/// @file text_encoding.hpp
/// @brief Windows 窄字符文件 API 的路径编码转换
///
/// GUI 内部统一 UTF-8；Windows 的 std::ifstream/ofstream 与文件对话框
/// 返回/接受本地代码页（ANSI）窄字符路径，读写文件前需转换。
/// 其余平台路径即 UTF-8，两个函数原样直通。
#pragma once

#include <string>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX // 防止 windows.h 的 min/max 宏污染 std::max/std::min
#endif
#include <windows.h>
#endif

namespace fluxcell::gui {

#ifdef _WIN32

/// UTF-8 → 本地代码页（Windows 窄字符 API 如 std::ifstream 以 ANSI 解析路径）
inline std::string utf8_to_acp(const std::string& utf8) {
    if (utf8.empty()) {
        return utf8;
    }
    const int wlen = MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(),
                                         static_cast<int>(utf8.size()),
                                         nullptr, 0);
    std::wstring wide(static_cast<std::size_t>(wlen), L'\0');
    MultiByteToWideChar(CP_UTF8, 0, utf8.c_str(),
                        static_cast<int>(utf8.size()), wide.data(), wlen);
    const int alen = WideCharToMultiByte(CP_ACP, 0, wide.c_str(), wlen,
                                         nullptr, 0, nullptr, nullptr);
    std::string ansi(static_cast<std::size_t>(alen), '\0');
    WideCharToMultiByte(CP_ACP, 0, wide.c_str(), wlen, ansi.data(), alen,
                        nullptr, nullptr);
    return ansi;
}

/// 本地代码页 → UTF-8（文件对话框返回的路径转 GUI 内部编码）
inline std::string acp_to_utf8(const std::string& ansi) {
    if (ansi.empty()) {
        return ansi;
    }
    const int wlen = MultiByteToWideChar(CP_ACP, 0, ansi.c_str(),
                                         static_cast<int>(ansi.size()),
                                         nullptr, 0);
    std::wstring wide(static_cast<std::size_t>(wlen), L'\0');
    MultiByteToWideChar(CP_ACP, 0, ansi.c_str(),
                        static_cast<int>(ansi.size()), wide.data(), wlen);
    const int ulen = WideCharToMultiByte(CP_UTF8, 0, wide.c_str(), wlen,
                                         nullptr, 0, nullptr, nullptr);
    std::string utf8(static_cast<std::size_t>(ulen), '\0');
    WideCharToMultiByte(CP_UTF8, 0, wide.c_str(), wlen, utf8.data(), ulen,
                        nullptr, nullptr);
    return utf8;
}

#else

inline std::string utf8_to_acp(const std::string& s) { return s; }
inline std::string acp_to_utf8(const std::string& s) { return s; }

#endif

} // namespace fluxcell::gui
