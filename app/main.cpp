/// @file main.cpp
/// @brief FluxCell 应用程序入口
#include <fluxcell/gui/application.hpp>
#include <fluxcell/gui/main_window.hpp>

#include <iostream>
#include <memory>
#include <string>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>

namespace {

/// UTF-8 消息转宽字符（MessageBoxW 需要，避免本地代码页乱码）
std::wstring utf8_to_wide(const std::string& text) {
    if (text.empty()) {
        return {};
    }
    const int len = MultiByteToWideChar(CP_UTF8, 0, text.c_str(),
                                        static_cast<int>(text.size()),
                                        nullptr, 0);
    std::wstring result(len, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, text.c_str(), static_cast<int>(text.size()),
                        result.data(), len);
    return result;
}

} // namespace
#endif

int main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;

    try {
        // 创建 GUI 应用
        fluxcell::gui::Application app("FluxCell — DFN Battery Simulator", 1600, 900);

        // 创建主窗口并注入渲染回调与退出回调（菜单“退出”→ Application）
        auto main_window = std::make_unique<fluxcell::gui::MainWindow>();
        main_window->set_exit_callback([&app]() { app.request_exit(); });
        app.set_render_callback([&main_window]() { main_window->render(); });

        // 运行主循环
        app.run();

    } catch (const std::exception& e) {
#ifdef _WIN32
        // GUI 子系统无控制台，致命错误弹窗提示
        MessageBoxW(nullptr, utf8_to_wide(e.what()).c_str(),
                    L"FluxCell 启动失败", MB_ICONERROR | MB_OK);
#else
        std::cerr << "Fatal error: " << e.what() << std::endl;
#endif
        return 1;
    }

    return 0;
}
