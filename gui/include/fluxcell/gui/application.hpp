/// @file application.hpp
/// @brief Dear ImGui 应用程序框架
#pragma once

#include <fluxcell/core/types.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

struct GLFWwindow;
struct ImGuiIO; // imgui.h 中定义，仅前向声明避免引依赖

namespace fluxcell::gui {

class MainWindow;

/// ImGui 应用程序生命周期管理
class Application {
public:
    Application(const std::string& title, int width, int height);
    ~Application();

    Application(const Application&) = delete;
    Application& operator=(const Application&) = delete;

    /// 运行主循环
    void run();

    /// 注册每帧渲染回调（替代默认 Demo 窗口）
    using RenderCallback = std::function<void()>;
    void set_render_callback(RenderCallback cb) { render_callback_ = std::move(cb); }

    /// 请求退出
    void request_exit();

    /// 窗口是否应该关闭
    bool should_close() const;

    GLFWwindow* window() { return window_; }

private:
    void setup_imgui();
    void load_fonts(ImGuiIO& io); ///< 加载系统中文字体（较大字号）
    void teardown_imgui();
    void begin_frame();
    void end_frame();

    GLFWwindow* window_ = nullptr;
    bool running_ = true;
    std::string title_;
    RenderCallback render_callback_;
};

} // namespace fluxcell::gui
