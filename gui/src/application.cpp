/// @file application.cpp
/// @brief Dear ImGui 应用程序框架实现
#include <fluxcell/gui/application.hpp>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <GLFW/glfw3.h>

#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace fluxcell::gui {

// ══════════════════════════════════════════════════════════════
//  GLFW 错误回调
// ══════════════════════════════════════════════════════════════

static void glfw_error_callback(int error, const char* description) {
    fprintf(stderr, "GLFW Error %d: %s\n", error, description);
}

// ══════════════════════════════════════════════════════════════
//  窗口图标（Win32 资源提取，来源 app.rc 中嵌入的 fluxcell.ico）
// ══════════════════════════════════════════════════════════════

#ifdef _WIN32

/// 从 exe 图标资源提取最接近 requested 尺寸的一档并转为 RGBA 像素
static bool load_icon_pixels(int requested, int& width, int& height,
                             std::vector<unsigned char>& pixels) {
    const HICON icon = static_cast<HICON>(LoadImageW(
        GetModuleHandleW(nullptr), MAKEINTRESOURCEW(1), IMAGE_ICON,
        requested, requested, LR_DEFAULTCOLOR));
    if (!icon) {
        return false;
    }

    bool ok = false;
    ICONINFO info{};
    if (GetIconInfo(icon, &info) && info.hbmColor) {
        BITMAP bm{};
        if (GetObjectW(info.hbmColor, sizeof(bm), &bm) > 0 &&
            bm.bmBitsPixel == 32 && bm.bmWidth > 0 && bm.bmHeight > 0) {
            BITMAPINFOHEADER bi{};
            bi.biSize        = sizeof(bi);
            bi.biWidth       = bm.bmWidth;
            bi.biHeight      = -bm.bmHeight; // top-down
            bi.biPlanes      = 1;
            bi.biBitCount    = 32;
            bi.biCompression = BI_RGB;
            pixels.resize(static_cast<size_t>(bm.bmWidth) * bm.bmHeight * 4);
            const HDC dc = GetDC(nullptr);
            if (dc) {
                ok = GetDIBits(dc, info.hbmColor, 0, bm.bmHeight, pixels.data(),
                               reinterpret_cast<BITMAPINFO*>(&bi),
                               DIB_RGB_COLORS) != 0;
                ReleaseDC(nullptr, dc);
            }
            if (ok) {
                // 位图为 BGRA，调换 R/B 通道给 OpenGL/GLFW 用
                for (size_t i = 0; i + 3 < pixels.size(); i += 4) {
                    std::swap(pixels[i], pixels[i + 2]);
                }
                width  = bm.bmWidth;
                height = bm.bmHeight;
            }
        }
    }
    if (info.hbmColor) {
        DeleteObject(info.hbmColor);
    }
    if (info.hbmMask) {
        DeleteObject(info.hbmMask);
    }
    DestroyIcon(icon);
    return ok;
}

/// 设置 GLFW 窗口/任务栏图标（多档尺寸交给系统按 DPI 挑选）
static void setup_window_icon(GLFWwindow* window) {
    static constexpr int kSizes[] = {16, 24, 32, 48};
    GLFWimage images[4];
    std::vector<std::vector<unsigned char>> buffers; // 数据存活至 glfwSetWindowIcon 拷贝
    buffers.reserve(std::size(kSizes));
    int count = 0;
    for (int size : kSizes) {
        std::vector<unsigned char> pixels;
        int w = 0;
        int h = 0;
        if (load_icon_pixels(size, w, h, pixels)) {
            images[count].width  = w;
            images[count].height = h;
            images[count].pixels = pixels.data();
            buffers.push_back(std::move(pixels));
            ++count;
        }
    }
    if (count > 0) {
        glfwSetWindowIcon(window, count, images);
    }
}

#endif // _WIN32

// ══════════════════════════════════════════════════════════════
//  构造/析构
// ══════════════════════════════════════════════════════════════

Application::Application(const std::string& title, int width, int height)
    : title_(title) {
    glfwSetErrorCallback(glfw_error_callback);

    if (!glfwInit()) {
        throw std::runtime_error("Failed to initialize GLFW");
    }

    // OpenGL 3.3 + Core Profile
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window_ = glfwCreateWindow(width, height, title.c_str(), nullptr, nullptr);
    if (!window_) {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }

    // 窗口/任务栏图标（exe 内嵌的 fluxcell.ico 资源）
#ifdef _WIN32
    setup_window_icon(window_);
#endif

    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // VSync

    setup_imgui();
}

Application::~Application() {
    teardown_imgui();
    if (window_) {
        glfwDestroyWindow(window_);
    }
    glfwTerminate();
}

// ══════════════════════════════════════════════════════════════
//  ImGui 初始化/清理
// ══════════════════════════════════════════════════════════════

void Application::setup_imgui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;    // 启用停靠
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // 键盘导航

    // 字体：加载系统中文字体并使用较大字号，保证中文显示与可读性
    load_fonts(io);

    // 样式：亮色主题
    ImGui::StyleColorsLight();

    // 后端初始化
    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL3_Init("#version 330");
}

void Application::load_fonts(ImGuiIO& io) {
    // 候选字体按优先级排列：微软雅黑 → 雅黑 Light → 黑体 → 宋体
    static const char* font_candidates[] = {
        "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/msyhl.ttc",
        "C:/Windows/Fonts/simhei.ttf",
        "C:/Windows/Fonts/simsun.ttc",
    };

    constexpr float kFontSize = 18.0f; // 默认 13px，加大以提升可读性
    namespace fs = std::filesystem;

    // 合并字形范围：基本拉丁（含 ° ± ×）+ 希腊字母（η φ κ α Δ 等物理符号）
    // + 数学运算符（− ≤ ≈ 等）+ 中文全量汉字。
    // 注意：ChineseSimplifiedCommon 仅含 ~2500 常用字，界面文案中的次常用字
    //（如“幂”）会渲染为问号，必须用 ChineseFull（全部 GB 汉字）。
    // ranges 在首次 atlas Build（第一次 NewFrame）前必须保持有效，
    // 因此用 static 存储而非局部变量。
    static ImVector<ImWchar> glyph_ranges;
    if (glyph_ranges.empty()) {
        static const ImWchar extra_ranges[] = {
            0x0370, 0x03FF, // 希腊字母
            0x2100, 0x214F, // Letterlike Symbols（℃ 等）
            0x2200, 0x22FF, // Mathematical Operators（− × ≈ ≤ 等）
            0x25A0, 0x25FF, // Geometric Shapes（▲ ▼ 排序按钮等）
            0,
        };
        ImFontGlyphRangesBuilder builder;
        builder.AddRanges(io.Fonts->GetGlyphRangesDefault());
        builder.AddRanges(extra_ranges);
        builder.AddRanges(io.Fonts->GetGlyphRangesChineseFull());
        builder.BuildRanges(&glyph_ranges);
    }

    ImFont* loaded = nullptr;
    for (const char* path : font_candidates) {
        if (!fs::exists(path)) {
            continue;
        }
        loaded = io.Fonts->AddFontFromFileTTF(path, kFontSize, nullptr, glyph_ranges.Data);
        if (loaded) {
            break;
        }
    }

    if (!loaded) {
        // 兜底：默认字体路径 + 同样的字形范围（系统中文字体均不可用时仅西文可显示）
        loaded = io.Fonts->AddFontFromFileTTF(
            "C:/Windows/Fonts/arial.ttf", kFontSize, nullptr, glyph_ranges.Data);
    }
    (void)loaded;
}

void Application::teardown_imgui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

// ══════════════════════════════════════════════════════════════
//  帧管理
// ══════════════════════════════════════════════════════════════

void Application::begin_frame() {
    glfwPollEvents();
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void Application::end_frame() {
    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(window_, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.84f, 0.84f, 0.86f, 1.0f); // 与亮色主题协调的浅灰背景
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window_);
}

// ══════════════════════════════════════════════════════════════
//  主循环
// ══════════════════════════════════════════════════════════════

void Application::run() {
    running_ = true;
    while (!should_close()) {
        begin_frame();

        if (render_callback_) {
            render_callback_();
        } else {
            ImGui::ShowDemoWindow();
        }

        end_frame();
    }
}

void Application::request_exit() {
    running_ = false;
}

bool Application::should_close() const {
    return !running_ || glfwWindowShouldClose(window_);
}

} // namespace fluxcell::gui
