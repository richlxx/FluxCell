/// @file plot_util.hpp
/// @brief ImPlot 绘图小工具（跨面板共享）
#pragma once

#include <fluxcell/core/types.hpp>

#include <implot.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace fluxcell::gui {

/// y 轴替代 AutoFit：各系列数据范围上下外扩 8% 后固定。
/// （ImPlot 的 AutoFit 轴紧贴数据极值，曲线会顶到图表上下边框）
template <typename... Series>
inline void fit_y_padded(const Series&... series) {
    const std::pair<const real_t*, int> all[] = {
        {series.data(), static_cast<int>(series.size())}...};
    real_t mn = std::numeric_limits<real_t>::max();
    real_t mx = std::numeric_limits<real_t>::lowest();
    for (const auto& [p, n] : all) {
        for (int i = 0; i < n; ++i) {
            const real_t v = p[i];
            if (!std::isfinite(v)) {
                continue;
            }
            mn = std::min(mn, v);
            mx = std::max(mx, v);
        }
    }
    if (mn > mx) {
        return; // 无有效数据，保留默认范围
    }
    real_t pad = (mx - mn) * 0.08f;
    if (pad <= 0.0f) {
        pad = std::max(std::abs(mx), static_cast<real_t>(1.0))
              * static_cast<real_t>(0.1); // 常量曲线兜底
    }
    ImPlot::SetupAxisLimits(ImAxis_Y1, mn - pad, mx + pad,
                            ImPlotCond_Always);
}

} // namespace fluxcell::gui
