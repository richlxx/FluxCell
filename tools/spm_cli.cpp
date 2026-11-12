/// @file spm_cli.cpp
/// @brief SPM 命令行仿真工具：恒流放电并导出 CSV（用于 PyBaMM 基准对比）
///
/// 用法: spm_cli <params.yml> <current_A> <dt_s> <output.csv> [max_time_s] [voltage_min]
#include <fluxcell/core/io/io.hpp>
#include <fluxcell/core/model/spm_model.hpp>

#include <cmath>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: spm_cli <params.yml> <current_A> <dt_s> <output.csv>"
                     " [max_time_s=7200] [voltage_min=2.5]\n";
        return 2;
    }

    try {
        const std::string params_path = argv[1];
        const fluxcell::real_t current = std::stod(argv[2]);
        const fluxcell::real_t dt = std::stod(argv[3]);
        const std::string output_path = argv[4];
        const fluxcell::real_t max_time = argc > 5 ? std::stod(argv[5]) : 7200.0;
        const fluxcell::real_t voltage_min = argc > 6 ? std::stod(argv[6]) : 2.5;

        auto params = fluxcell::io::load_spm_parameters(params_path);
        fluxcell::model::SPMModel model(params, current);

        std::vector<fluxcell::SolutionSnapshot> snapshots;
        snapshots.push_back(model.snapshot(current));

        while (model.time() < max_time && model.voltage() > voltage_min) {
            model.step(current, dt);
            snapshots.push_back(model.snapshot(current));
        }

        fluxcell::CellGeometry geom; // CSV 导出不使用几何字段
        fluxcell::io::export_csv(output_path, snapshots, geom);

        std::cout << "steps: " << snapshots.size()
                  << "  t_end: " << model.time() << " s"
                  << "  V_end: " << model.voltage() << " V"
                  << "  -> " << output_path << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
