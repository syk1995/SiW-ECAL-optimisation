#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <iostream>

namespace py = pybind11;
py::tuple Event_Selection(std::vector<std::vector<float>>& E,
                          std::vector<std::vector<float>>& PosX,
                          std::vector<std::vector<float>>& PosY,
                          std::vector<std::vector<float>>& PosZ,
                          std::vector<float>& E_truth) {

    size_t total = E.size();
    std::cout << "[Event_Selection] Total events before filtering: " << total << std::endl;
    std::vector<std::vector<float>> E_new, X_new, Y_new, Z_new;
    std::vector<float> E_truth_new;
    for (size_t i = 0; i < total; ++i) {
        if (!E[i].empty()) {
            E_new.push_back(std::move(E[i]));
            X_new.push_back(std::move(PosX[i]));
            Y_new.push_back(std::move(PosY[i]));
            Z_new.push_back(std::move(PosZ[i]));
            E_truth_new.push_back(E_truth[i]);
        }
    }
    std::cout << "[Event_Selection] After filtering: " << E_new.size() << " events remain." << std::endl;
    std::cout << "[Event_Selection] Aligning Z positions..." << std::endl;
    #pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < static_cast<long>(Z_new.size()); ++i) {
        if (Z_new[i].empty()) continue;
        float min_z = *std::min_element(Z_new[i].begin(), Z_new[i].end());
        for (auto& z : Z_new[i]) z -= min_z;
    }

    std::cout << "[Event_Selection] Z alignment completed." << std::endl;
    return py::make_tuple(E_new, X_new, Y_new, Z_new, E_truth_new);
}

// pybind11 模块
PYBIND11_MODULE(event_selection, m) {
    m.doc() = "Event selection: remove empty events + Z alignment, return ExyzEtruth";
    m.def("Event_Selection", &Event_Selection, py::arg("E"), py::arg("PosX"),
          py::arg("PosY"), py::arg("PosZ"), py::arg("E_truth"));
}
