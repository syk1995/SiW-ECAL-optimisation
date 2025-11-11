// ReadRootEvents.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <algorithm>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "nanoflann.hpp"

namespace py = pybind11;

struct Event {
    std::vector<float> energy;
    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<float> pos_z;
    std::vector<std::pair<int, int>> edge_index;
    float E_truth;
};

// ---- 自定义点云结构 ----
template <typename T>
struct PointCloud {
    struct Point {
        T x, y, z;
    };
    std::vector<Point> pts;

    // nanoflann接口: 返回点的数量
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    // nanoflann接口: 返回第idx个点的第dim维坐标
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].x;
        else if (dim == 1) return pts[idx].y;
        else return pts[idx].z;
    }

    // 可选: bounding-box优化 (不需要时可返回 false)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using namespace nanoflann;

std::vector<std::pair<int,int>> build_knn_graph(const Event &ev, int k=3) {
    PointCloud<float> cloud;
    for (size_t i=0; i<ev.pos_x.size(); ++i)
        cloud.pts.push_back({ev.pos_x[i], ev.pos_y[i], ev.pos_z[i]});
    typedef KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<float, PointCloud<float>>,
        PointCloud<float>,
        3 /*维度*/
    > kd_tree_t;
    kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();
    std::vector<std::pair<int,int>> edges;
    std::vector<unsigned int> idx(k+1);
    std::vector<float> dist(k+1);
    for (size_t i=0; i<cloud.pts.size(); ++i) {
        index.knnSearch(&cloud.pts[i].x, k+1, &idx[0], &dist[0]);
        for (int j=1; j<=k; ++j) {  // 跳过自身
            edges.emplace_back(i, idx[j]);
        }
    }
    return edges;
}

// Sort nodes of Event by pos_z
void sort_event_by_z(Event &ev) {
    std::vector<size_t> idx(ev.pos_z.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return ev.pos_z[a] < ev.pos_z[b];
    });
    std::vector<float> energy_new, x_new, y_new, z_new;
    for (size_t i : idx) {
        energy_new.push_back(ev.energy[i]);
        x_new.push_back(ev.pos_x[i]);
        y_new.push_back(ev.pos_y[i]);
        z_new.push_back(ev.pos_z[i]);
    }
    ev.energy = std::move(energy_new);
    ev.pos_x = std::move(x_new);
    ev.pos_y = std::move(y_new);
    ev.pos_z = std::move(z_new);
}


std::vector<Event> ReadRoot(const std::string &filename, float threshold = 0.0,int KNN_K=8) {
    std::vector<Event> events;

    TFile file(filename.c_str(), "READ");
    if (!file.IsOpen()) throw std::runtime_error("Cannot open ROOT file");

    TTree *tree = (TTree*)file.Get("events");
    if (!tree) throw std::runtime_error("Cannot find TTree 'events'");
    // Branch variables
    std::vector<float> *energy = nullptr;
    std::vector<float> *pos_x  = nullptr;
    std::vector<float> *pos_y  = nullptr;
    std::vector<float> *pos_z  = nullptr;
    std::vector<float> *E_truth = nullptr;
    tree->SetBranchAddress("simplecaloRO.energy", &energy);
    tree->SetBranchAddress("simplecaloRO.position.x", &pos_x);
    tree->SetBranchAddress("simplecaloRO.position.y", &pos_y);
    tree->SetBranchAddress("simplecaloRO.position.z", &pos_z);
    tree->SetBranchAddress("MCParticles.p0", &E_truth);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; ++i) {
        tree->GetEntry(i);
        Event ev;
        if(E_truth->size() != 1) {
            throw std::runtime_error("E_truth size is not 1 for event " + std::to_string(i));
            continue;
        }
        ev.E_truth = (*E_truth)[0];
        for (size_t j=0; j<energy->size(); ++j) {
            if ((*energy)[j] * 1000 > threshold) {// Convert GeV to MeV
                ev.energy.push_back((*energy)[j]);
                ev.pos_x.push_back((*pos_x)[j]);
                ev.pos_y.push_back((*pos_y)[j]);
                ev.pos_z.push_back((*pos_z)[j]);
            }
        }
        if (ev.energy.empty()) continue; // skip empty events
        float zmin = *std::min_element(ev.pos_z.begin(), ev.pos_z.end());
        for (auto &z : ev.pos_z) z -= zmin;
        //precision 
        for (size_t i = 0; i < ev.pos_z.size(); ++i) {
            ev.pos_z[i] = std::round(ev.pos_z[i] * 1000.0) / 1000.0;
            ev.E_truth = std::round(ev.E_truth * 1000.0) / 1000.0;
        }
        sort_event_by_z(ev);
        ev.edge_index = build_knn_graph(ev, KNN_K);
        events.push_back(ev);
    }
    file.Close();
    return events;
}

PYBIND11_MODULE(ReadRoot, m) {
    py::class_<Event>(m, "Event")
        .def_readonly("energy", &Event::energy)
        .def_readonly("pos_x", &Event::pos_x)
        .def_readonly("pos_y", &Event::pos_y)
        .def_readonly("pos_z", &Event::pos_z)
        .def_readonly("E_truth", &Event::E_truth)
        .def_readonly("edge_index", &Event::edge_index);
    m.def("read_root", &ReadRoot, py::arg("filename"), py::arg("threshold")=0.0, py::arg("KNN_K")=8);
}
