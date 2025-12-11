// ReadRootClusters.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "nanoflann.hpp"
#include "iostream"
namespace py = pybind11;

// ---------------- Event Struct -----------------
struct Event {
    std::vector<float> energy;
    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<float> pos_z;
    float E_truth;
};
struct Event_Clusters {
    std::vector<std::vector<float>> features;
    float E_truth;
};

// ---------------- KDTree PointCloud -----------------
template <typename T>
struct PointCloud {
    struct Point { T x, y, z; };
    std::vector<Point> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].x;
        else if (dim == 1) return pts[idx].y;
        else return pts[idx].z;
    }
    template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

using namespace nanoflann;

// ---------------- Compute Energy Weighted Barycenter -----------------
std::tuple<float,float,float> compute_energy_barycenter(const Event &ev) {
    float totalE = std::accumulate(ev.energy.begin(), ev.energy.end(), 0.0f);
    float cx = 0, cy = 0, cz = 0;
    for(size_t i=0;i<ev.energy.size();++i) {
        cx += ev.energy[i] * ev.pos_x[i];
        cy += ev.energy[i] * ev.pos_y[i];
        cz += ev.energy[i] * ev.pos_z[i];
    }
    cx /= totalE; cy /= totalE; cz /= totalE;
    return {cx, cy, cz};
}

// ---------------- Find First Hit by Z + Energy Weighted for Multiple Minimum -----------------
std::tuple<float,float,float> compute_first_point(const Event &ev) {
    float min_z = *std::min_element(ev.pos_z.begin(), ev.pos_z.end());
    float totalE = 0, cx=0, cy=0, cz=0;
    for(size_t i=0;i<ev.energy.size();++i){
        if(ev.pos_z[i] == min_z){
            totalE += ev.energy[i];
            cx += ev.energy[i] * ev.pos_x[i];
            cy += ev.energy[i] * ev.pos_y[i];
            cz += ev.energy[i] * ev.pos_z[i];
        }
    }
    cx/=totalE; cy/=totalE; cz/=totalE;
    return {cx, cy, cz};
}

void sort_hits_by_particle_direction(Event &ev)
{
    auto [sx, sy, sz] = compute_first_point(ev);
    auto [cx, cy, cz] = compute_energy_barycenter(ev);
    float vx = cx - sx, vy = cy - sy, vz = cz - sz;
    if(std::abs(vz) < 1e-6f) vz = (vz >= 0 ? 1e-6f : -1e-6f);
    size_t N = ev.energy.size();
    std::vector<float> newx(N), newy(N), newz(N);
    std::vector<float> r2(N);

    for(size_t i=0;i<N;i++){
        float dz = ev.pos_z[i] - sz;
        float t = dz / vz;
        float px = sx + t*vx;
        float py = sy + t*vy;
        newx[i] = ev.pos_x[i]-px;
        newy[i] = ev.pos_y[i]-py;
        newz[i] = dz;
        r2[i] = newx[i]*newx[i] + newy[i]*newy[i];
    }
    std::vector<size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
        if(newz[a] != newz[b])
            return newz[a] < newz[b];
        return r2[a] < r2[b];
    });

    std::vector<float> x_out,y_out,z_out,e_out;
    x_out.reserve(N); y_out.reserve(N); z_out.reserve(N); e_out.reserve(N);

    for(size_t i: idx){
        e_out.push_back(ev.energy[i]);
        x_out.push_back(newx[i]);
        y_out.push_back(newy[i]);
        z_out.push_back(newz[i]);
    }

    ev.energy = std::move(e_out);
    ev.pos_x = std::move(x_out);
    ev.pos_y = std::move(y_out);
    ev.pos_z = std::move(z_out);
}

Event_Clusters build_cluster(Event &ev, int cluster_size)
{
    sort_hits_by_particle_direction(ev);
    PointCloud<float> cloud;
    cloud.pts.reserve(ev.pos_x.size());
    for(size_t i=0;i<ev.pos_x.size();++i)
        cloud.pts.push_back({ev.pos_x[i], ev.pos_y[i], ev.pos_z[i]});

    using kd_tree_t = KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<float, PointCloud<float>>,
        PointCloud<float>, 3
    >;

    kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::vector<std::vector<int>> clusters_idx(ev.pos_x.size());
    std::vector<unsigned int> idx(cluster_size);
    std::vector<float> dist(cluster_size);

    for(size_t i=0;i<cloud.pts.size();++i){
        size_t n = index.knnSearch(&cloud.pts[i].x, cluster_size, &idx[0], &dist[0]);
        for(size_t j=0;j<n;++j)
            clusters_idx[i].push_back(idx[j]);
        std::sort(clusters_idx[i].begin(), clusters_idx[i].end());
    }
    std::vector<std::vector<int>> unique_idx;
    std::set<std::vector<int>> seen;
    for (auto& ids : clusters_idx) {
        if (seen.insert(ids).second) {
            unique_idx.push_back(ids); // 保留插入顺序
        }
    }
    Event_Clusters result;
    result.E_truth = ev.E_truth;
    int count =0;
    for(auto const& ids : unique_idx){
        std::vector<float> feat;
        feat.reserve(cluster_size * 4 + 6);// 4 features per hit + 6 summary features: total E, N hits, E/N,cx, cy, cz
        float total_energy = 0.0f;
        float cx=0.0f, cy=0.0f, cz=0.0f;
        int nonzero_hits = ids.size();
        for(int id : ids){
            float e = ev.energy[id];
            feat.push_back(ev.pos_x[id]);
            feat.push_back(ev.pos_y[id]);
            feat.push_back(ev.pos_z[id]);
            feat.push_back(e);
            total_energy += e;
            cx += ev.pos_x[id] * e;
            cy += ev.pos_y[id] * e;
            cz += ev.pos_z[id] * e;
        }
        while(feat.size() < (cluster_size * 4)){
            feat.push_back(0.0f);
        }
        if (nonzero_hits < 1) continue;
        cx /= total_energy;
        cy /= total_energy;
        cz /= total_energy;
        feat.push_back(total_energy);
        feat.push_back(static_cast<float>(nonzero_hits));
        feat.push_back(total_energy / nonzero_hits);
        feat.push_back(cx);
        feat.push_back(cy);
        feat.push_back(cz);
        if(feat.size() != cluster_size *4 +6)
            throw std::runtime_error("Feature size " + std::to_string(feat.size()) + " mismatch");
        result.features.push_back(std::move(feat));
        count++;
    }
    return result;
}

std::vector<Event_Clusters> ReadRootClusters(const std::string &filename,
                                             float threshold,
                                             int cluster_size,
                                             int max_events=-1)
{
    std::vector<Event_Clusters> all_event_clusters;
    TFile file(filename.c_str(), "READ");
    if (!file.IsOpen()) throw std::runtime_error("Cannot open ROOT file");

    TTree *tree = (TTree*)file.Get("events");
    if (!tree) throw std::runtime_error("Cannot find TTree 'events'");

    std::vector<float> *energy=nullptr, *pos_x=nullptr, *pos_y=nullptr, *pos_z=nullptr, *E_truth=nullptr;
    tree->SetBranchAddress("simplecaloRO.energy", &energy);
    tree->SetBranchAddress("simplecaloRO.position.x", &pos_x);
    tree->SetBranchAddress("simplecaloRO.position.y", &pos_y);
    tree->SetBranchAddress("simplecaloRO.position.z", &pos_z);
    tree->SetBranchAddress("MCParticles.p0", &E_truth);

    if(max_events == -1) max_events = tree->GetEntries();
    std::cout << "Reading ROOT file: " << filename
              << " max_events=" << max_events
              << " threshold=" << threshold << std::endl;
    for(Long64_t i=0; i<max_events; ++i){
        tree->GetEntry(i);
        if(E_truth->size() != 1)
            throw std::runtime_error("E_truth size !=1 at event "+std::to_string(i));

        float truth = std::round((*E_truth)[0] * 1000.0f) / 1000.0f;

        std::vector<float> e,x,y,z;
        for(size_t j=0;j<energy->size();++j){
            if((*energy)[j]*1000 > threshold){
                e.push_back((*energy)[j]);
                x.push_back((*pos_x)[j]);
                y.push_back((*pos_y)[j]);
                z.push_back((*pos_z)[j]);
            }
        }
        if(e.empty()) continue;
        Event ev;
        ev.energy = std::move(e);
        ev.pos_x  = std::move(x);
        ev.pos_y  = std::move(y);
        ev.pos_z  = std::move(z);
        ev.E_truth = truth;

        Event_Clusters clusters = build_cluster(ev, cluster_size);
        all_event_clusters.push_back(std::move(clusters));
    }

    file.Close();
    return all_event_clusters;
}

PYBIND11_MODULE(ReadRootClusters, m) {

    py::class_<Event_Clusters>(m, "Event_Clusters")
        .def_readonly("features", &Event_Clusters::features)
        .def_readonly("E_truth", &Event_Clusters::E_truth);

    m.def("read_root_clusters", &ReadRootClusters,
          py::arg("filename"),
          py::arg("threshold") = 0.0,
          py::arg("cluster_size") = 4,
          py::arg("max_events") = -1);
}

