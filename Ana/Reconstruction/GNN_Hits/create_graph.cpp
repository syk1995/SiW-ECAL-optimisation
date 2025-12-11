// create_graph_single.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <torch/torch.h>
#include <torch/script.h>
#include "TFile.h"
#include "TTree.h"
#include "nanoflann.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace py = pybind11;
using namespace nanoflann;

struct Event {
    std::vector<float> energy;
    std::vector<float> pos_x;
    std::vector<float> pos_y;
    std::vector<float> pos_z;
    std::vector<std::pair<int,int>> edge_index;
    float E_truth;
};

template <typename T>
struct PointCloud {
    struct Point { T x, y, z; };
    std::vector<Point> pts;
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
        if(dim==0) return pts[idx].x;
        else if(dim==1) return pts[idx].y;
        return pts[idx].z;
    }
    template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
};

std::vector<std::pair<int,int>> build_knn_graph(const Event &ev, int k=3) {
    PointCloud<float> cloud;
    for(size_t i=0;i<ev.pos_x.size();++i)
        cloud.pts.push_back({ev.pos_x[i], ev.pos_y[i], ev.pos_z[i]});
    typedef KDTreeSingleIndexAdaptor<
        L2_Simple_Adaptor<float, PointCloud<float>>,
        PointCloud<float>, 3> kd_tree_t;
    kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();

    std::vector<std::pair<int,int>> edges;
    std::vector<unsigned int> idx(k+1);
    std::vector<float> dist(k+1);
    for(size_t i=0;i<cloud.pts.size();++i){
        index.knnSearch(&cloud.pts[i].x, k+1, &idx[0], &dist[0]);
        for(int j=1;j<=k;++j) edges.emplace_back(i, idx[j]);
    }
    return edges;
}

void CreateGraph(const std::string &filename, const std::string &pt_file,
                 float threshold=0.0, float KNN_ratio=0.5, int KNN_max=16)
{
    TFile file(filename.c_str(), "READ");
    if(!file.IsOpen()) throw std::runtime_error("Cannot open ROOT file");

    TTree *tree = (TTree*)file.Get("events");
    if(!tree) throw std::runtime_error("Cannot find TTree 'events'");

    std::vector<float> *energy=nullptr, *pos_x=nullptr, *pos_y=nullptr, *pos_z=nullptr, *E_truth=nullptr;
    tree->SetBranchAddress("simplecaloRO.energy", &energy);
    tree->SetBranchAddress("simplecaloRO.position.x", &pos_x);
    tree->SetBranchAddress("simplecaloRO.position.y", &pos_y);
    tree->SetBranchAddress("simplecaloRO.position.z", &pos_z);
    tree->SetBranchAddress("MCParticles.p0", &E_truth);

    Long64_t nentries = tree->GetEntries();

    std::vector<torch::Tensor> x_list, edge_list, y_list;

    for(Long64_t i=0;i<nentries;++i){
        if(i%1000==0)
            std::cout << "Processing event " << i << "/" << nentries << std::endl;
        tree->GetEntry(i);
        if(E_truth->size()!=1) continue;

        Event ev;
        ev.E_truth = (*E_truth)[0];
        for(size_t j=0;j<energy->size();++j){
            if((*energy)[j]*1000>threshold){
                ev.energy.push_back((*energy)[j]);
                ev.pos_x.push_back((*pos_x)[j]);
                ev.pos_y.push_back((*pos_y)[j]);
                ev.pos_z.push_back((*pos_z)[j]);
            }
        }
        if(ev.energy.empty()) continue;

        float zmin = *std::min_element(ev.pos_z.begin(), ev.pos_z.end());
        for(auto &z: ev.pos_z) z -= zmin;
        for(size_t j=0;j<ev.pos_z.size();++j)
            ev.pos_z[j] = std::round(ev.pos_z[j]*1000.0)/1000.0;
        ev.E_truth = std::round(ev.E_truth*1000.0)/1000.0;
        int KNN_K = std::min(KNN_max, std::max(1, (int)(ev.energy.size()*KNN_ratio)));
        ev.edge_index = build_knn_graph(ev, KNN_K);

        torch::Tensor x = torch::zeros({(int64_t)ev.energy.size(),4});
        for(size_t j=0;j<ev.energy.size();++j){
            x[j][0] = ev.energy[j];
            x[j][1] = ev.pos_x[j];
            x[j][2] = ev.pos_y[j];
            x[j][3] = ev.pos_z[j];
        }

        torch::Tensor edge = torch::zeros({2,(int64_t)ev.edge_index.size()}, torch::kLong);
        for(size_t j=0;j<ev.edge_index.size();++j){
            edge[0][j] = ev.edge_index[j].first;
            edge[1][j] = ev.edge_index[j].second;
        }

        torch::Tensor y = torch::tensor({ev.E_truth});

        x_list.push_back(x);
        edge_list.push_back(edge);
        y_list.push_back(y);
    }

    // 保存为一个pt文件，字典形式
    torch::serialize::OutputArchive archive;
    archive.write("x_list", torch::stack(x_list));
    archive.write("edge_list", torch::stack(edge_list));
    archive.write("y_list", torch::stack(y_list));
    archive.save_to(pt_file);

    file.Close();
}

PYBIND11_MODULE(create_graph, m){
    m.def("create_graph", &CreateGraph,
          py::arg("filename"), py::arg("pt_file"),
          py::arg("threshold")=0.0, py::arg("KNN_ratio")=0.5,py::arg("KNN_max")=16);
}
