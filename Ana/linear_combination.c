// #include "/home/llr/ilc/ritzmann/work/root_macros/sum_hits_test.c"
#include "/home/llr/ilc/ritzmann/work/root_macros/sum_hits_weighted.c"
#include "/home/llr/ilc/ritzmann/work/root_macros/sum_energy_weighted2.c"


TFile* open_tree(float input_energy){
    // cout << "file with energy: " << input_energy << " openend!" << endl;
    string argument1, argument2, beam_energy;
    if(input_energy < 1){
        argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        argument2 = "keV.root"; 
        beam_energy = to_string(int(1000*input_energy));
    }

    else if((input_energy >= 1) && (input_energy <= 150)){
        argument1 = "/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_";
        argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        beam_energy = to_string(int(input_energy));
    }

    else {
        argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        argument2 = "GeV.root";
        beam_energy = to_string(int(input_energy));
    }

    string argument = argument1 + beam_energy + argument2;
    TFile *input = new TFile((char*)argument.c_str(), "read");
    return input;
}

TSpline3* spline_hits = sum_hits_weighted();
TSpline3* spline_energy = sum_energy_weighted2();

Double_t energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

// function to be minimized
double chi_squared(float energy, int entries, double energy_mean, double energy_val[5000], double hit_mean, double hit_val[5000], const double *alpha){
    double chi = 0;

    for (int i = 0; i < entries; i++){
        double nom = alpha[0]*(energy_val[i] - energy_mean) + (1. - alpha[0])*(hit_mean - hit_val[i]);
        double denom = 17/TMath::Sqrt(energy);
        chi += TMath::Power(nom/denom, 2);
    }
    return chi;
} 

void NumericalMinimization(double_t (&alpha)[16], double_t (&alpha_err)[16]){
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);

    double energy_mean[16];
    double std_energy[16];

    double hit_mean[16];
    double std_hit[16];

    double energy_val[5000];
    double hit_val[5000];

    for (int i = 0; i < 16; i++){
        cout << "energy: " << energies[i] << endl;

        TFile* input = open_tree(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");
        int entries = tree->GetEntries();


        TCanvas *c1 = new TCanvas(("hit_energy" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str());
        tree->Draw("Sum$(hit_energy*f2(hit_slab))");
        TH1F* hist_ene0 = (TH1F*)gPad->GetPrimitive("htemp");

        double xmin_ene = spline_energy->Eval(hist_ene0->GetXaxis()->GetXmin());
        double xmax_ene = spline_energy->Eval(hist_ene0->GetXaxis()->GetXmax());
        int nbins_ene;

        if (energies[i] < 0.5){
            nbins_ene = 30;
        }
        else if (energies[i] < 1){
            nbins_ene = 50;
        }
        else if (energies[i] < 5){
            nbins_ene = 60;
        }
        else if (energies[i] == 60 or energies[i] == 80){
            nbins_ene = 100;
        }
        else {
            nbins_ene = 80;
        }  

        TH1F* hist_hitene = new TH1F(("hit_energy" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str(), nbins_ene , xmin_ene, xmax_ene);

        
        TCanvas *c2 = new TCanvas(("wnhits" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str());
        tree->Draw("Sum$(f2(hit_slab))");
        TH1F* hist0 = (TH1F*)gPad->GetPrimitive("htemp");
    
        double xmin = hist0->GetXaxis()->GetXmin();
        double xmax = hist0->GetXaxis()->GetXmax();

        double factor;
        if (energies[i] <= 1){
            factor = 1.2;
        }
        else if (energies[i] <= 40 && energies[i] > 1){
        factor = 2;
        }
        else if (energies[i] > 40){
            factor = energies[i]/20;
        }

        int nbins = ceil((xmax-xmin)/factor);
        // cout << "nbins: " << nbins << endl;
        Double_t xbins[nbins + 1];
        xbins[0] = spline_hits->Eval(xmin);
            
        for (int i = 0; i < nbins+1; i++){
            double val = spline_hits->Eval(xmin+factor*(i+1));
            xbins[i+1] = val;
        }
        TH1F* hist_wnhits = new TH1F((to_string(energies[i])).c_str(), ("wnhits " + to_string(energies[i])).c_str(), nbins , xbins);

        double energy_val[5000];
        double hit_val[5000];

        int nhit_len = 0;
        tree->SetBranchAddress("nhit_len", &nhit_len);

        vector<int> *hit_slab = 0;
        tree->SetBranchAddress("hit_slab", &hit_slab);

        vector<int> *hit_energy = 0;
        tree->SetBranchAddress("hit_energy", &hit_energy);

        for (int j = 0; j < entries; j++){
            tree->GetEntry(j);
            int *slice = hit_slab->data();
            int *energy = hit_energy->data();
            double wnhits = 0;
            double whitene = 0;
            for (int k = 0; k < nhit_len; k++) {
                int layer = *slice;
                int ene = *energy;
                wnhits = wnhits + f2(layer);
                whitene = whitene + f2(layer)*ene;
                slice++;
                energy++;
            }
            hist_hitene->Fill(spline_energy->Eval(whitene));
            energy_val[j] = spline_energy->Eval(whitene);
            hist_wnhits->Fill(spline_hits->Eval(wnhits));
            hit_val[j] = spline_hits->Eval(wnhits);
        }
        hit_mean[i] = hist_wnhits->GetMean();
        std_hit[i] = hist_wnhits->GetStdDev();
        energy_mean[i] = hist_hitene->GetMean();
        std_energy[i] = hist_hitene->GetStdDev();

        std::function<double(const double*)> boundF = std::bind(chi_squared, energies[i], entries, energy_mean[i], energy_val, hit_mean[i], hit_val, std::placeholders::_1);

        ROOT::Math::Functor f(boundF, 1);
        double step[9] = {0.001};

        TRandom2 r(1);
        double variable[2] = {r.Uniform(-20,20)};

        min->SetFunction(f);

        min->SetVariable(0, "alpha", variable[0], step[0]);

        min->Minimize();

        const double *xs = min->X();
        const double *errors = min->Errors();

        alpha[i] = xs[0];
        alpha_err[i] = errors[0];
        cout << "minimum: chi_squ(" << xs[0] << ") = " << boundF(xs) << endl;
        input->Close();
    }
}

void linear_combination(){
    double_t alpha[16];
    double_t alpha_err[16];

    NumericalMinimization(alpha, alpha_err);

    ofstream linear_file;
    linear_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_linear_mean.txt");
    linear_file << "energy,mean,sigma" << "\n";

    TGraphErrors *alpha_ene = new TGraphErrors(16, energies, alpha, 0, alpha_err);
    TCanvas *alpha_val = new TCanvas("alpha_vals", "alpha_vals");
    alpha_ene->Draw();

    float xmin, xmax;
    int nbins;

    for (int i = 0; i < 16; i++){
        if (energies[i] < 1){
            xmin = 0; 
            xmax = 2;
            nbins = 15;
        }

        else if (energies[i] < 10){
            xmin = 0.4 * energies[i];
            xmax = 2 * energies[i];
            nbins = 50;
        }

        else if (energies[i] < 100){
            xmin = 0.5 * energies[i];
            xmax = 1.5 * energies[i];
            nbins = 80;
        }
        else {
            xmin = 0.7 * energies[i];
            xmax = 1.8 * energies[i];
            nbins = 80;
        }

        TCanvas *c_ene = new TCanvas(("ene_dist" + to_string(energies[i])).c_str(), ("ene_dist" + to_string(energies[i])).c_str());

        TH1F* hist_ene = new TH1F(("hist_ene" + to_string(energies[i])).c_str(), ("hist_ene" + to_string(energies[i])).c_str(), nbins , xmin, xmax);

        TFile* input = open_tree(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");
        int entries = tree->GetEntries();

        int nhit_len = 0;
        tree->SetBranchAddress("nhit_len", &nhit_len);

        vector<int> *hit_slab = 0;
        tree->SetBranchAddress("hit_slab", &hit_slab);

        vector<int> *hit_energy = 0;
        tree->SetBranchAddress("hit_energy", &hit_energy);

        for (int j = 0; j < entries; j++){
            tree->GetEntry(j);
            int *slice = hit_slab->data();
            int *energy = hit_energy->data();
            double wnhits = 0;
            double whitene = 0;
            for (int k = 0; k < nhit_len; k++) {
                int layer = *slice;
                int ene = *energy;
                wnhits = wnhits + f2(layer);
                whitene = whitene + f2(layer)*ene;
                slice++;
                energy++;
            }
            double hit_val = spline_hits->Eval(wnhits);
            double ene_val = spline_energy->Eval(whitene);
            hist_ene->Fill(alpha[i]*ene_val + (1. - alpha[i])*hit_val);
        }

        hist_ene->Draw();
        linear_file << energies[i] << "," << hist_ene->GetMean() << "," << hist_ene->GetStdDev() << "\n";
    }
}