// weighted sum of hits, where the weighting is according to the energy, divided into two areas
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
    // TTree *tree = (TTree*)input->Get("ecal");

    return input;
}


// double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3
float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

double energy_rec(TTree* tree, const double *par, int entry){

    vector<int> *hit_energy = 0;
    tree->SetBranchAddress("hit_energy", &hit_energy);

    // int entries = tree->GetEntries();

    int N1 = 0; 
    int N2 = 0;

    // for (int i = 0; i < entries; i++){
    tree->GetEntry(entry);
    int n_hit = hit_energy->size();
    int *ene = hit_energy->data();
    for (int j = 0; j < n_hit; j++){
        if (*ene >= 0.5 && *ene < 15){
            N1 = N1 + 1;
        }
        if (*ene >= 15){
            N2 = N2 + 1;
        }
        ene++;
    }
    int n_hits = N1 + N2;

    return par[0] * N1 + par[1] * N2;
}

double hits_weighted(TTree* tree, double alpha_par[3], double beta_par[3], int entry){

    vector<int> *hit_energy = 0;
    tree->SetBranchAddress("hit_energy", &hit_energy);

    int N1 = 0; 
    int N2 = 0;

    tree->GetEntry(entry);
    int n_hit = hit_energy->size();
    int *ene = hit_energy->data();
    for (int j = 0; j < n_hit; j++){
        if (*ene >= 0.5 && *ene < 3){
            N1 = N1 + 1;
        }
        if (*ene >= 3){
            N2 = N2 + 1;
        }
        ene++;
    }
    int n_hits = N1 + N2;
    double alpha = alpha_par[2] * TMath::Power(n_hits,2) + alpha_par[1] * n_hits + alpha_par[0];
    double beta = beta_par[2] * TMath::Power(n_hits,2) + beta_par[1] * n_hits + beta_par[0];

    return alpha * N1 + beta * N2;
}


double chi_squ(float energy, const double *par){
    double chi = 0; 
    TFile* input = open_tree(energy);
    TTree *tree = (TTree*)input->Get("ecal");
    int entries = tree->GetEntries();

    for (int i = 0; i < entries; i++){
        double ene_rec = energy_rec(tree, par, i);
        chi += TMath::Power((energy - ene_rec),2)/energy;
    }
    input->Close();
    return chi;
    }


void NumericalMinimization(double_t (&alpha)[16], double_t (&beta)[16], double_t (&alpha_par)[3], double_t (&beta_par)[3]){
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    // min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);

    double mean_nhits[16];
    double std_nhits[16];

    double alpha_errors[16];
    double beta_errors[16];

    for (int i = 0; i < 16; i++){
        cout << "energy: " << energies[i] << endl;

        // calculate the minimization for alpha and beta for one individual energy value
        std::function<double(const double*)> boundF = std::bind(chi_squ, energies[i], std::placeholders::_1);
        ROOT::Math::Functor f(boundF, 2);
        double step[9] = {0.01, 0.01};

        TRandom2 r(1);
        double variable[2] = {r.Uniform(-20,20), r.Uniform(-20,20)};

        min->SetFunction(f);

        min->SetVariable(0, "alpha", variable[0], step[0]);
        min->SetVariable(1, "beta", variable[1], step[1]);

        min->Minimize();

        const double *xs = min->X();
        const double *errors = min->Errors();

        cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << endl;
        alpha[i] = xs[0];
        beta[i] = xs[1];
        cout << "beta: " << beta[i] << endl;
        alpha_errors[i] = errors[0];
        beta_errors[i] = errors[1];

        // calculate the number of hits
        TFile* input = open_tree(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");
        int entries = tree->GetEntries();
        tree->Draw("nhit_len>>nhits");
        TH1F* h_nhits = (TH1F*)gPad->GetPrimitive("nhits");

        int xmax = h_nhits->GetXaxis()->GetXmax();

        int nbins = xmax;
        Double_t xbins[nbins+1];
        for (int i = 0; i < nbins+1; i++){
            xbins[i] = i + 0.5;
        }
        TH1F* hits = new TH1F("nhits", "number of hits", nbins, xbins);
        int nhit_len = 0;
        tree->SetBranchAddress("nhit_len", &nhit_len);

        for (int i = 0; i < entries; i++){
            tree->GetEntry(i);
            if (nhit_len!=0){
                hits->Fill(nhit_len);
            }
        }
        mean_nhits[i] = hits->GetMean();
        std_nhits[i] = hits->GetStdDev();
    }

    //parametrize alpha and beta as 2nd order poly of nhits
    TGraphErrors *alpha_hits = new TGraphErrors(16, mean_nhits, alpha, std_nhits, alpha_errors);
    TF1 *a_fit = (TF1 *)(gROOT->GetFunction("pol2"));
    alpha_hits->Fit(a_fit, "WQ");
    alpha_hits->Fit(a_fit, "Q","", mean_nhits[5]);
    TCanvas *alpha_fit = new TCanvas("alpha_polynomial", "alpha_polynomial");
    alpha_hits->Draw("APL*");
    TF1 *alpha_pol = alpha_hits->GetFunction("pol2");
    double a0 = alpha_pol->GetParameter(0);
    double a1 = alpha_pol->GetParameter(1);
    double a2 = alpha_pol->GetParameter(2);
    cout << "reduced Chi-square for alpha fit: " << alpha_pol->GetChisquare() / alpha_pol->GetNDF() << endl;
    alpha_par[0] = a0, alpha_par[1] = a1, alpha_par[2] = a2; 


    TGraphErrors *beta_hits = new TGraphErrors(16, mean_nhits, beta, std_nhits, beta_errors);
    TF1 *b_fit = (TF1 *)(gROOT->GetFunction("pol2"));
    beta_hits->Fit(b_fit, "WQ");
    beta_hits->Fit(b_fit, "Q", "", mean_nhits[5]);
    TCanvas *beta_fit = new TCanvas("beta_polynomial", "beta_polynomial");
    beta_hits->Draw("APL*");
    TF1 *beta_pol = beta_hits->GetFunction("pol2");
    double b0 = beta_pol->GetParameter(0);
    double b1 = beta_pol->GetParameter(1);
    double b2 = beta_pol->GetParameter(2);
    cout << "reduced Chi-square for beta fit: " << beta_pol->GetChisquare() / beta_pol->GetNDF() << endl;
    beta_par[0] = b0, beta_par[1] = b1, beta_par[2] = b2;
    }




void sum_hits_weighted_energy2(){
    // TTree *tree = open_tree(250);
    // cout << tree->GetEntries() << endl;

    double_t alpha[16];
    double_t beta[16];

    double_t alpha_par[3];
    double_t beta_par[3];

    double_t alpha_pol[3] = {0.0340894, 0.000194128, 0};
    double_t beta_pol[3] = {0.0411688, 0.000248029, 0};

    NumericalMinimization(alpha, beta, alpha_par, beta_par);

    // // for (int i = 0; i < 16; i++){
    // //     cout << energies[i] << " " << "alpha: " << alpha[i] << ", beta: " << beta[i] << endl;
    // // }
    // // const double par[9] = {-2.66907e-08, 0.000159691, 0.032614, -1.56747e-08, 0.000137639, 0.0365951, -1.14926e-08, 0.000124234, 0.0384119};

    // ofstream wnhits_file;
    // wnhits_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_weighted_energy_mean_test.txt");
    // wnhits_file << "energy,mean,sigma" << "\n";
    // int nbins, min_bin, max_bin;
    // for (int i = 0; i < 16; i++){
    //     TCanvas *c1 = new TCanvas(("weighted hits" + to_string(energies[i])).c_str());

    //     if (i < 5){
    //         nbins = 25;
    //         min_bin = 0;
    //         max_bin = 3 * energies[i];

    //     }

    //     if (i < 8 && i > 4){
    //         nbins = 50;
    //         min_bin = 0;
    //         max_bin = 3 * energies[i];

    //     }

    //     if (i >= 8){
    //         nbins = 70;
    //         min_bin = max(0.0, energies[i]*0.2);
    //         max_bin = energies[i]*1.8;
    //     }
    //     // TH1F *hist = new TH1F((to_string(energies[i])).c_str(), (to_string(energies[i])).c_str(), 50, max(0.0, energies[i]*0.2), energies[i]*1.8);
    // TH1F *hist = new TH1F((to_string(energies[i])).c_str(), (to_string(energies[i])).c_str(), nbins, min_bin, max_bin);
    // TFile* input = open_tree(energies[i]);
    // TTree *tree = (TTree*)input->Get("ecal");       
    // int entries = tree->GetEntries();
    // // tree->Draw("nhit_len>>nhits");
    // // TH1F* nhits = TH1F* hist_init = (TH1F*)gPad->GetPrimitive("nhits");
    // // int minbin = nhits->GetXaxis()->GetXmin()
    // // int maxbin = nhits->GetXaxis()->GetXmax()

    //     for (int j = 0; j < entries; j++){
    //         double ene_rec = hits_weighted(tree, alpha_par, beta_par, j);
    //         hist->Fill(ene_rec);
    //     }
    //     input->Close();
    //     hist->Draw();
    //     double_t mean = hist->GetMean();
    //     double_t sigma = hist->GetStdDev();
    //     // hist->Fit("gaus");
    //     // TF1 *fit = hist->GetFunction("gaus");
        
    //     // double_t mean = fit->GetParameter(1);
    //     // double sigma = fit->GetParameter(2);
    //     // double goodness = fit->GetChisquare() / fit->GetNDF();
    //     // cout << "goodness-of-fit:" << goodness << endl;
    //     // gStyle->SetOptFit(11111);
    //     // gSystem->ProcessEvents();
    //     // gSystem->ProcessEvents();
    //     // hist->Draw();

    //     wnhits_file << energies[i] << "," << mean << "," << sigma << "\n";

        
    // }
}


