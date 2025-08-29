// weighted sum of hits, where the weighting is according to the energy, divided into three areas
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


float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

double energy_rec(TTree* tree, const double *par, int entry){

    vector<int> *hit_energy = 0;
    tree->SetBranchAddress("hit_energy", &hit_energy);

    // int entries = tree->GetEntries();

    int N1 = 0; 
    int N2 = 0;
    int N3 = 0; 

    // for (int i = 0; i < entries; i++){
    tree->GetEntry(entry);
    int n_hit = hit_energy->size();
    int *ene = hit_energy->data();
    for (int j = 0; j < n_hit; j++){
        if (*ene >= 0.5 && *ene < 5){
            N1 = N1 + 1;
        }
        if (*ene >= 5 && *ene < 15){
            N2 = N2 + 1;
        }
        if (*ene >= 15){
            N3 = N3 + 1;
        }
    }
    int n_hits = N1 + N2 + N3;

    // double alpha = a1 * TMath::Power(n_hits,2) + a2*n_hits+a3;
    // double beta = b1 * TMath::Power(n_hits,2) + b2 * n_hits + b3;
    // double gamma = c1 * TMath::Power(n_hits,2) + c2 * n_hits + c3;  

    double alpha = par[0] * TMath::Power(n_hits,2.) + par[1]*n_hits+par[2];
    double beta = par[3] * TMath::Power(n_hits,2.) + par[4] * n_hits + par[5];
    double gamma = par[6] * TMath::Power(n_hits,2.) + par[7] * n_hits + par[8];   

    return alpha * N1 + beta * N2 + gamma * N3;
}

double chi_squ(const double *par){
    double chi = 0; 

    for (int i = 0; i < 16; i++){
        TFile* input = open_tree(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");       
        int entries = tree->GetEntries();
        for (int j = 0; j < entries; j++){
            double ene_rec = energy_rec(tree, par, j);
            chi += TMath::Power((energies[i] - ene_rec),2)/energies[i];
        }
        input->Close();
    }
    return chi;
}

void NumericalMinimization(){
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");
    min->SetMaxFunctionCalls(100000); // for Minuit/Minuit2
    // min->SetMaxIterations(10000);  // for GSL
    min->SetTolerance(0.001);
    min->SetPrintLevel(1);


    ROOT::Math::Functor f(&chi_squ, 9);
    double step[9] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

    TRandom2 r(1);
    double variable[9] = {r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20), r.Uniform(-20,20)};
    // double variable[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}

    min->SetFunction(f);

    min->SetVariable(0, "a2", variable[0], step[0]);
    min->SetVariable(1, "a1", variable[1], step[1]);
    min->SetVariable(2, "a0", variable[2], step[2]);
    min->SetVariable(3, "b2", variable[3], step[3]);
    min->SetVariable(4, "b1", variable[4], step[4]);
    min->SetVariable(5, "b0", variable[5], step[5]);
    min->SetVariable(6, "c2", variable[6], step[6]);
    min->SetVariable(7, "c1", variable[7], step[7]);
    min->SetVariable(8, "c0", variable[8], step[8]);


    min->Minimize();

    const double *xs = min->X();


    cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "," << xs[3] << "," << xs[4] << "," << xs[5] << "," << xs[6] << "," << xs[7] << "," << xs[8] << "): " << min->MinValue()  << endl;
   
}




void sum_hits_weighted_energy(){
    // TTree *tree = open_tree(250);
    // cout << tree->GetEntries() << endl;

    NumericalMinimization();
    const double par[9] = {-2.66907e-08, 0.000159691, 0.032614, -1.56747e-08, 0.000137639, 0.0365951, -1.14926e-08, 0.000124234, 0.0384119};

    ofstream wnhits_file;
    wnhits_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_weighted_energy_gauss_test.txt");
    wnhits_file << "energy,mean,sigma" << "\n";
    int nbins, min_bin, max_bin;
    for (int i = 0; i < 16; i++){
        TCanvas *c1 = new TCanvas(("weighted hits" + to_string(energies[i])).c_str());

        if (i < 5){
            nbins = 25;
            min_bin = 0;
            max_bin = 3 * energies[i];

        }

        if (i < 8 && i > 4){
            nbins = 50;
            min_bin = 0;
            max_bin = 3 * energies[i];

        }

        if (i >= 8){
            nbins = 70;
            min_bin = max(0.0, energies[i]*0.2);
            max_bin = energies[i]*1.8;
        }
        // TH1F *hist = new TH1F((to_string(energies[i])).c_str(), (to_string(energies[i])).c_str(), 50, max(0.0, energies[i]*0.2), energies[i]*1.8);
        TH1F *hist = new TH1F((to_string(energies[i])).c_str(), (to_string(energies[i])).c_str(), nbins, min_bin, max_bin);
        TFile* input = open_tree(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");       
        int entries = tree->GetEntries();
        for (int j = 0; j < entries; j++){
            double ene_rec = energy_rec(tree, par, j);
            hist->Fill(ene_rec);
        }
        input->Close();
        hist->Fit("gaus");
        TF1 *fit = hist->GetFunction("gaus");
        // double_t mean = hist->GetMean();
        // double_t sigma = hist->GetStdDev();
        double_t mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        double goodness = fit->GetChisquare() / fit->GetNDF();
        cout << "goodness-of-fit:" << goodness << endl;
        gStyle->SetOptFit(11111);
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();
        hist->Draw();

        wnhits_file << energies[i] << "," << mean << "," << sigma << "\n";

        
    }
}


