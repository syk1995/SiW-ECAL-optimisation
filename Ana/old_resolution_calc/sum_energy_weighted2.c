
// double f1(int layer){
//     // X0 stands for the radiation length in the corresponding material

// 	double X0W = 3.504;
// 	double X0Si = 93.70*1000; //microns!
// 	double Ecal_WThickness = double(0.7)/X0W;
// 	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
// 	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
// 	// double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
//     double f0 = 1/(d_SiX0[layer] + d_WX0[layer]);
	
// 	return 1.0/f0;
// }
#include "/home/llr/ilc/ritzmann/work/root_macros/functions_resolution_calc.c"

void sum_energy_weighted_notcorr2(char* filepath, float input_energy, double_t &mean, double_t &error_val, ofstream &file_notcorr);
void sum_energy_weighted_corr2(char* filepath, TSpline3* spline_fit, float input_energy, ofstream &mean_file, ofstream &median_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms68_file, ofstream &rms90_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit, ofstream &gamma_file);
void calc_median__2(TTree *tree, TH1F *h, TSpline3* spline_fit, ofstream &median_file, float input_energy);

// weighted sum of the hit energies, weighting according to the tungsten thickness
TSpline3* sum_energy_weighted2(){

    // energy values of the simulation
    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};
    double_t energies_plot[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    // arrays to store the mean values and corresponding errors
    double_t mean_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t error_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t y_errors[16] = {10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7};

    // txt file for mean number of hits values
    ofstream file_notcorr;
    file_notcorr.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_mean_test2.txt");
    // file_notcorr << "energy,mean,mean_err,sigma,sigma_err, goodness" << "\n";
    file_notcorr << "energy,mean,sigma" << "\n";

    // read the root files and calculate the mean number of hits and corresponding error
    for(int i = 0; i < 16; i++) {
        string argument1, argument2, beam_energy;
        if(i < 3){
            argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
            argument2 = "keV.root"; 
            beam_energy = to_string(1000*energies[i]);
        }

        else if((i >= 3) && (i < 14)){
            argument1 = "/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_";
            argument2 = "GeV_5kevt_-42_-42_build_masked.root";
            beam_energy = to_string(energies[i]);
        }

        else {
            argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
            argument2 = "GeV.root";
            beam_energy = to_string(energies[i]);
        }
        // string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val, error_val;
        sum_energy_weighted_notcorr2((char*)argument.c_str(), energies[i], mean_val, error_val, file_notcorr);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }
    file_notcorr.close();

    // spline interpolation
    TSpline3 *spline_fit = new TSpline3("spline_interpol", mean_values, energies_plot, 16);

    // polynomial fit
    // TCanvas *c_fit = new TCanvas("polynomial_fit", "polynomial_fit");
    // TGraphErrors *hits = new TGraphErrors(9, mean_values, energies_plot, error_values, 0);
    // TF1 *gr_fit = ((TF1 *)(gROOT->GetFunction("pol4")));
    // hits->Fit(gr_fit, "WQ"); // "initial pre-fit"
    // hits->Fit(gr_fit, "Q"); // "final fit"   
    // hits->Draw("APL*");
    // gSystem->ProcessEvents();
    // gSystem->ProcessEvents();
    // TF1 *fit_pol = hits->GetFunction("pol4");
    // double a0 = fit_pol->GetParameter(0);
    // double a1 = fit_pol->GetParameter(1);
    // double a2 = fit_pol->GetParameter(2);
    // double a3 = fit_pol->GetParameter(3);
    // double a4 = fit_pol->GetParameter(4);
    // double fit_params[5] = {a0, a1, a2, a3, a4};


    // txt files, which store the predicted energy values for the different calculation methods
    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_mean_test2.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream median_file;
    median_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_median_test2.txt");
    median_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_gauss_test2.txt");
    gauss_file << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream gauss_file2s;
    gauss_file2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_gauss2s_test2.txt");
    gauss_file2s << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream rms68_file;
    rms68_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_rms68_test2.txt");
    rms68_file << "energy,mean,sigma" << "\n";

    ofstream rms90_file;
    rms90_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_rms90_test2.txt");
    rms90_file << "energy,mean,sigma" << "\n";

    ofstream myfile_3sig;
    myfile_3sig.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy_weighted_test2.txt");
    myfile_3sig << "energy,ratio" << "\n";

    ofstream myfile_3sig_fit;
    myfile_3sig_fit.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy_weighted_fitted_test2.txt");
    myfile_3sig_fit << "energy,ratio" << "\n";

    ofstream gamma_file;
    gamma_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_gamma_test2.txt");
    gamma_file << "energy,mean,sigma,chisq,ndf" << "\n";

    for(int i = 0; i < 16; i++) {
        string argument1, argument2, beam_energy;
        if(i < 3){
            argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
            argument2 = "keV.root"; 
            beam_energy = to_string(1000*energies[i]);
        }

        else if((i >= 3) && (i < 14)){
            argument1 = "/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_";
            argument2 = "GeV_5kevt_-42_-42_build_masked.root";
            beam_energy = to_string(energies[i]);
        }

        else {
            argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
            argument2 = "GeV.root";
            beam_energy = to_string(energies[i]);
        }
        // string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_weighted_corr2((char*)argument.c_str(), spline_fit, energies[i], mean_file, median_file, gauss_file, gauss_file2s, rms68_file,  rms90_file, myfile_3sig, myfile_3sig_fit, gamma_file);

    }
    mean_file.close();
    median_file.close();
    gauss_file.close();
    gauss_file2s.close();
    rms68_file.close();
    rms90_file.close();
    myfile_3sig.close();
    myfile_3sig_fit.close();
    gamma_file.close();

    return spline_fit;
}


// calculate the initial estimate, could be mean or also median
void sum_energy_weighted_notcorr2(char* filepath, float input_energy, double_t &mean, double_t &error_val, ofstream &file_notcorr){
   
    TFile *input = new TFile(filepath, "read");
    TTree *tree = (TTree*)input->Get("ecal");
    double_t entries = tree->GetEntries();

    TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

    tree->Draw("Sum$(hit_energy*f2(hit_slab))"); 
    TH1F* hist_init = (TH1F*)gPad->GetPrimitive("htemp");
    // TH1F* hist = (TH1F*)gPad->GetPrimitive("htemp");
    double xmin = hist_init->GetBinLowEdge(1);
    double xmax = hist_init->GetBinLowEdge(hist_init->GetNbinsX()+1);

    int nbins;
    if (input_energy < 0.5){
        nbins = 30;
    }
    else if (input_energy < 1){
        nbins = 50;
    }
    else if (input_energy < 5){
        nbins = 60;
    }
    else if (input_energy == 60 or input_energy == 80){
        nbins = 100;
    }
    else {
        nbins = 80;
    }  

    TH1F* hist = new TH1F(("hits" + to_string(input_energy)).c_str(), to_string(input_energy).c_str(), nbins , xmin, xmax);
    
    vector<int> *hit_energy = 0;
    vector<int> *hit_slab = 0;

    tree->SetBranchAddress("hit_energy", &hit_energy);
    tree->SetBranchAddress("hit_slab", &hit_slab);

    // iterating over all events
    for (int i = 0; i < entries; i++) {
        tree->GetEntry(i);
        int n_hit = hit_energy->size(); //number of hits per event
        int *energy = hit_energy->data(); //energy of each hit
        int *slice = hit_slab->data(); //layer number of each hit
        double tot_energy = 0;
        double sumE = 0;
        for (int j = 0; j < n_hit; j++) {
            int layer = *slice;
            tot_energy = tot_energy + *energy * f2(layer);
            sumE = sumE + *energy;
            energy++;
            slice++;
        }
        hist->Fill(tot_energy);
        // cout << (sumE == sum_energy) << endl;
     }

    hist->Draw();
    
    mean = hist->GetMean();
    double_t std_hist = hist->GetStdDev();
    error_val = std_hist / TMath::Sqrt(entries);

    gSystem->ProcessEvents();
    gSystem->ProcessEvents();

    file_notcorr << input_energy << "," << mean << "," << std_hist << "\n";
}

void sum_energy_weighted_corr2(char* filepath, TSpline3* spline_fit, float input_energy, ofstream &mean_file, ofstream &median_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms68_file, ofstream &rms90_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit, ofstream &gamma_file){
  
    TFile *input = new TFile(filepath, "read");
    TTree *tree = (TTree*)input->Get("ecal");
    double_t entries = tree->GetEntries(); 
    cout << "Number of events:" << entries << endl;

    TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

    tree->Draw("Sum$(hit_energy*f2(hit_slab))");
    TH1F* hist0 = (TH1F*)gPad->GetPrimitive("htemp");
    double_t mean0 = hist0->GetMean();
    double_t sigma0 = hist0->GetStdDev();

    double bin_width = 3.49*sigma0/TMath::Power(entries, 1/3);

    double xmin = spline_fit->Eval(hist0->GetBinLowEdge(1));
    double xmax = spline_fit->Eval(hist0->GetBinLowEdge(hist0->GetNbinsX()+1));

    int nbins;

    if (input_energy < 0.5){
        nbins = 30;
    }
    else if (input_energy < 1){
        nbins = 50;
    }
    else if (input_energy < 5){
        nbins = 60;
    }
    else if (input_energy == 60 or input_energy == 80){
        nbins = 100;
    }
    else {
        nbins = 80;
    }  

    TH1F* hist_corr = new TH1F(("corr_hits" + to_string(input_energy)).c_str(), to_string(input_energy).c_str(), nbins , xmin, xmax);
    
    float_t sum_energy = 0;
    tree->SetBranchAddress("sum_energy", &sum_energy);

    vector<int> *hit_energy = 0;
    vector<int> *hit_slab = 0;

    tree->SetBranchAddress("hit_energy", &hit_energy);
    tree->SetBranchAddress("hit_slab", &hit_slab);

    // iterating over all events
    for (int i = 0; i < entries; i++) {
        tree->GetEntry(i);
        int n_hit = hit_energy->size(); //number of hits per event
        int *energy = hit_energy->data(); //energy of each hit
        int *slice = hit_slab->data(); //layer number of each hit
        double tot_energy = 0;
        double sumE = 0;
        for (int j = 0; j < n_hit; j++) {
            int layer = *slice;
            tot_energy = tot_energy + *energy * f2(layer);
            sumE = sumE + *energy;
            energy++;
            slice++;
        }
        hist_corr->Fill(spline_fit->Eval(tot_energy));
        // cout << (sumE == sum_energy) << endl;
     }

    hist_corr->Draw();
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    
    TH1F* hist2_sig = (TH1F*)hist_corr->Clone();

    double ratio_3sig = calc_ratio(hist_corr, entries);
    myfile_3sig << input_energy << "," << ratio_3sig << "\n";

    double mean_hist = hist_corr->GetMean();
    double std_hist = hist_corr->GetStdDev();
    mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

    double mean_rms68, rms_val68;
    rms68(hist_corr, mean_rms68, rms_val68);
    // correction value for rms68
    rms68_file << input_energy << "," << mean_rms68 << "," << 1.769 * rms_val68 << "\n";

    double mean_rms90, rms_val90;
    rms90(hist_corr, mean_rms90, rms_val90);
    // correction value for rms90
    rms90_file << input_energy << "," << mean_rms90 << "," << 1.232 * rms_val90 << "\n";

    calc_median__2(tree, hist_corr, spline_fit, median_file, input_energy);
        
    fit_gauss(hist_corr, gauss_file, input_energy, myfile_3sig_fit, entries);
    fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);

    gamma_fit(hist_corr, input_energy, gamma_file);

    }

    // median and corresponding error for weighted sum energy 2
    void calc_median__2(TTree *tree, TH1F *h, TSpline3* spline_fit, ofstream &median_file, float input_energy){
        double median, MAD;
        double percentile = 0.5;
        h->ComputeIntegral();
        h->GetQuantiles(1, &median, &percentile);

        double_t entries = tree->GetEntries(); 
        cout << "Number of events:" << entries << endl;

        double xmin_mad = h->GetBinLowEdge(1) - median;
        double xmax_mad = h->GetBinLowEdge(h->GetNbinsX()+1) - median;

        TH1F* hist_mad = new TH1F(("hist_mad" + to_string(input_energy)).c_str(), ("hist_mad"+to_string(input_energy)).c_str(), 100 , xmin_mad, xmax_mad);
        
        vector<int> *hit_energy = 0;
        vector<int> *hit_slab = 0;
        tree->SetBranchAddress("hit_energy", &hit_energy);
        tree->SetBranchAddress("hit_slab", &hit_slab);

        // iterating over all events
        for (int i = 0; i < entries; i++) {
            tree->GetEntry(i);
            int n_hit = hit_energy->size(); //number of hits per event
            int *energy = hit_energy->data(); //energy of each hit
            int *slice = hit_slab->data(); //layer number of each hit
            double tot_energy = 0;
            for (int j = 0; j < n_hit; j++) {
                int layer = *slice;
                tot_energy = tot_energy + *energy * f2(*slice);
                energy++;
                slice++;
            }
            hist_mad->Fill(abs(spline_fit->Eval(tot_energy) - median));

        }
        hist_mad->ComputeIntegral();
        hist_mad->GetQuantiles(1, &MAD, &percentile);

        median_file << input_energy << "," << median << "," << 1.4826 * MAD << "\n";
    }


