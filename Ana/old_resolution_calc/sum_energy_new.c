
void sum_energy_notcorr(char* filename, float input_energy, double_t &mean, double_t &error_val){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_mean.txt", std::ios_base::app);

    // double_t mean;
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");

if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");
        double_t entries = tree->GetEntries();


        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);


        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }

    else if ((input_energy >= 1) && (input_energy < 200)){   
        int int_input_energy = (int)input_energy; 
        cout << "small " << int_input_energy << endl;
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");
        double_t entries = tree->GetEntries();

        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);


        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }

    else {
        cout << "big" << endl;
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 

        TTree *tree = (TTree*)input->Get("ecal");
        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());
        double_t entries = tree->GetEntries();

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);


        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }
    // return mean;
}


double b0 = 0;
double b1 = 0;
double b2 = 0;
double b3 = 0;
double b4 = 0;

double eesum(float energy_sum){
    
    double val = b4*pow(energy_sum,4)+b3*pow(energy_sum,3)+b2*pow(energy_sum,2)+b1*energy_sum+b0;
    return val;
}

double_t sum_energy_notcorr(char* filename, int input_energy, double_t &mean, double_t &error_val);

TSpline3 *spline_fit_e;

double_t spline_hits_e(double energy_sum){
    double_t val = spline_fit_e->Eval(energy_sum);
    return val;
}



void sum_energy_corr(char* filename, float input_energy, ofstream &mean_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit){

    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        //get the number of events
        int entries = tree->GetEntries();
        cout << "Number of events:" << entries << endl;

        TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());
        // tree->Draw("eesum(Sum$(hit_energy))"); 
        tree->Draw("spline_hits_e(Sum$(hit_energy))");

        TH1F* hist_gauss = (TH1F*)gPad->GetPrimitive("htemp");

        double ratio_3sig = calc_ratio(hist_gauss, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";

        double mean_hist = hist_gauss->GetMean();
        double std_hist = hist_gauss->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        // rms68(hist_gauss, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";

        rms90(hist_gauss, mean_rms, rms_val);
        rms_file << input_energy << "," << mean_rms << "," << 1.23172 * rms_val << "\n";

        TH1F* hist_gauss_sig = (TH1F*)hist_gauss->Clone();
        
        fit_gauss(hist_gauss, gauss_file, input_energy, myfile_3sig_fit, entries);
        fit_gauss_2sig(hist_gauss_sig, gauss_file2s, mean_hist, std_hist, input_energy);
    }
    
    
    
    else if ((input_energy >= 1) && (input_energy < 200)){  
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 
        
        TTree *tree = (TTree*)input->Get("ecal");

        //get the number of events
        int entries = tree->GetEntries();
        cout << "Number of events:" << entries << endl;

        TCanvas *c = new TCanvas(("corr_c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        // tree->Draw("eesum(Sum$(hit_energy))"); 
        tree->Draw("spline_hits_e(Sum$(hit_energy))");
        TH1F* hist2 = (TH1F*)gPad->GetPrimitive("htemp");

        double ratio_3sig = calc_ratio(hist2, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";

        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        // rms68(hist2, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";

        rms90(hist2, mean_rms, rms_val);
        rms_file << input_energy << "," << mean_rms << "," << 1.23172 * rms_val << "\n";

        TH1F* hist2_sig = (TH1F*)hist2->Clone();

        fit_gauss(hist2, gauss_file, input_energy, myfile_3sig_fit, entries);
        fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);
    }

    else {
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 


        TTree *tree = (TTree*)input->Get("ecal");

        int entries = tree->GetEntries();


        TCanvas *c = new TCanvas(("corr_c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        // tree->Draw("eesum(Sum$(hit_energy))"); 
        tree->Draw("spline_hits_e(Sum$(hit_energy))");
        TH1F* hist2 = (TH1F*)gPad->GetPrimitive("htemp");

        double ratio_3sig = calc_ratio(hist2, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";

        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        // rms68(hist2, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";

        rms90(hist2, mean_rms, rms_val);
        rms_file << input_energy << "," << mean_rms << "," << 1.23172 * rms_val << "\n";

        TH1F* hist2_sig = (TH1F*)hist2->Clone();

        fit_gauss(hist2, gauss_file, input_energy, myfile_3sig_fit, entries);
        fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);
    }
}


void sum_energy_new(){

    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t energies_plot[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    // double_t energies_plot[13] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t mean_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t error_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t y_errors[16] = {10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7};



    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_mean.txt");
    myfile_events << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events.close();

    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val, error_val;
        sum_energy_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }
    
    for (int i = 3; i < 14; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val, error_val;
        sum_energy_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }

    for (int i = 14; i < 16; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val, error_val;
        sum_energy_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }

    TCanvas *c_fit = new TCanvas("gaussian_fit", "gaussian_fit");

    // TGraph *hits = new TGraph(13, mean_values, energies_plot);
    // hits->Fit("pol4");
    spline_fit_e = new TSpline3("spline_interpol", mean_values, energies_plot, 16);

    TGraphErrors *hits = new TGraphErrors(16, mean_values, energies_plot, error_values, 0);
    TF1 *gr_fit = ((TF1 *)(gROOT->GetFunction("pol4")));
    hits->Fit(gr_fit, "WQ"); // "initial pre-fit"
    hits->Fit(gr_fit, "Q"); // "final fit"
    
    hits->Draw("APL*");
    // hits->GetXaxis()->SetLimits(0,400);
    TF1 *fit_pol = hits->GetFunction("pol4");
    b0 = fit_pol->GetParameter(0);
    b1 = fit_pol->GetParameter(1);
    b2 = fit_pol->GetParameter(2);
    b3 = fit_pol->GetParameter(3);
    b4 = fit_pol->GetParameter(4);


    for (const auto& e : mean_values) {
        std::cout << e << " " << eesum(e) << std::endl;
        }

    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr_mean.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr_gauss.txt");
    gauss_file << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream gauss_file2s;
    gauss_file2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr_gauss2s.txt");
    gauss_file2s << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream rms_file;
    rms_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr_rms90.txt");
    rms_file << "energy,mean,sigma" << "\n";

    ofstream myfile_3sig;
    myfile_3sig.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy.txt");
    myfile_3sig << "energy,ratio" << "\n";

    ofstream myfile_3sig_fit;
    myfile_3sig_fit.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy_fitted.txt");
    myfile_3sig_fit << "energy,ratio" << "\n";
    
    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }
    
    
    for (int i = 3; i < 14; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }

    for (int i = 14; i < 16; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }  
    mean_file.close();
    gauss_file.close();
    gauss_file2s.close();
    rms_file.close();
    myfile_3sig.close();
    myfile_3sig_fit.close();
}

