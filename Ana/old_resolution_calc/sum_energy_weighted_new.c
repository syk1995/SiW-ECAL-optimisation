
double f(int layer){
    // X0 stands for the radiation length in the corresponding material

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
    // double f0 = 1/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}

void sum_energy_weighted_notcorr(char* filename, float input_energy, double_t &mean, double_t &error_val){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_mean.txt", std::ios_base::app);

    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");

    // double_t mean; //mean value of number of hits

    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");
        TTree *tree = (TTree*)input->Get("ecal");
        double_t entries = tree->GetEntries();


        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("Sum$(hit_energy*f(hit_slab))"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);
        cout << mean << " " << error_val << endl;

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


        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("Sum$(hit_energy*f(hit_slab))"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);
        cout << mean << " " << error_val << endl;

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
        double_t entries = tree->GetEntries();

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("Sum$(hit_energy*f(hit_slab))"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double_t std_hist = hist->GetStdDev();
        // double_t bins = hist->GetNbinsX();
        error_val = std_hist / TMath::Sqrt(entries);
        cout << mean << " " << error_val << endl;

        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }
    // return mean;
}


double c0 = 0;
double c1 = 0;
double c2 = 0;
double c3 = 0;
double c4 = 0;

double ewesum(float energy_sum){
    
    double val = c4*pow(energy_sum,4)+c3*pow(energy_sum,3)+c2*pow(energy_sum,2)+c1*energy_sum+c0;
    return val;
}

TSpline3 *spline_fit_we1;

double_t spline_hits_we1(double energy_sum){
    double_t val = spline_fit_we1->Eval(energy_sum);
    return val;
}

double_t sum_energy_weighted_notcorr(char* filename, int input_energy, double_t &mean, double_t &error_val);

// void rms68(TH1F *h, double &mean, double &rms_val);
// void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy);
// void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy);




void sum_energy_weighted_corr(char* filename, float input_energy, ofstream &mean_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit){
    
    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        //get the number of events
        int entries = tree->GetEntries();
        cout << "Number of events:" << entries << endl;

        TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());
        // tree->Draw("ewesum(Sum$(hit_energy*f(hit_slab)))"); 
        tree->Draw("spline_hits_we1(Sum$(hit_energy*f(hit_slab)))");
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

        // tree->Draw("ewesum(Sum$(hit_energy*f(hit_slab)))"); 
        tree->Draw("spline_hits_we1(Sum$(hit_energy*f(hit_slab)))");

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

        // tree->Draw("ewesum(Sum$(hit_energy*f(hit_slab)))"); 
        tree->Draw("spline_hits_we1(Sum$(hit_energy*f(hit_slab)))");
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


void sum_energy_weighted_new(){

    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t energies_plot[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    double_t mean_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t error_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t y_errors[16] = {10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7};


    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_mean.txt");
    myfile_events << "energy,mean,sigma" << "\n";
    myfile_events.close();

    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;
        cout << argument << endl;

        double_t mean_val, error_val;
        sum_energy_weighted_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
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
        sum_energy_weighted_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
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
        sum_energy_weighted_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }

    TCanvas *c_fit = new TCanvas("gaussian_fit", "gaussian_fit");

    spline_fit_we1 = new TSpline3("spline_interpol", mean_values, energies_plot, 16);


    TGraphErrors *hits = new TGraphErrors(16, mean_values, energies_plot, error_values, 0);
    TF1 *gr_fit = ((TF1 *)(gROOT->GetFunction("pol4")));
    hits->Fit(gr_fit, "WQ"); // "initial pre-fit"
    hits->Fit(gr_fit, "Q"); // "final fit"
    
    // TGraph *hits = new TGraph(13, mean_values, energies_plot);
    // hits->Fit("pol4");

    hits->Draw("APL*");
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    // hits->GetXaxis()->SetLimits(0,400);
    TF1 *fit_pol = hits->GetFunction("pol4");
    c0 = fit_pol->GetParameter(0);
    c1 = fit_pol->GetParameter(1);
    c2 = fit_pol->GetParameter(2);
    c3 = fit_pol->GetParameter(3);
    c4 = fit_pol->GetParameter(4);

    //cout << p0 << p1<< p2 << p3 << endl;

    for (const auto& e : mean_values) {
    std::cout << e << " " << ewesum(e) << std::endl;
    }

    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_mean.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_gauss.txt");
    gauss_file << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream gauss_file2s;
    gauss_file2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_gauss2s.txt");
    gauss_file2s << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream rms_file;
    rms_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_corr_rms90.txt");
    rms_file << "energy,mean,sigma" << "\n";

    ofstream myfile_3sig;
    myfile_3sig.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy_weighted1.txt");
    myfile_3sig << "energy,ratio" << "\n";

    ofstream myfile_3sig_fit;
    myfile_3sig_fit.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_energy_weighted1_fitted.txt");
    myfile_3sig_fit << "energy,ratio" << "\n";
    
    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_weighted_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }
    
    
    for (int i = 3; i < 14; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_weighted_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }

    for (int i = 14; i < 16; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_weighted_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, myfile_3sig, myfile_3sig_fit);
    }  
    mean_file.close();
    gauss_file.close();
    gauss_file2s.close();
    rms_file.close();
    myfile_3sig.close();
    myfile_3sig_fit.close();
    
    
}

// void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy){
//     h->Fit("gaus");
//     TF1 *fit = h->GetFunction("gaus");
//     double_t mean = fit->GetParameter(1);
//     double sigma = fit->GetParameter(2);
//     double mean_err = fit->GetParError(1);
//     double sigma_err = fit->GetParError(2);
//     float res = sigma/mean;
//     double goodness = fit->GetChisquare() / fit->GetNDF();
//     cout << "resolution:" << res << endl;
//     cout << "goodness-of-fit:" << goodness << endl;
//     h->Draw();
//     gStyle->SetOptFit(11111);
//     gSystem->ProcessEvents();
//     gSystem->ProcessEvents();
//     gauss_file << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
// }

// void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy){
//     h->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
//     TF1 *fit = h->GetFunction("gaus");
//     double_t mean = fit->GetParameter(1);
//     double sigma = fit->GetParameter(2);
//     double mean_err = fit->GetParError(1);
//     double sigma_err = fit->GetParError(2);
//     float res = sigma/mean;
//     double goodness = fit->GetChisquare() / fit->GetNDF();
//     cout << "resolution:" << res << endl;
//     cout << "goodness-of-fit:" << goodness << endl;
//     h->Draw();
//     gStyle->SetOptFit(11111);
//     gSystem->ProcessEvents();
//     gSystem->ProcessEvents();
//     gauss_file2s << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
// }

// void rms68(TH1F *h, double &mean, double &rms_val) {
//     TAxis *axis = h->GetXaxis();
//     Int_t nbins = axis->GetNbins();
//     Int_t imean = axis->FindBin(h->GetMean());
//     Double_t entries = 0.68*h->GetEntries();
//     Double_t w = h->GetBinContent(imean);
//     Double_t x = h->GetBinCenter(imean);
//     Double_t sumw = w;
//     Double_t sumwx = w*x;
//     Double_t sumwx2 = w*x*x;
//     for (Int_t i=1;i<nbins;i++) {
//     if (i> 0) {
//         w = h->GetBinContent(imean-i);
//         x = h->GetBinCenter(imean-i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (i<= nbins) {
//         w = h->GetBinContent(imean+i);
//         x = h->GetBinCenter(imean+i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (sumw > entries) break;
//     }
//     mean = sumwx/sumw; //this would be the mean value of the central 90 percent
//     // 
//     Double_t rms2 = TMath::Abs(sumwx2/sumw);
//     rms_val = TMath::Sqrt(rms2);
//     printf("RMS of central 90 percent = %g, RMS total = %g, Standard deviation = %g\n",rms_val,h->GetRMS(),h->GetStdDev());
// }



