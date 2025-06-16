
void sum_hits_notcorr(char* filename, float input_energy, double_t &mean, double_t &error_val){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_mean.txt", std::ios_base::app);


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

        tree->Draw("nhit_len"); 
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


        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("nhit_len"); 
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
        double_t entries = tree->GetEntries();

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("nhit_len"); 
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


double_t a0 = 0;
double_t a1 = 0;
double_t a2 = 0;
double_t a3 = 0;
double_t a4 = 0;

TSpline3 *spline_fit;

double spline_hits(double number_hits){
    double val = spline_fit->Eval(number_hits);
    return val;
}

double_t enhits(double number_hits){

    // a0 = -9.64060835e-02;
    // a1 = 3.42225914e-02;
    // a2 = 1.95823398e-04;
    // a3 = -2.35185160e-07;
    // a4 = 2.23654647e-10;

    // cout << a0 << a1 << a2 << a3 << a4 << endl;
    
    // double_t val = a3*pow(number_hits,3)+a2*pow(number_hits,2)+a1*number_hits+a0;
    double_t val = a4*pow(number_hits,4)+a3*pow(number_hits,3)+a2*pow(number_hits,2)+a1*number_hits+a0;
    return val;
}

double_t sum_hits_notcorr(char* filename, int input_energy, double_t &mean, double_t &error_val);

void rms68(TH1F *h, double &mean, double &rms_val);
void rms90(TH1F *h, double &mean, double &rms_val);
void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries);
void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy);
double calc_ratio(TH1F * hist, int entries);
double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries);



void sum_hits_corr(char* filename, float input_energy, ofstream &mean_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms_file, TSpline3* spline_fit, ofstream &myfile_3sig, ofstream &myfile_3sig_fit, ofstream &median_file){
    //write the fit parameters to a txt file
    //for the e- files from fabricio
    
    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        //get the number of events
        int entries = tree->GetEntries();
        cout << "Number of events:" << entries << endl;

        TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());
        // tree->Draw("enhits(nhit_len)"); 
        tree->Draw("spline_hits(nhit_len)");
        TH1F* hist2 = (TH1F*)gPad->GetPrimitive("htemp");

        cout << hist2->GetNbinsX() << endl;

        TH1F* hist2_sig = (TH1F*)hist2->Clone();

        double ratio_3sig = calc_ratio(hist2, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";

        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        rms90(hist2, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";
        rms_file << input_energy << "," << mean_rms << "," << 1.23172 * rms_val << "\n";


        double median, MAD, percentile;
        percentile = 0.5;
        hist2->ComputeIntegral();
        hist2->GetQuantiles(1, &median, &percentile);

        tree->Draw("(spline_hits(nhit_len) - median)>>minus_median");
        
        TH1F *minus_median = (TH1F*)gPad->GetPrimitive("minus_median");

        median_file << input_energy << "," << median << "\n";

        
        
        fit_gauss(hist2, gauss_file, input_energy, myfile_3sig_fit, entries);
        fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);

        // TCanvas *rms = new TCanvas(("rms_68" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());
        
    }
    
    
    
    else if ((input_energy >= 1) && (input_energy < 200)){  
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 
        
        TTree *tree = (TTree*)input->Get("ecal");

        //get the number of events
        int entries = tree->GetEntries();
        cout << "Number of events:" << entries << endl;

        TCanvas *c = new TCanvas(("corr_c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        // tree->Draw("enhits(nhit_len)"); 
        tree->Draw("spline_hits(nhit_len)");
        TH1F* hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        double ratio_3sig = calc_ratio(hist2, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";
        cout << hist2->GetNbinsX() << endl;

        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        rms90(hist2, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";
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

        // tree->Draw("enhits(nhit_len)"); 
        tree->Draw("spline_hits(nhit_len)");
        TH1F* hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        double ratio_3sig = calc_ratio(hist2, entries);
        myfile_3sig << input_energy << "," << ratio_3sig << "\n";
        cout << hist2->GetNbinsX() << endl;

        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();
        mean_file << input_energy << "," << mean_hist << "," << std_hist << "\n";

        double mean_rms, rms_val;

        rms90(hist2, mean_rms, rms_val);
        // rms_file << input_energy << "," << mean_rms << "," << 1.76 * rms_val << "\n";
        rms_file << input_energy << "," << mean_rms << "," << 1.23172 * rms_val << "\n";

        TH1F* hist2_sig = (TH1F*)hist2->Clone();

        fit_gauss(hist2, gauss_file, input_energy, myfile_3sig_fit, entries);
        fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);
    }
}

void sum_hits_new(){

    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t energies_plot[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    double_t mean_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t error_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t y_errors[16] = {10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7};

    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_mean.txt");
    // myfile_events << "energy,mean,mean_err,sigma,sigma_err, goodness" << "\n";
    myfile_events << "energy,mean,sigma" << "\n";
    myfile_events.close();

    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;
        // cout << argument << endl;
        double_t mean_val, error_val;
        sum_hits_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
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
        sum_hits_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
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
        sum_hits_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
    }
    
    TCanvas *c_fit = new TCanvas("gaussian_fit", "gaussian_fit");

    spline_fit = new TSpline3("spline_interpol", mean_values, energies_plot, 16);
    cout << spline_fit->Eval(1287.74) << endl;

    TGraphErrors *hits = new TGraphErrors(9, mean_values, energies_plot, error_values, 0);
    TF1 *gr_fit = ((TF1 *)(gROOT->GetFunction("pol4")));
    hits->Fit(gr_fit, "WQ"); // "initial pre-fit"
    hits->Fit(gr_fit, "Q"); // "final fit"
    
    // TGraph *hits = new TGraph(9, mean_values, energies_plot);
    // hits->Fit("pol4");
    
    hits->Draw("APL*");
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    // hits->GetXaxis()->SetLimits(0,400);
    TF1 *fit_pol = hits->GetFunction("pol4");
    a0 = fit_pol->GetParameter(0);
    a1 = fit_pol->GetParameter(1);
    a2 = fit_pol->GetParameter(2);
    a3 = fit_pol->GetParameter(3);
    a4 = fit_pol->GetParameter(4);

    // cout << a0 << a1<< a2 << a3 << a4 << endl;

    for (const auto& e : mean_values) {
    std::cout << e << " " << enhits(e) << std::endl;
    std::cout << e << " " << spline_hits(e) << std::endl;
    }
    cout << a0 << a1<< a2 << a3 << a4 << endl;

    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss.txt");
    gauss_file << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream gauss_file2s;
    gauss_file2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss2s.txt");
    gauss_file2s << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream rms_file;
    rms_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_rms90.txt");
    rms_file << "energy,mean,sigma" << "\n";

    ofstream median_file;
    median_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_median.txt");
    median_file << "energy,mean,sigma" << "\n";

    ofstream myfile_3sig;
    myfile_3sig.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_hits.txt");
    myfile_3sig << "energy,ratio" << "\n";

    ofstream myfile_3sig_fit;
    myfile_3sig_fit.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_hits_fitted.txt");
    myfile_3sig_fit << "energy,ratio" << "\n";
    
    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, spline_fit, myfile_3sig, myfile_3sig_fit, median_file);
    }
    
    
    for (int i = 3; i < 14; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, spline_fit, myfile_3sig, myfile_3sig_fit, median_file);
    }

    for (int i = 14; i < 16; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i], mean_file, gauss_file, gauss_file2s, rms_file, spline_fit, myfile_3sig, myfile_3sig_fit, median_file);
    }  
    mean_file.close();
    gauss_file.close();
    gauss_file2s.close();
    rms_file.close();
    myfile_3sig.close();
    myfile_3sig_fit.close();
    median_file.close();
}



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

void rms68(TH1F *h, double &mean, double &rms_val) {
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t maxbin = h->GetMaximumBin();
    Int_t imean = axis->FindBin(h->GetMean());
    Double_t entries = 0.68*h->GetEntries();
    Double_t w = h->GetBinContent(maxbin);
    Double_t x = h->GetBinCenter(maxbin);
    Double_t sumw = w;
    Double_t sumwx = w*x;
    Double_t sumwx2 = w*x*x;
    Int_t range;
    for (Int_t i = 1; i<nbins;i++) {
    if (i> 0) {
        w = h->GetBinContent(maxbin-i);
        x = h->GetBinCenter(maxbin-i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (i<= nbins) {
        w = h->GetBinContent(maxbin+i);
        x = h->GetBinCenter(maxbin+i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (sumw > entries) {
        range = i;
        break;
    }
    }
    cout << "maxbin: " << maxbin << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(maxbin-range);
    Double_t x2 = h->GetBinLowEdge(maxbin+range+1);

    TH1F* hist68 = new TH1F("h68", "h68", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist68->GetNbinsX(); j++){
        hist68->SetBinContent(j, h->GetBinContent(maxbin-range + j));
    }

    // hist68->Draw();
    mean = hist68->GetMean();
    rms_val = hist68->GetStdDev(); //this would be the mean value of the central 90 percent
    // 
    // Double_t rms2 = TMath::Abs(sumwx2/sumw);
    // rms_val = TMath::Sqrt(rms2);
    printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}

void rms90(TH1F *h, double &mean, double &rms_val) {
    double median, percentile;
    percentile = 0.5;
    h->ComputeIntegral();
    h->GetQuantiles(1, &median, &percentile);
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t imedian = axis->FindBin(median);
    Double_t entries = 0.9*h->GetEntries();
    Double_t w = h->GetBinContent(imedian);
    Double_t x = h->GetBinCenter(imedian);
    Double_t sumw = w;
    Double_t sumwx = w*x;
    Double_t sumwx2 = w*x*x;
    Int_t range;
    for (Int_t i = 1; i<nbins;i++) {
    if (i> 0) {
        w = h->GetBinContent(imedian-i);
        x = h->GetBinCenter(imedian-i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (i<= nbins) {
        w = h->GetBinContent(imedian+i);
        x = h->GetBinCenter(imedian+i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (sumw > entries) {
        range = i;
        break;
    }
    }
    cout << "median bin: " << imedian << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(imedian-range);
    Double_t x2 = h->GetBinLowEdge(imedian+range+1);

    TH1F* hist90 = new TH1F("h90", "h90", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist90->GetNbinsX(); j++){
        hist90->SetBinContent(j, h->GetBinContent(imedian-range + j));
    }

    // hist68->Draw();
    mean = hist90->GetMean();
    rms_val = hist90->GetStdDev(); //this would be the mean value of the central 90 percent
    // 
    // Double_t rms2 = TMath::Abs(sumwx2/sumw);
    // rms_val = TMath::Sqrt(rms2);
    printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}

void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries){
    h->Fit("gaus");
    TF1 *fit = h->GetFunction("gaus");
    double_t mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double mean_err = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    float res = sigma/mean;
    double goodness = fit->GetChisquare() / fit->GetNDF();
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    h->Draw();
    gStyle->SetOptFit(11111);
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    gauss_file << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";

    double ratio = calc_ratio_fitted(h, mean, sigma, entries);
    myfile_3sig_fit << input_energy << "," << ratio << "\n";
}

void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy){
    h->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
    TF1 *fit = h->GetFunction("gaus");
    double_t mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double mean_err = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    float res = sigma/mean;
    double chisq = fit->GetChisquare();
    double ndf = fit->GetNDF();
    double goodness = chisq / ndf;
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    h->Draw();
    gStyle->SetOptFit(11111);
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    gauss_file2s << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
}

double calc_ratio(TH1F * hist, int entries){
    double_t mean = hist->GetMean();
    double_t sigma = hist->GetStdDev();
    
    int nbin_left = hist->GetXaxis()->FindBin(mean - 3*sigma);
    int nbin_right = hist->GetXaxis()->FindBin(mean + 3*sigma) + 1;

    double_t ex_entries = 0;

    for(int i = 0; i<nbin_left; i++){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    return ex_entries / entries;
}

double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries){
    int nbin_left = hist->GetXaxis()->FindBin(mean_fit - 3*sigma_fit);
    int nbin_right = hist->GetXaxis()->FindBin(mean_fit + 3*sigma_fit) + 1;

    double_t ex_entries = 0;

    for(int i = 0; i<nbin_left; i++){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    return ex_entries / entries;
}





