// I've tried to put all the functions which are used repeadetly into another macro. I don't really know how this is usually done.
#include "/home/llr/ilc/ritzmann/work/root_macros/functions_resolution_calc.c"
void sum_hits_notcorr(char* filepath, float input_energy, double_t &mean, double_t &error_val, ofstream &file_notcorr);
void sum_hits_notcorr_low(ofstream &file_notcorr, double_t (&mean_low)[24], double_t (&error_low)[24]);
void sum_hits_corr(char* filepath, TSpline3* spline_fit, float input_energy, ofstream &mean_file, ofstream &median_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms68_file, ofstream &rms90_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit, ofstream &gamma_file);
void calc_median(TTree *tree, TH1F *h, TSpline3* spline_fit, ofstream &median_file, float input_energy);
// void gamma_fit(TH1F *h, float input_energy, ofstream &gamma_file);

// void rms68(TH1F *h, double &mean, double &rms_val);
// void rms90(TH1F *h, double &mean, double &rms_val);
// void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries);
// void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy);
// double calc_ratio(TH1F * hist, int entries);
// double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries);


TSpline3* sum_hits(){
    // energy values of the simulation
    double_t energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    // energy values in case we would also include the values from the continuous simulation
    double_t energies_plot[24] = {0.02, 0.05, 0.08, 0.11, 0.125, 0.15, 0.175, 0.2, 0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};
    // arrays to store the mean values and corresponding errors
    // double_t mean_values[24] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // double_t error_values[24] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t mean_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double_t error_values[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    
    double_t y_errors[16] = {10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7, 10e-7};

    // txt file for mean number of hits values
    ofstream file_notcorr;
    file_notcorr.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_mean_test.txt");
    file_notcorr << "energy,mean,sigma" << "\n";

    // this would additionally consider the points from the continuous energy simulation
    // sum_hits_notcorr_low(file_notcorr, mean_values, error_values);

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
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val, error_val;
        // calculates the mean number of hits for every energy
        sum_hits_notcorr((char*)argument.c_str(), energies[i], mean_val, error_val, file_notcorr);
        mean_values[i] = mean_val;
        error_values[i] = error_val;
        // mean_values[i+8] = mean_val;
        // error_values[i+8] = error_val;
    }
    file_notcorr.close();

    // spline interpolation
    TSpline3 *spline_fit = new TSpline3("Spline Interpolation", mean_values, energies, 16);
    spline_fit->SetLineWidth(2);
    spline_fit->SetLineColor(kRed);
    spline_fit->Draw();

    // polynomial fit
    TCanvas *c_fit = new TCanvas("polynomial_fit", "polynomial_fit");
    // c_fit->DrawFrame(0, 120, 0, 7);
    TGraphErrors *hits = new TGraphErrors(16, mean_values, energies, error_values, 0);
    hits->SetTitle("number of hits;number of hits;energy");
    // spline_fit->Draw("C P");

    // hits->GetYAxis()->SetTitle("energy");
    hits->GetXaxis()->SetRangeUser(0, 100);
    hits->GetYaxis()->SetRangeUser(0, 7);
    // spline_point->Draw();
    spline_fit->Draw();
    hits->Draw("*same");


    TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->AddEntry(spline_fit, "spline", "L");
    legend->AddEntry(hits, "number of hits", "P");
    legend->Draw();

    // this was from the polynomial fit we tried in the beginning
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
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_test.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream median_file;
    median_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_median_test.txt");
    median_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss_test.txt");
    gauss_file << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream gauss_file2s;
    gauss_file2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss2s_test.txt");
    gauss_file2s << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";

    ofstream rms68_file;
    rms68_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_rms68_test.txt");
    rms68_file << "energy,mean,sigma" << "\n";

    ofstream rms90_file;
    rms90_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_rms90_test.txt");
    rms90_file << "energy,mean,sigma" << "\n";

    ofstream myfile_3sig;
    myfile_3sig.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_hits_test.txt");
    myfile_3sig << "energy,ratio" << "\n";

    ofstream myfile_3sig_fit;
    myfile_3sig_fit.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas_sum_hits_fitted_test.txt");
    myfile_3sig_fit << "energy,ratio" << "\n";

    ofstream gamma_file;
    gamma_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gamma_test.txt");
    gamma_file << "energy,mean,sigma,chisq,ndf" << "\n";

    // iterate over all the energies
    for(int i = 0; i < 16; i++) {
        string argument1, argument2, beam_energy;

        // use the correct path
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

        // calculated the energy from the number of hits
        sum_hits_corr((char*)argument.c_str(), spline_fit, energies[i], mean_file, median_file, gauss_file, gauss_file2s, rms68_file, rms90_file, myfile_3sig, myfile_3sig_fit, gamma_file);

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

// to include the values from the continuous energy simulation
void sum_hits_notcorr_low(ofstream &file_notcorr, double_t (&mean_low)[24], double_t (&error_low)[24]){
    TFile *input_ecal = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_0-6GeV_uniform_106.root");
    TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform_106.root");
    
    TTree *tree = (TTree*)input_ecal->Get("ecal");
    TTree *tpart = (TTree*)input_particle->Get("MyLCTuple");

    tree->Draw("nhit_len>>nhits");
    TH1F* h_nhits = (TH1F*)gPad->GetPrimitive("nhits");

    double xmin_hits = h_nhits->GetXaxis()->GetXmin();
    double xmax_hits = h_nhits->GetXaxis()->GetXmax();

    // integer binning with the bin center at the integer number
    int nbins = int(ceil(xmax_hits));
    Double_t xbins[nbins + 1];
    // xbins[0] = -0.5;
    for (int i = 0; i < nbins+1; i++){
        double val = i-0.5;
        xbins[i] = val;
        // cout << i << " " << xbins[i] << endl;
    }

    TH2F* hist = new TH2F("hist", "E(sim) vs. nhits", nbins, xbins, 100, 0, 2);

    int nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    vector<int> *hit_slab = 0;
    tree->SetBranchAddress("hit_slab", &hit_slab);

    Float_t mcene[50];
    tpart->SetBranchAddress("mcene", &mcene);

    int entries_p = tpart->GetEntries();
    int entries_d = tree->GetEntries();
    cout << entries_p << endl;

    for (int i = 0; i < entries_d; i++){
        // cout << i << endl;
        tpart->GetEntry(i);
        tree->GetEntry(i);
        // hist->Fill(spline_fit_hits->Eval(nhit_len),mcene[0]);
        if(nhit_len!=0){
            hist->Fill(nhit_len,mcene[0]);
        }
        // hist->Fill(mcene[0], nhit_len);  
    }

    double energies_low[9] = {0.02, 0.05, 0.08, 0.11, 0.125, 0.15, 0.175, 0.2, 0.25};
    for (int i = 0; i < 9; i++){
        TCanvas *c2 = new TCanvas(("new" + to_string(energies_low[i])).c_str(), to_string(energies_low[i]).c_str());
        int bin = hist->GetYaxis()->FindBin(energies_low[i]);
        // cout << bin << endl;
        TH1D* histy = hist->ProjectionX((to_string(energies_low[i]) + "_py").c_str(), bin, bin);
        // histy->Draw();
        mean_low[i] = histy->GetMean();
        error_low[i] = histy->GetStdDev();

        file_notcorr << energies_low[i] << "," << histy->GetMean() << "," << histy->GetStdDev() << "\n";

        // histy->Fit("gaus");
        // TF1 *fit = histy->GetFunction("gaus");
        // double mean = fit->GetParameter(1);
        // double sigma = fit->GetParameter(2);
        // cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        // gauss_file << energies[i] << "," << mean << "," << sigma << "\n";
        histy->Draw();
    }
}


// calculate the initial estimate, could be mean or also median
void sum_hits_notcorr(char* filepath, float input_energy, double_t &mean, double_t &error_val, ofstream &file_notcorr){
   
    TFile *input = new TFile(filepath, "read");
    TTree *tree = (TTree*)input->Get("ecal");
    double_t entries = tree->GetEntries();

    TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

    tree->Draw("nhit_len");
    TH1F* hist = (TH1F*)gPad->GetPrimitive("htemp");
    mean = hist->GetMean();
    double_t std_hist = hist->GetStdDev();
    error_val = std_hist / TMath::Sqrt(entries);

    gSystem->ProcessEvents();
    gSystem->ProcessEvents();

    int xmin = hist->GetXaxis()->GetXmin();
    int xmax = hist->GetXaxis()->GetXmax();

    double factor;
    if (input_energy <= 10){
        factor = 1;
    }

    // introduce factor for higher energies, such that we don't have too many bins
    if (input_energy > 10){
        factor = input_energy/20;
    }    

    int nbins = int(ceil((xmax-xmin)/factor));
    Double_t xbins[nbins + 1];

    xbins[0] = xmin-0.5;

    for(int i = 0; i < nbins+1; i++){
        double val = xmin + factor*(i + 1);
        xbins[i+1] = val-0.5;
        // cout << i << " " << xbins[i+1] << endl;
    }

    TH1F* hist_not_corr = new TH1F((to_string(input_energy)).c_str(), ("uncorrected hits " + to_string(input_energy)).c_str(), nbins , xbins);

    int nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    for (int i = 0; i < entries; i++){
        tree->GetEntry(i);
        // int nhits = nhit_len->data();
        // exclude the zeros
        if (nhit_len != 0){
            hist_not_corr->Fill(nhit_len);
        }   
    }
    TCanvas *c2 = new TCanvas(("new_bin" + to_string(input_energy)).c_str(), ("new_bin" + to_string(input_energy)).c_str());
    hist_not_corr->Draw();

    file_notcorr << input_energy << "," << mean << "," << std_hist << "\n";
}

// calculate the energy value from the number of hits
void sum_hits_corr(char* filepath, TSpline3* spline_fit, float input_energy, ofstream &mean_file, ofstream &median_file, ofstream &gauss_file, ofstream &gauss_file2s, ofstream &rms68_file, ofstream &rms90_file, ofstream &myfile_3sig, ofstream &myfile_3sig_fit, ofstream &gamma_file){
  
    TFile *input = new TFile(filepath, "read");
    TTree *tree = (TTree*)input->Get("ecal");
    double_t entries = tree->GetEntries(); 
    cout << "Number of events:" << entries << endl;

    TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

    tree->Draw("nhit_len");
    TH1F* hist0 = (TH1F*)gPad->GetPrimitive("htemp");
    double_t mean0 = hist0->GetMean();
    double_t sigma0 = hist0->GetStdDev();

    int xmin = hist0->GetXaxis()->GetXmin();
    int xmax = hist0->GetXaxis()->GetXmax();

    cout << xmin << " " << xmax << endl;
    cout << spline_fit->Eval(xmin) << " " << spline_fit->Eval(xmax) << endl;

    // define a new binning according to the spline interpolation
    double factor;
    if (input_energy <= 10){
        factor = 1;
    }

    if (input_energy > 10){
        factor = input_energy/20;
    }

    int nbins = ceil((xmax-xmin)/factor);
    cout << "nbins: " << nbins << endl;
    Double_t xbins[nbins + 1];
    xbins[0] = spline_fit->Eval(xmin);
        
    for (int i = 0; i < nbins+1; i++){
        double val = spline_fit->Eval(xmin+factor*(i+1)-0.5);
        xbins[i+1] = val;
        // cout << i << " " << xbins[i+1] << endl;
        }

    TH1F* hist_corr = new TH1F((to_string(input_energy)).c_str(), ("corrected hits " + to_string(input_energy)).c_str(), nbins , xbins);
    
    int nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    for (int i = 0; i < entries; i++){
        tree->GetEntry(i);
        if (nhit_len!=0){
            hist_corr->Fill(spline_fit->Eval(nhit_len));
        } 
    }
    hist_corr->GetXaxis()->SetTitle("Energy [GeV]");
    hist_corr->GetYaxis()->SetTitle("number of events");
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

    calc_median(tree, hist_corr, spline_fit, median_file, input_energy);
        
    fit_gauss(hist_corr, gauss_file, input_energy, myfile_3sig_fit, entries);
    fit_gauss_2sig(hist2_sig, gauss_file2s, mean_hist, std_hist, input_energy);

    gamma_fit(hist_corr, input_energy, gamma_file);
    }

    void calc_median(TTree *tree, TH1F *h, TSpline3* spline_fit, ofstream &median_file, float input_energy){
        double median, MAD;
        double percentile = 0.5;
        h->ComputeIntegral();
        h->GetQuantiles(1, &median, &percentile);

        double_t entries = tree->GetEntries(); 
        cout << "Number of events:" << entries << endl;

        double xmin_mad = h->GetBinLowEdge(1) - median;
        double xmax_mad = h->GetBinLowEdge(h->GetNbinsX()+1) - median;

        TH1F* hist_mad = new TH1F(("hist_mad" + to_string(input_energy)).c_str(), ("hist_mad"+to_string(input_energy)).c_str(), 50 , xmin_mad, xmax_mad);
        
        int nhit_len = 0;
        tree->SetBranchAddress("nhit_len", &nhit_len);

        for (int i = 0; i < entries; i++){
            tree->GetEntry(i);
            hist_mad->Fill(abs(spline_fit->Eval(nhit_len) - median));
        }

        hist_mad->ComputeIntegral();
        hist_mad->GetQuantiles(1, &MAD, &percentile);

        median_file << input_energy << "," << median << "," << 1.4826 * MAD << "\n";
    }

//     void gamma_fit(TH1F *h, float input_energy, ofstream &gamma_file){

//         TCanvas *c = new TCanvas(("gamma" + to_string(input_energy)).c_str(), ("gamma" + to_string(input_energy)).c_str());

//         TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-TMath::Power([1]/[2],2)*((x-[1])/[1]-TMath::Log(x/[1])))", 0, h->GetBinLowEdge(h->GetNbinsX()) + 20); // "xmin" = 0, "xmax" = 100
//         f1->SetParameters(h->GetMaximum(), h->GetMean(), h->GetStdDev()); // you MUST set non-zero initial values for parameters
//         f1->SetParNames("A", "x0", "s0");
//         h->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
//         double A = f1->GetParameter("A");
//         double x0 = f1->GetParameter("x0");
//         double s0 = f1->GetParameter("s0");
//         gStyle->SetOptFit(11111);
//         double chisq = f1->GetChisquare();
//         double ndf = f1->GetNDF();

//         gamma_file << input_energy << "," << x0 + TMath::Power(s0,2)/x0 << "," << s0*TMath::Sqrt(1 + TMath::Power(s0,2)/TMath::Power(x0,2)) << "," << chisq << "," << ndf << "\n";

//         // cout << x0 + TMath::Power(s0,2)/x0 << " " << h->GetMean() << endl;
//         // cout << s0*TMath::Sqrt(1 + TMath::Power(s0,2)/TMath::Power(x0,2)) << " " << h->GetStdDev() << endl;
//         h->Draw();
        
//     }

// void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries){
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

//     double ratio = calc_ratio_fitted(h, mean, sigma, entries);
//     myfile_3sig_fit << input_energy << "," << ratio << "\n";
// }

// void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy){
//     h->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
//     TF1 *fit = h->GetFunction("gaus");
//     double_t mean = fit->GetParameter(1);
//     double sigma = fit->GetParameter(2);
//     double mean_err = fit->GetParError(1);
//     double sigma_err = fit->GetParError(2);
//     float res = sigma/mean;
//     double chisq = fit->GetChisquare();
//     double ndf = fit->GetNDF();
//     double goodness = chisq / ndf;
//     cout << "resolution:" << res << endl;
//     cout << "goodness-of-fit:" << goodness << endl;
//     h->Draw();
//     gStyle->SetOptFit(11111);
//     gSystem->ProcessEvents();
//     gSystem->ProcessEvents();
//     gauss_file2s << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
// }

// double calc_ratio(TH1F * hist, int entries){
//     double_t mean = hist->GetMean();
//     double_t sigma = hist->GetStdDev();
    
//     int nbin_left = hist->GetXaxis()->FindBin(mean - 3*sigma);
//     int nbin_right = hist->GetXaxis()->FindBin(mean + 3*sigma) + 1;

//     double_t ex_entries = 0;

//     for(int i = 0; i<nbin_left; i++){
//         ex_entries = ex_entries + hist->GetBinContent(i);
//     }

//     for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
//         ex_entries = ex_entries + hist->GetBinContent(i);
//     }

//     return ex_entries / entries;
// }

// double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries){
//     int nbin_left = hist->GetXaxis()->FindBin(mean_fit - 3*sigma_fit);
//     int nbin_right = hist->GetXaxis()->FindBin(mean_fit + 3*sigma_fit) + 1;

//     double_t ex_entries = 0;

//     for(int i = 0; i<nbin_left; i++){
//         ex_entries = ex_entries + hist->GetBinContent(i);
//     }

//     for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
//         ex_entries = ex_entries + hist->GetBinContent(i);
//     }

//     return ex_entries / entries;
// }

// void rms90(TH1F *h, double &mean, double &rms_val) {
//     double median, percentile;
//     percentile = 0.5;
//     h->ComputeIntegral();
//     h->GetQuantiles(1, &median, &percentile);
//     TAxis *axis = h->GetXaxis();
//     Int_t nbins = axis->GetNbins();
//     Int_t imedian = axis->FindBin(median);
//     Double_t entries = 0.9*h->GetEntries();
//     Double_t w = h->GetBinContent(imedian);
//     Double_t x = h->GetBinCenter(imedian);
//     Double_t sumw = w;
//     Double_t sumwx = w*x;
//     Double_t sumwx2 = w*x*x;
//     Int_t range;
//     for (Int_t i = 1; i<nbins;i++) {
//     if (i> 0) {
//         w = h->GetBinContent(imedian-i);
//         x = h->GetBinCenter(imedian-i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (i<= nbins) {
//         w = h->GetBinContent(imedian+i);
//         x = h->GetBinCenter(imedian+i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (sumw > entries) {
//         range = i;
//         break;
//     }
//     }
//     cout << "median bin: " << imedian << " " << "range: " << range << endl;
//     Double_t x1 = h->GetBinLowEdge(imedian-range);
//     Double_t x2 = h->GetBinLowEdge(imedian+range+1);

//     TH1F* hist90 = new TH1F("h90", "h90", 2*range, x1, x2);
//     for(Int_t j = 1; j <= hist90->GetNbinsX(); j++){
//         hist90->SetBinContent(j, h->GetBinContent(imedian-range + j));
//     }

//     // hist68->Draw();
//     mean = hist90->GetMean();
//     rms_val = hist90->GetStdDev(); //this would be the mean value of the central 90 percent
//     // 
//     // Double_t rms2 = TMath::Abs(sumwx2/sumw);
//     // rms_val = TMath::Sqrt(rms2);
//     printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
// }

// void rms68(TH1F *h, double &mean, double &rms_val) {
//     double median, percentile;
//     percentile = 0.5;
//     h->ComputeIntegral();
//     h->GetQuantiles(1, &median, &percentile);
//     TAxis *axis = h->GetXaxis();
//     Int_t nbins = axis->GetNbins();
//     Int_t imedian = axis->FindBin(median);
//     Double_t entries = 0.68*h->GetEntries();
//     Double_t w = h->GetBinContent(imedian);
//     Double_t x = h->GetBinCenter(imedian);
//     Double_t sumw = w;
//     Double_t sumwx = w*x;
//     Double_t sumwx2 = w*x*x;
//     Int_t range;
//     for (Int_t i = 1; i<nbins;i++) {
//     if (i> 0) {
//         w = h->GetBinContent(imedian-i);
//         x = h->GetBinCenter(imedian-i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (i<= nbins) {
//         w = h->GetBinContent(imedian+i);
//         x = h->GetBinCenter(imedian+i);
//         sumw += w;
//         sumwx += w*x;
//         sumwx2 += w*x*x;
//     }
//     if (sumw > entries) {
//         range = i;
//         break;
//     }
//     }
//     cout << "median bin: " << imedian << " " << "range: " << range << endl;
//     Double_t x1 = h->GetBinLowEdge(imedian-range);
//     Double_t x2 = h->GetBinLowEdge(imedian+range+1);

//     TH1F* hist90 = new TH1F("h90", "h90", 2*range, x1, x2);
//     for(Int_t j = 1; j <= hist90->GetNbinsX(); j++){
//         hist90->SetBinContent(j, h->GetBinContent(imedian-range + j));
//     }

//     // hist68->Draw();
//     mean = hist90->GetMean();
//     rms_val = hist90->GetStdDev(); //this would be the mean value of the central 90 percent
//     // 
//     // Double_t rms2 = TMath::Abs(sumwx2/sumw);
//     // rms_val = TMath::Sqrt(rms2);
//     printf("RMS of central 68 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
// }



