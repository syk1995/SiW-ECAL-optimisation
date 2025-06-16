
double_t sum_hits_notcorr(char* filename, float input_energy){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_mean.txt", std::ios_base::app);

    //input of the files fabricio gave me
    // TFile *input = new TFile("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_80GeV_5kevt_-42_-42_build_masked.root", "read");
    //TFile *input;

    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");


    double_t mean;

    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("nhit_len"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double std_hist = hist->GetStdDev();


        // //get the number of events
        // int entries = tree->GetEntries();
        
        // cout << "Number of events:" << entries << endl;
        
        // // hist->Fit("gaus");
        // hist->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        // TF1 *fit = hist->GetFunction("gaus");
        // mean = fit->GetParameter(1);
        // double sigma = fit->GetParameter(2);
        // double mean_err = fit->GetParError(1);
        // double sigma_err = fit->GetParError(2);
        // float res = sigma/mean;
        // double chi_squ = fit->GetChisquare();
        // int bins = hist->GetNbinsX();
        // double goodness = chi_squ / (bins - 3);
        // cout << "resolution:" << res << endl;
        // cout << "goodness-of-fit:" << goodness << endl;
        // gStyle->SetOptFit(11111);
        // hist->Draw();
        // // c->Update();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 

    }


    else if ((input_energy >= 1) && (input_energy < 200)){   
        int int_input_energy = (int)input_energy; 
        cout << "small " << int_input_energy << endl;
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("nhit_len"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double std_hist = hist->GetStdDev();

        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }


    else {
        cout << "big" << endl;
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("nhit_len"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        mean = hist->GetMean();
        double std_hist = hist->GetStdDev();

        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << std_hist << "\n";
        myfile.close(); 
    }
    return mean;
}


double_t a0 = 0;
double_t a1 = 0;
double_t a2 = 0;
double_t a3 = 0;
double_t a4 = 0;

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

double_t sum_hits_notcorr(char* filename, int input_energy);

void rms68(TH1F *h, double &mean, double &rms_val);




void sum_hits_corr(char* filename, float input_energy){
    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile_mean;
    myfile_mean.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean.txt", std::ios_base::app);

    ofstream myfile_gauss;
    myfile_gauss.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss.txt", std::ios_base::app);
    
    ofstream myfile_gauss2s;
    myfile_gauss2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss2s.txt", std::ios_base::app);
    
    ofstream myfile_gauss2s;
    myfile_gauss2s.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss2s.txt", std::ios_base::app);
    

    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("corr_c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());


        tree->Draw("enhits(nhit_len)"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        // double mean_hist = hist2->GetMean();
        // double std_hist = hist2->GetStdDev();

        // //get the number of events
        // int entries = tree->GetEntries();
        
        // cout << "Number of events:" << entries << endl;

        // // hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        // TF1 *fit2 = hist2->GetFunction("gaus");
        // double_t mean = fit2->GetParameter(1);
        // double sigma = fit2->GetParameter(2);
        // double mean_err = fit2->GetParError(1);
        // double sigma_err = fit2->GetParError(2);
        // float res = sigma/mean;
        // double chi_squ = fit2->GetChisquare();
        // int bins = hist2->GetNbinsX();
        // double goodness = chi_squ / (bins - 3);
        // cout << "resolution:" << res << endl;
        // cout << "goodness-of-fit:" << goodness << endl;
        // hist2->Draw();
        // gStyle->SetOptFit(11111);
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        double mean, rms_val;

        rms68(hist2, mean, rms_val);

        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << rms_val << "\n";
        myfile.close(); 

    }
    
    
    
    else if ((input_energy >= 1) && (input_energy < 200)){  
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 
        
        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("corr_c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("enhits(nhit_len)"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        // double mean_hist = hist2->GetMean();
        // double std_hist = hist2->GetStdDev();

        // //get the number of events
        // int entries = tree->GetEntries();
        
        // cout << "Number of events:" << entries << endl;

        // // hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        // TF1 *fit2 = hist2->GetFunction("gaus");
        // double_t mean = fit2->GetParameter(1);
        // double sigma = fit2->GetParameter(2);
        // double mean_err = fit2->GetParError(1);
        // double sigma_err = fit2->GetParError(2);
        // float res = sigma/mean;
        // double chi_squ = fit2->GetChisquare();
        // int bins = hist2->GetNbinsX();
        // double goodness = chi_squ / (bins - 3);
        // cout << "resolution:" << res << endl;
        // cout << "goodness-of-fit:" << goodness << endl;
        // hist2->Draw();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        double mean, rms_val;

        rms68(hist2, mean, rms_val);



        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << rms_val << "\n";
        myfile.close(); 
    }

    else {
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("corr_c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("enhits(nhit_len)"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        // double mean_hist = hist2->GetMean();
        // double std_hist = hist2->GetStdDev();


        // //get the number of events
        // int entries = tree->GetEntries();
        
        // cout << "Number of events:" << entries << endl;

        // // hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        // TF1 *fit2 = hist2->GetFunction("gaus");
        // double_t mean = fit2->GetParameter(1);
        // double sigma = fit2->GetParameter(2);
        // double mean_err = fit2->GetParError(1);
        // double sigma_err = fit2->GetParError(2);
        // float res = sigma/mean;
        // double chi_squ = fit2->GetChisquare();
        // int bins = hist2->GetNbinsX();
        // double goodness = chi_squ / (bins - 3);
        // cout << "resolution:" << res << endl;
        // cout << "goodness-of-fit:" << goodness << endl;
        // hist2->Draw();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        double mean, rms_val;

        rms68(hist2, mean, rms_val);


        // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
        myfile << input_energy << "," << mean << "," << rms_val << "\n";
        myfile.close(); 
    }
    
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");




    //tree->Draw("a3*TMath::Power(nhit_len,3)+a2*TMath::Power(nhit_len,2)+a1*nhit_len+a0"); 

    // double_t k[4] = {2.90894e-09, 0.000107585, 0.0467945, -0.573716};

    // auto random = [a3, a2, a1, a0](int number_hits){
    //     return a3*pow(number_hits,3)+a2*pow(number_hits,2)+a1*number_hits+a0;
    // };

    // cout << random(50) << endl;

    // tree->Draw("enhits(nhit_len)"); 
    // auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");

    // //get the number of events
    // int entries = tree->GetEntries();
    
    // cout << "Number of events:" << entries << endl;

    // hist2->Fit("gaus");
    // TF1 *fit2 = hist2->GetFunction("gaus");
    // double_t mean = fit2->GetParameter(1);
    // double sigma = fit2->GetParameter(2);
    // double mean_err = fit2->GetParError(1);
    // double sigma_err = fit2->GetParError(2);
    // float res = sigma/mean;
    // double chi_squ = fit2->GetChisquare();
    // int bins = hist2->GetNbinsX();
    // double goodness = chi_squ / (bins - 3);
    // cout << "resolution:" << res << endl;
    // cout << "goodness-of-fit:" << goodness << endl;
    // hist2->Draw();

    // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    // myfile.close(); 
}




void sum_hits_new(){

    float energies[15] = {0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t energies_plot[15] = {0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    // double_t energies_plot[13] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t mean_values[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // double_t mean_values[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_mean.txt");
    // myfile_events << "energy,mean,mean_err,sigma,sigma_err, goodness" << "\n";
    myfile_events << "energy,mean,sigma" << "\n";
    myfile_events.close();

    for (int i = 0; i < 2; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        cout << argument << endl;

        double_t mean_val = sum_hits_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }
    
    
    for (int i = 2; i < 13; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_hits_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }

    for (int i = 13; i < 15; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_hits_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }
    
    TCanvas *c_fit = new TCanvas("gaussian_fit", "gaussian_fit");


    TGraph *hits = new TGraph(9, mean_values, energies_plot);

    hits->Fit("pol4");
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
    // a0 = -9.64060835e-02;
    // a1 = 3.42225914e-02;
    // a2 = 1.95823398e-04;
    // a3 = -2.35185160e-07;
    // a4 = 2.23654647e-10;

    cout << a0 << a1<< a2 << a3 << a4 << endl;

    for (const auto& e : mean_values) {
    std::cout << e << " " << enhits(e) << std::endl;
    }
    cout << a0 << a1<< a2 << a3 << a4 << endl;


    double_t mean_values_new[11];

    //just summing the events, without applying masking
    ofstream myfile_events2;
    myfile_events2.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_rms68.txt");
    // myfile_events2 << "energy,mean,mean_err,sigma,sigma_err,goodness" << "\n";
    myfile_events2 << "energy,mean,sigma" << "\n";
    myfile_events2.close();

    for (int i = 0; i < 2; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i]);

    }
    
    
    for (int i = 2; i < 13; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i]);

    }

    for (int i = 13; i < 15; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_hits_corr((char*)argument.c_str(), energies[i]);

    }  
}

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
    for (Int_t i=1;i<nbins;i++) {
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
    if (sumw > entries) break;
    }
    cout << i << endl;
    Double_t x1 = GetBinLowEdge(maxbin-i);
    Double_t x2 = GetBinLowEdge(maxbin+i+1)
    Int_t number_bins = ceil(0.68 * nbins);

    TH1F *hist68("h68", "h68", 2*i, x1, x2);
    mean = hist68->GetMean();
    rms_val = hist68->GetStdDev(); //this would be the mean value of the central 90 percent
    // 
    // Double_t rms2 = TMath::Abs(sumwx2/sumw);
    // rms_val = TMath::Sqrt(rms2);
    // printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}




