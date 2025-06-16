
double_t sum_energy_notcorr(char* filename, float input_energy){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_notcorr.txt", std::ios_base::app);

    //input of the files fabricio gave me
    // TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

    double_t mean;
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");

if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        double mean_hist = hist->GetMean();
        double std_hist = hist->GetStdDev();

        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;

        // hist->Fit("gaus");
        hist->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit = hist->GetFunction("gaus");
        mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        double mean_err = fit->GetParError(1);
        double sigma_err = fit->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit->GetChisquare();
        int deg_free = fit->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist->Draw();
        // c->Update();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 


    }


    else if ((input_energy >= 1) && (input_energy < 200)){   
        int int_input_energy = (int)input_energy; 
        cout << "small " << int_input_energy << endl;
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
        double mean_hist = hist->GetMean();
        double std_hist = hist->GetStdDev();

        cout << mean_hist << " " << std_hist << endl;

        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;

        // hist->Fit("gaus");
        hist->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit = hist->GetFunction("gaus");
        mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        double mean_err = fit->GetParError(1);
        double sigma_err = fit->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit->GetChisquare();
        int deg_free = fit->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist->Draw();
        // c->Update();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << int_input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 
    }


    else {
        cout << "big" << endl;
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        int int_input_energy = (int)input_energy; 

        TTree *tree = (TTree*)input->Get("ecal");

        TCanvas *c = new TCanvas(("c_" + to_string(int_input_energy)).c_str(), to_string(int_input_energy).c_str());

        tree->Draw("Sum$(hit_energy)"); 
        auto hist = (TH1F*)gPad->GetPrimitive("htemp");
         double mean_hist = hist->GetMean();
        double std_hist = hist->GetStdDev();


        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;
        hist->Fit("gaus");
        // hist->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit = hist->GetFunction("gaus");
        mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        double mean_err = fit->GetParError(1);
        double sigma_err = fit->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit->GetChisquare();
        int deg_free = fit->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist->Draw();
        // c->Update();
        gSystem->ProcessEvents();
        gSystem->ProcessEvents();

        myfile << int_input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 
    }

    // TTree *tree = (TTree*)input->Get("ecal");

    // tree->Draw("Sum$(hit_energy)"); 
    // auto hist = (TH1F*)gPad->GetPrimitive("htemp");

    // //get the number of events
    // int entries = tree->GetEntries();
    
    // cout << "Number of events:" << entries << endl;

    // hist->Fit("gaus");
    // TF1 *fit = hist->GetFunction("gaus");
    // double_t mean = fit->GetParameter(1);
    // double sigma = fit->GetParameter(2);
    // double mean_err = fit->GetParError(1);
    // double sigma_err = fit->GetParError(2);
    // float res = sigma/mean;
    // double chi_squ = fit->GetChisquare();
    // int deg_free = fit->GetNDF();
    // double goodness = chi_squ / (deg_free);
    // cout << "resolution:" << res << endl;
    // cout << "goodness-of-fit:" << goodness << endl;
    // hist->Draw();

    // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    // myfile.close(); 

    return mean;
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

double_t sum_energy_notcorr(char* filename, int input_energy);




void sum_energy_corr(char* filename, float input_energy, double a0, double a1, double a2, double a3){
    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr.txt", std::ios_base::app);


    if (input_energy < 1) {
        cout << "very small" << endl;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        // tree->Draw("eesum(Sum$(hit_energy))");
        tree->Draw("Sum$(eesum(hit_energy))"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();

        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;

        hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit2 = hist2->GetFunction("gaus");
        double_t mean = fit2->GetParameter(1);
        double sigma = fit2->GetParameter(2);
        double mean_err = fit2->GetParError(1);
        double sigma_err = fit2->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit2->GetChisquare();
        int deg_free = fit2->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist2->Draw();

        myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 

    }
    
    
    
    else if ((input_energy >= 1) && (input_energy < 200)){  
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        // tree->Draw("eesum(Sum$(hit_energy))"); 
        tree->Draw("Sum$(eesum(hit_energy))"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();

        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;

        hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit2 = hist2->GetFunction("gaus");
        double_t mean = fit2->GetParameter(1);
        double sigma = fit2->GetParameter(2);
        double mean_err = fit2->GetParError(1);
        double sigma_err = fit2->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit2->GetChisquare();
        int deg_free = fit2->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist2->Draw();

        myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 
    }

    else {
        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        // tree->Draw("eesum(Sum$(hit_energy))"); 
        tree->Draw("Sum$(eesum(hit_energy))"); 
        auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");
        double mean_hist = hist2->GetMean();
        double std_hist = hist2->GetStdDev();


        //get the number of events
        int entries = tree->GetEntries();
        
        cout << "Number of events:" << entries << endl;

        hist2->Fit("gaus");
        // hist2->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
        TF1 *fit2 = hist2->GetFunction("gaus");
        double_t mean = fit2->GetParameter(1);
        double sigma = fit2->GetParameter(2);
        double mean_err = fit2->GetParError(1);
        double sigma_err = fit2->GetParError(2);
        float res = sigma/mean;
        double chi_squ = fit2->GetChisquare();
        int deg_free = fit2->GetNDF();
        double goodness = chi_squ / (deg_free);
        cout << "resolution:" << res << endl;
        cout << "goodness-of-fit:" << goodness << endl;
        hist2->Draw();

        myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
        myfile.close(); 
    }

    //input of the files fabricio gave me
    // TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

    
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");



    // TTree *tree = (TTree*)input->Get("ecal");

    //tree->Draw("a3*TMath::Power(nhit_len,3)+a2*TMath::Power(nhit_len,2)+a1*nhit_len+a0"); 

    // double_t k[4] = {2.90894e-09, 0.000107585, 0.0467945, -0.573716};

    // auto random = [a3, a2, a1, a0](int number_hits){
    //     return a3*pow(number_hits,3)+a2*pow(number_hits,2)+a1*number_hits+a0;
    // };

    // cout << random(50) << endl;

    // tree->Draw("eesum(Sum$(hit_energy))"); 
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
    // int deg_free = fit2->GetNDF();
    // double goodness = chi_squ / (deg_free);
    // cout << "resolution:" << res << endl;
    // cout << "goodness-of-fit:" << goodness << endl;
    // hist2->Draw();

    // myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    // myfile.close(); 
}


void sum_energy_new(){

    float energies[15] = {0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t energies_plot[15] = {0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60, 80, 100, 150, 200, 250};

    // double_t energies_plot[13] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    double_t mean_values[15] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_notcorr.txt");
    myfile_events << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events.close();

    for (int i = 0; i < 2; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_energy_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }
    
    
    for (int i = 2; i < 13; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_energy_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }

    for (int i = 13; i < 15; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_energy_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }

    TCanvas *c_fit = new TCanvas("gaussian_fit", "gaussian_fit");

    TGraph *hits = new TGraph(9, mean_values, energies_plot);

    hits->Fit("pol4");
    hits->Draw("APL*");
    // hits->GetXaxis()->SetLimits(0,400);
    TF1 *fit_pol = hits->GetFunction("pol4");
    // double p0 = fit_pol->GetParameter(0);
    // double p1 = fit_pol->GetParameter(1);
    // double p2 = fit_pol->GetParameter(2);
    // double p3 = fit_pol->GetParameter(3);
    b0 = fit_pol->GetParameter(0);
    b1 = fit_pol->GetParameter(1);
    b2 = fit_pol->GetParameter(2);
    b3 = fit_pol->GetParameter(3);
    b4 = fit_pol->GetParameter(4);

    //cout << p0 << p1<< p2 << p3 << endl;

    for (const auto& e : mean_values) {
        std::cout << e << " " << eesum(e) << std::endl;
        }

    double_t mean_values_new[11];

    //just summing the events, without applying masking
    ofstream myfile_events2;
    myfile_events2.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_corr.txt");
    myfile_events2 << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events2.close();

    for (int i = 0; i < 2; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], a0, a1, a2, a3);

    }
    
    
    for (int i = 2; i < 13; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], a0, a1, a2, a3);

    }

    for (int i = 13; i < 15; i++) {

        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "GeV.root";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        sum_energy_corr((char*)argument.c_str(), energies[i], a0, a1, a2, a3);

    }    

    
}




