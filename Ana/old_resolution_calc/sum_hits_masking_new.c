
double_t sum_hits_masking_notcorr(char* filename, int input_energy){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_masking_notcorr.txt", std::ios_base::app);

    //input of the files fabricio gave me
    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

    
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");


    TTree *tree = (TTree*)input->Get("ecal");

    tree->Draw("Sum$(hit_isMasked==0)"); 
    auto hist = (TH1F*)gPad->GetPrimitive("htemp");

    //get the number of events
    int entries = tree->GetEntries();
    
    cout << "Number of events:" << entries << endl;

    hist->Fit("gaus");
    TF1 *fit = hist->GetFunction("gaus");
    double_t mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double mean_err = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    float res = sigma/mean;
    double chi_squ = fit->GetChisquare();
    int bins = hist->GetNbinsX();
    double goodness = chi_squ / (bins - 3);
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    hist->Draw();

    myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    myfile.close(); 

    return mean;
}


double d0 = 0;
double d1 = 0;
double d2 = 0;
double d3 = 0;
double d4 = 0;

double enhits_mask(int number_hits){
    
    double val = d4*pow(number_hits,4)+d3*pow(number_hits,3)+d2*pow(number_hits,2)+d1*number_hits+d0;
    return val;
}

double_t sum_hits_masking_notcorr(char* filename, int input_energy);




void sum_hits_masking_corr(char* filename, float input_energy, double a0, double a1, double a2, double a3){
    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_masking_corr.txt", std::ios_base::app);

    //input of the files fabricio gave me
    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s", filename), "read");

    
    // //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    // //input for the photon files
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");



    TTree *tree = (TTree*)input->Get("ecal");


    tree->Draw("enhits_mask(Sum$(hit_isMasked==0))"); 
    auto hist2 = (TH1F*)gPad->GetPrimitive("htemp");

    //get the number of events
    int entries = tree->GetEntries();
    
    cout << "Number of events:" << entries << endl;

    hist2->Fit("gaus");
    TF1 *fit2 = hist2->GetFunction("gaus");
    double_t mean = fit2->GetParameter(1);
    double sigma = fit2->GetParameter(2);
    double mean_err = fit2->GetParError(1);
    double sigma_err = fit2->GetParError(2);
    float res = sigma/mean;
    double chi_squ = fit2->GetChisquare();
    int bins = hist2->GetNbinsX();
    double goodness = chi_squ / (bins - 3);
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    hist2->Draw();

    myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    myfile.close(); 
}


void sum_hits_masking_new(){

    int energies[11] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150};

    double_t energies_plot[11] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150};

    double_t mean_values[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_masking_notcorr.txt");
    myfile_events << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events.close();

    for (int i = 0; i < 11; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2;

        double_t mean_val = sum_hits_masking_notcorr((char*)argument.c_str(), energies[i]);

        mean_values[i] = mean_val;

    }

    TGraph *hits = new TGraph(11, mean_values, energies_plot);

    hits->Fit("pol4");
    hits->Draw("APL*");
    // hits->GetXaxis()->SetLimits(0,400);
    TF1 *fit_pol = hits->GetFunction("pol4");
    // double p0 = fit_pol->GetParameter(0);
    // double p1 = fit_pol->GetParameter(1);
    // double p2 = fit_pol->GetParameter(2);
    // double p3 = fit_pol->GetParameter(3);
    d0 = fit_pol->GetParameter(0);
    d1 = fit_pol->GetParameter(1);
    d2 = fit_pol->GetParameter(2);
    d3 = fit_pol->GetParameter(3);
    d4 = fit_pol->GetParameter(4);

    //cout << p0 << p1<< p2 << p3 << endl;


    double_t mean_values_new[11];

    //just summing the events, without applying masking
    ofstream myfile_events2;
    myfile_events2.open ("/home/llr/ilc/ritzmann/work/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_masking_corr.txt");
    myfile_events2 << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events2.close();

    for (int i = 0; i < 11; i++) {

        string argument1 = "ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2;

        sum_hits_masking_corr((char*)argument.c_str(), energies[i], d0, d1, d2, d3);

    }

    
}




