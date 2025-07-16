
void sum_events(char* filename, float input_energy){

    //write the fit parameters to a txt file
    //for the e- files from fabricio
    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_events.txt", std::ios_base::app);

    //input of the files fabricio gave me
    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s.root", filename), "read");

    
    //for photon files
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt", std::ios_base::app);

    //input for the photon files
    //TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");


    TTree *tree = (TTree*)input->Get("ecal");

    tree->Draw("nhit_len");
    auto hist = (TH1F*)gPad->GetPrimitive("htemp");

    //TH1F *hist = new TH1F("hist", "Histogram", 5000, 0, 5000);

    //get the number of events
    int entries = tree->GetEntries();
    
    cout << "Number of events:" << entries << endl;


    //iterating over all events
    // for (int i = 0; i < entries; i++) {
    //     tree->GetEntry(i);
    //     ULong_t n_hit = hit_isMasked->size();  //number of hits per event
    //     //cout << "Number of hits:" << n_hit << endl;
    //     hist->Fill(n_hit);
    // }
     //hist->Fit("gaus", "C", "M", 130, 210);

    // TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]/[2]),(x/[2]))*(TMath::Exp(-([1]/[2])))/TMath::Gamma((x/[2])+1.)", 0, 10); // "xmin" = 0, "xmax" = 10
    // //TF1 *f1 = new TF1("f1","[0]*TMath::Power(([1]),(std::floor(x)/1.0))*(TMath::Exp(-([1])))/TMath::Gamma((std::floor(x)/1.0)+1.)", 0, 10); // "xmin" = 0, "xmax" = 10
    // f1->SetParameters(1, 1, 1); // you MUST set non-zero initial values for parameters
    // f1->SetParNames("p0", "p1", "p2");
    // hist->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
    // double p0 = f1->GetParameter("p0");
    // double p1 = f1->GetParameter("p1");
    // double p2 = f1->GetParameter("p2");
    // cout << sqrt(p1/p2) / (p1/p2) << endl;
    // cout << p2 << endl;
    // f1->SetParameters(1, 1);
    // f1->SetParNames("p0", "p1");
    // hist->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
    // double p0 = f1->GetParameter("p0");
    // double p1 = f1->GetParameter("p1");
    // cout << sqrt(p1) / p1 << endl;
    hist->Fit("gaus");
    TF1 *fit = hist->GetFunction("gaus");
    double mean = fit->GetParameter(1);
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
}


