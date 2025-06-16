
void sum_energy(char* filename, float input_energy){

    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy.txt", std::ios_base::app);

    //input for the files fabricio gave me
    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s.root", filename), "read");
    

    //open files for photons
    // ofstream myfile;
    // myfile.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_energy.txt", std::ios_base::app);

    // //input for the photons
    // TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s.root", filename), "read");

    TTree *tree = (TTree*)input->Get("ecal");

    //TH1F *hist = new TH1F("hist", "Histogram", 55000, 0, 55000);

    tree -> Draw("Sum$(hit_energy)");
    auto hist = (TH1F*)gPad->GetPrimitive("htemp");

    //number of events
    int entries = tree->GetEntries();
    
    cout << "Number of events:" << entries << endl;

    vector<int> *hit_energy = 0;
    vector<int> *hit_isMasked = 0;

    tree->SetBranchAddress("hit_energy", &hit_energy);

    //iterating over all events
    // for (int i = 0; i < entries; i++) {
    //     tree->GetEntry(i);
    //     ULong_t n_hit = hit_energy->size(); //number of hits per event
    //     int *energy = hit_energy->data(); //energy of each hit
    //     int *mask = hit_isMasked->data(); //if hit is masked or not
    //     int tot_energy = 0;

    //     for (int j = 0; j < n_hit; j++) {
    //         //respecting masking
    //         // if (*mask == 0)
    //         // {
    //         //     tot_energy = tot_energy + *energy;
            
    //         // //cout << i << " : " << j << " : " << *energy << endl;
    //         // }
    //         //not respecting masking
    //         tot_energy = tot_energy + *energy;

    //         mask++;
    //         energy++;
    //     }
    //     //cout << "Total Energy:" << tot_energy << endl;
    //     hist->Fill(tot_energy);
    //  }
     hist->Draw();
     //hist->Fit("gaus", "C", "M", 130, 210);
     hist->Fit("gaus");
     TF1 *fit = hist->GetFunction("gaus");
     double mean = fit->GetParameter(1);
     double sigma = fit->GetParameter(2);
     double mean_err = fit->GetParError(1);
     double sigma_err = fit->GetParError(2);
     double chi_squ = fit->GetChisquare();
     int bins = hist->GetNbinsX();
     double goodness = chi_squ / (bins - 3);
     cout << "resolution:" << sigma/mean << endl;
     cout << "goodness-of-fit:" << goodness << endl;
     hist->Draw();

    //  myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
    //  myfile.close(); 
}