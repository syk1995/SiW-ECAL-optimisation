
void sum_events_masking(char* filename, int input_energy){

    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_events_masking.txt", std::ios_base::app);

    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s.root", filename), "read");

    TTree *tree = (TTree*)input->Get("ecal");

    tree -> Draw("Sum$(hit_isMasked==0)");
    auto hist = (TH1F*)gPad->GetPrimitive("htemp");

    //TH1F *hist = new TH1F("hist", "Histogram", 5000, 0, 5000);

    //get the number of events
    int entries = tree->GetEntries();
    
    cout << "Number of events:" << entries << endl;

    vector<int> *hit_slab = 0;

    vector<int> *hit_isMasked = 0;

    tree->SetBranchAddress("hit_slab", &hit_slab);

    tree->SetBranchAddress("hit_isMasked", &hit_isMasked);
    
    //iterating over all events
    // for (int i = 0; i < entries; i++) {
    //     tree->GetEntry(i);
    //     ULong_t n_hit = hit_isMasked->size(); //number of hits per event
    //     //cout << "Number of hits:" << n_hit << endl;
    //     int *mask = hit_isMasked->data(); //if hit is masked or not
    //     //exclude all the hits which were masked 
    //     int n_hit_unmasked = 0;
    //     for (int j = 0; j < n_hit; j++) {
    //         if (*mask == 0)
    //         {
    //             n_hit_unmasked++;
                
    //         }
            
    //         //cout << i << " : " << j << " : " << *mask << endl;
    //         mask++;
    //     }
    //     hist->Fill(n_hit_unmasked);

    // }
     //hist->Fit("gaus", "C", "M", 130, 210);
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