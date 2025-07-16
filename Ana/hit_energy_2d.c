// makes a 2d histogram of hit energy vs. input energy of the simulation

void hit_energy_2d(){
    // TSpline3 *spline_fit_hits = sum_hits_test();

    // read in the two files and add one as a friend
    TFile *input_ecal = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_0-6GeV_uniform_106.root");
    TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform_106.root");
    
    TTree *tree = (TTree*)input_ecal->Get("ecal");
    TTree *tpart = (TTree*)input_particle->Get("MyLCTuple");

    tree->AddFriend("MyLCTuple", "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform.root");
    
    tree->Draw("hit_energy:MyLCTuple.mcene[0]", "", "colz");

    TH2F* hist = new TH2F("hist", "E(sim) vs. hit energy", 250, 0, 250, 100, 0, 6.6);

    Float_t mcene[50];
    tpart->SetBranchAddress("mcene", &mcene);

    vector<int> *hit_energy = 0;
    tree->SetBranchAddress("hit_energy", &hit_energy);

    int entries_p = tpart->GetEntries();
    int entries_d = tree->GetEntries();
    // cout << entries_p << endl;

    for (int i = 0; i < min(entries_d,entries_p); i++){
        tpart->GetEntry(i);
        tree->GetEntry(i);
        int n_hit = hit_energy->size();
        int *energy = hit_energy->data();
        for (int j = 0; j < n_hit; j++){
            hist->Fill(*energy,mcene[0]);
            energy++;
        }
    }
    TCanvas *c = new TCanvas("new", "new");
    hist->Draw("COLZ");
}



