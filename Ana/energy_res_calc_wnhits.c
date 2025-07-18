#include "/home/llr/ilc/ritzmann/work/root_macros/sum_hits_test.c"
#include "/home/llr/ilc/ritzmann/work/root_macros/sum_hits_weighted_test.c"

// this is exactly the same as energy_res_calc_nhits.c but now for the weighted number of hits

TSpline3 *spline_fit_hits = sum_hits_weighted_test();

void energy_res_calc_wnhits(){

    // read in the two files and add one as a friend
    TFile *input_ecal = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_0-6GeV_uniform_106.root");
    TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform_106.root");
    
    TTree *tree = (TTree*)input_ecal->Get("ecal");
    TTree *tpart = (TTree*)input_particle->Get("MyLCTuple");

    tree->AddFriend("MyLCTuple", "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform_106.root");
    
    tree->Draw("nhit_len", "MyLCTuple.mcene[0]<0.5 && MyLCTuple.mcene[0]>0.3");

    // apply the correction to nhit_len and fill it into a 2d histogram

    tree->Draw("Sum$(f2(hit_slab))>>nhits");
    TH1F* h_nhits = (TH1F*)gPad->GetPrimitive("nhits");

    tpart->Draw("mcene[0]>>part_ene");
    TH1F* part_ene = (TH1F*)gPad->GetPrimitive("part_ene");

    double xmin_hits = h_nhits->GetXaxis()->GetXmin();
    double xmax_hits = h_nhits->GetXaxis()->GetXmax();

    cout << xmin_hits << " " << xmax_hits << endl;

    double xmin_ene = part_ene->GetXaxis()->GetXmin();
    double xmax_ene = part_ene->GetXaxis()->GetXmax();

    cout << xmin_ene << " " << xmax_ene << endl;

    double factor = 1.2;
    int nbins = ceil((xmax_hits-xmin_hits)/factor);
    cout << "nbins: " << nbins << endl;
    Double_t xbins[nbins + 1];
    xbins[0] = spline_fit_hits->Eval(xmin_hits);
        
    for (int i = 0; i < nbins+1; i++){
        double val = spline_fit_hits->Eval(xmin_hits+factor*(i+1));
        xbins[i+1] = val;
        cout << i << " " << xbins[i+1] << endl;
        }

    TH2F* hist = new TH2F("hist", "E(sim) vs. E(nhits)", nbins, xbins, 100, 0, 6);

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
        int n_hit = hit_slab->size(); //number of hits per event
        int *slice = hit_slab->data(); //layer number of each hit
        double wnhits = 0;
        for (int j = 0; j < n_hit; j++) {
            int layer = *slice;
            wnhits = wnhits + f2(layer);
            slice++;
        }
        hist->Fill(spline_fit_hits->Eval(wnhits), mcene[0]);
        // hist->Fill(nhit_len,mcene[0]);
    }
    TCanvas *c = new TCanvas("new", "new");
    hist->Draw("col");

    double energies[11] = {0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    
    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_weighted_corr_mean_corrected.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_weighted_corr_gauss_corrected.txt");
    gauss_file << "energy,mean,sigma" << "\n";

    for (int i = 0; i < 11; i++){
        TCanvas *c2 = new TCanvas(("new" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str());
        int bin = hist->GetXaxis()->FindBin(energies[i]);
        TH1D* histy = hist->ProjectionY((to_string(energies[i]) + "_py").c_str(), bin, bin);
        // histy->Draw();

        mean_file << energies[i] << "," << histy->GetMean() << "," << histy->GetStdDev() << "\n";

        histy->Fit("gaus");
        TF1 *fit = histy->GetFunction("gaus");
        double mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        gauss_file << energies[i] << "," << mean << "," << sigma << "\n";

        histy->Draw();
    }
}

