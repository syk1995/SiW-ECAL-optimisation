void continuous_nhits_energy(){
    // 0-6 GeV energy simulation
    // TFile *input_ecal = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_0-6GeV_uniform_106.root");
    // TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-6_uniform_106.root");
    
    // 0-17 GeV energy simulation
    TFile *input_ecal = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_0-17GeV_uniform.root");
    TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-17_uniform.root");
    
    TTree *tree = (TTree*)input_ecal->Get("ecal");
    TTree *tpart = (TTree*)input_particle->Get("MyLCTuple");

    // create 2d histogram with the correct binning, so that every integer corresponds to a bin center
    tree->Draw("nhit_len>>nhits");
    TH1F* h_nhits = (TH1F*)gPad->GetPrimitive("nhits");

    int xmin_hits = h_nhits->GetXaxis()->GetXmin();
    int xmax_hits = h_nhits->GetXaxis()->GetXmax();

    int nbins = xmax_hits - xmin_hits;

    Double_t xbins[nbins+1];
    for(int i = 0; i < nbins + 1; i++){
        xbins[i] = i-0.5;
    }

    TH2F* hist = new TH2F("hist", "E(sim) vs. nhits", nbins, xbins, 300, 0, 17);

    int nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    vector<int> *hit_slab = 0;
    tree->SetBranchAddress("hit_slab", &hit_slab);

    Float_t mcene[10000];
    tpart->SetBranchAddress("mcene", &mcene);

    int entries_p = tpart->GetEntries();
    int entries_d = tree->GetEntries();
    cout << entries_d << endl;

    for (int i = 0; i < entries_d; i++){
        tpart->GetEntry(i);
        tree->GetEntry(i);
        if (nhit_len!=0){
        hist->Fill(nhit_len,mcene[0]);
        }
    }

    TCanvas *correction = new TCanvas("correction", "nhits");
    hist->Draw();

    // determine realtion between E(sim) and nhtis by projecting the input energy bin onto the x-axis -> mean hit count
    
    // 298 bins because we have 300 energy bins and I wanted to exclude the first two
    double input_energy[298];
    double mean_hits[298];
    double std_hits[298];


    for (int i = 2; i < 298 + 2; i++){
        TH1D *slice = hist->ProjectionX((to_string(i) + "_px_nhits").c_str(), i, i);
        input_energy[i-2] = hist->GetYaxis()->GetBinCenter(i);
        mean_hits[i-2] = slice->GetMean();
        std_hits[i-2] = slice->GetStdDev();
        // cout << number_hits[i-2] << " " << mean_energy[i-2]<< " " << number_hits[0] << endl;
    }

    TGraphErrors *enehits = new TGraphErrors(298, mean_hits, input_energy, std_hits, 0);
    TCanvas *data = new TCanvas("data", "nhits vs. energy");
    enehits->Draw();

    // spline interpolation
    TSpline3 *spline_hits = new TSpline3("Spline Interpolation", mean_hits, input_energy, 298);
    spline_hits->SetLineWidth(2);
    spline_hits->SetLineColor(kRed);
    spline_hits->Draw("same");

    // adjust the binning according to the spline
    // int nhit_bins = hist->GetNbinsX();
    int nhit_bins = mean_hits[297] + 1;
    Double_t xbins_sp[nhit_bins + 1];
    int m = 0;
    for (int i = 0; i < nhit_bins+1; i++){
        double val = spline_hits->Eval((i-0.5));
        if (val > xbins_sp[i-1]){
            xbins_sp[i] = val;
            // cout << 1 << endl;
        }
        else {
            xbins_sp[i] = xbins_sp[i-1] + 0.03;
            // xbins_sp[i-1] = val;
            cout << 0 << " " << xbins_sp[i] << " " << xbins_sp[i-1] << endl;
        }
        
        // xbins_sp[i] = val;
        cout << i << " " << xbins_sp[i] << endl;
    }

    TCanvas *histo = new TCanvas("histo", "E(sim) vs. E(nhits)");

    // 
    TH2F* hist_spline = new TH2F("hist", "E(sim) vs. E(nhits)", nhit_bins, xbins_sp, 300, 0, 17);

    for (int i = 0; i < entries_d; i++){
        tpart->GetEntry(i);
        tree->GetEntry(i);
        if (nhit_len!=0 && nhit_len <= mean_hits[297] && i != 174){
            hist_spline->Fill(spline_hits->Eval(nhit_len),mcene[0]);
        }
    }
    hist_spline->Draw();

    // I only took energies up to 4 GeV because afterwards the histogram isn't as nice anymore because of the spline
    double energies[17] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    
    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_yproj.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss_energy_corrected_yproj.txt");
    gauss_file << "energy,mean,sigma" << "\n";

    ofstream mean_file_x;
    mean_file_x.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_xproj.txt");
    mean_file_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean_x;
    every_bin_mean_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_mean_xproj.txt");
    every_bin_mean_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss_x;
    every_bin_gauss_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_gauss_xproj.txt");
    every_bin_gauss_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean;
    every_bin_mean.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_mean_yproj.txt");
    every_bin_mean << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss;
    every_bin_gauss.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_gauss_yproj.txt");
    every_bin_gauss << "energy,mean,sigma" << "\n";

   

    for (int i = 0; i < 17; i++){
        // look for energy in the reconstructed energy bins and then project onto the input energies
        TCanvas *c2 = new TCanvas(("new" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str());
        int bin = hist_spline->GetXaxis()->FindBin(energies[i]);
        TH1D* histy = hist_spline->ProjectionY((to_string(energies[i]) + "_py").c_str(), bin, bin);
        // histy->Draw();

        mean_file << energies[i] << "," << histy->GetMean() << "," << histy->GetStdDev() << "\n";

        histy->Fit("gaus");
        TF1 *fit = histy->GetFunction("gaus");
        double mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        gauss_file << energies[i] << "," << mean << "," << sigma << "\n";

        histy->Draw();

        // look for energy in the input energies and then project onto the reconstructed energies
        int bin_y = hist_spline->GetYaxis()->FindBin(energies[i]);
        TH1D* histy_x = hist_spline->ProjectionX((to_string(energies[i]) + "px").c_str(), bin_y, bin_y);
        histy_x->Draw();

        mean_file_x << energies[i] << "," << histy_x->GetMean() << "," << histy_x->GetStdDev() << "\n";
    }

    // calculate resolution for every bin center starting from reconstructed energy and projecting onto input energies
    int num_xbins = int(ceil(hist_spline->GetNbinsX() / 3.));
    // int num_xbins = nhit_bins;
    cout << "number of xbins: " << num_xbins << endl;
    for (int i = 2; i < num_xbins; i++){
        double energy_bin = hist_spline->GetXaxis()->GetBinCenter(i);
        cout << "number of hits: " << i << " bin energy: " << energy_bin << endl;
        TCanvas *c2 = new TCanvas(("new" + to_string(energy_bin)).c_str(), to_string(energy_bin).c_str());
        TH1D* histy_bins = hist_spline->ProjectionY((to_string(energy_bin) + "_py").c_str(), i, i);
        // histy_bins->Draw();
        every_bin_mean << energy_bin << "," << histy_bins->GetMean() << "," << histy_bins->GetStdDev() << "\n";
        
        histy_bins->Fit("gaus");
        TF1* fit_gauss = histy_bins->GetFunction("gaus");
        histy_bins->Draw();
        // fit_gauss->Draw();
        double mean = fit_gauss->GetParameter(1);
        double sigma = fit_gauss->GetParameter(2);
        // cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        every_bin_gauss << energy_bin << "," << mean << "," << sigma << "\n";
    }

    // calculate resolution for every bin center starting from input energy and projecting onto reconstructed energy
    int num_ybins = int(ceil(hist_spline->GetNbinsY()));
    cout << "number of ybins: " << num_ybins << endl;
    for (int i = 2; i < num_ybins; i++){
        double energy_bin = hist_spline->GetYaxis()->GetBinCenter(i);
        if (energy_bin <= 8){
            // TCanvas *c2 = new TCanvas(("new" + to_string(energy_bin)).c_str(), to_string(energy_bin).c_str());
            TH1D* histy_bins = hist_spline->ProjectionX((to_string(energy_bin) + "_px ").c_str(), i, i);
            // histy_bins->Draw();
            every_bin_mean_x << energy_bin << "," << histy_bins->GetMean() << "," << histy_bins->GetStdDev() << "\n";
            
            histy_bins->Fit("gaus");
            TF1* fit_gauss = histy_bins->GetFunction("gaus");
            // histy_bins->Draw();
            double mean = fit_gauss->GetParameter(1);
            double sigma = fit_gauss->GetParameter(2);
            // // cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
            every_bin_gauss_x << energy_bin << "," << mean << "," << sigma << "\n";
        }
    }

    mean_file.close();
    mean_file_x.close();
    every_bin_mean.close();
    every_bin_mean_x.close();
    gauss_file.close();
    every_bin_gauss.close();
    every_bin_gauss_x.close();





}