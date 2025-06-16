void continuous_nhits_energy_toy(){
    TFile *input = new TFile("/home/llr/ilc/ritzmann/ToySim.root");
    TFile *input_particle = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/LCTuple-MConly_0-17_uniform.root");
    
    TTree *tree = (TTree*)input->Get("tree");

    tree->Draw("nhit_len>>nhits");
    TH1F* h_nhits = (TH1F*)gPad->GetPrimitive("nhits");

    int xmin_hits = h_nhits->GetXaxis()->GetXmin();
    int xmax_hits = h_nhits->GetXaxis()->GetXmax();

    int nbins = xmax_hits - xmin_hits;

    Double_t xbins[nbins+1];
    for(int i = 0; i < nbins + 1; i++){
        xbins[i] = i-0.5;
    }

    TH2F* hist = new TH2F("hist", "E(sim) vs. nhits", nbins, xbins, 200, 0, 20);

    double_t nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    vector<int> *hit_slab = 0;
    tree->SetBranchAddress("hit_slab", &hit_slab);

    double_t e_i;
    tree->SetBranchAddress("e_i", &e_i);

    int entries_d = tree->GetEntries();
    cout << entries_d << endl;

    for (int i = 0; i < entries_d; i++){
        tree->GetEntry(i);
        if (nhit_len!=0){
        hist->Fill(nhit_len,e_i);
        }
    }
    TCanvas *correction = new TCanvas("correction", "nhits");
    hist->Draw();

    double input_energy[198];
    double mean_hits[198];
    double std_hits[198];


    for (int i = 1; i < 198 + 1; i++){
        TH1D *slice = hist->ProjectionX((to_string(i) + "_px_nhits").c_str(), i, i);
        input_energy[i-1] = hist->GetYaxis()->GetBinCenter(i);
        mean_hits[i-1] = slice->GetMean();
        std_hits[i-1] = slice->GetStdDev();
    }

    TGraphErrors *enehits = new TGraphErrors(198, mean_hits, input_energy, std_hits, 0);
    TCanvas *data = new TCanvas("data", "nhits vs. energy");
    enehits->Draw();

    TSpline3 *spline_hits = new TSpline3("Spline Interpolation", mean_hits, input_energy, 198);
    spline_hits->SetLineWidth(2);
    spline_hits->SetLineColor(kRed);
    spline_hits->Draw("same");

    // int nhit_bins = hist->GetNbinsX();
    int nhit_bins = int(ceil(mean_hits[196]));
    Double_t xbins_sp[nhit_bins + 1-36];
    for (int i = 36; i < nhit_bins+1; i++){
        double val = spline_hits->Eval((i-0.5));
        // cout << i << ": " << val << endl;
        if (val < xbins_sp[i-37]){
            cout << 0 << endl;
            val = xbins_sp[i-37] + 0.02;
        }
        xbins_sp[i-36] = val;
        cout << i-36 << " " << xbins_sp[i-36] << endl;
    }

    TCanvas *data1 = new TCanvas("data1", "E(nhits) vs. energy");


    TH2F* hist_spline = new TH2F("hist", "E(sim) vs. E(nhits)", nhit_bins-36, xbins_sp, 200, 0, 20);

    for (int i = 0; i < entries_d; i++){
        tree->GetEntry(i);
        if (nhit_len>36 && nhit_len < mean_hits[197] && e_i < input_energy[197]){
        hist_spline->Fill(spline_hits->Eval(nhit_len),e_i);
        }
    }
    hist_spline->Draw();


    TCanvas *data2 = new TCanvas("data2", "energy vs. E(nhits)");

    // rotated histogram for profile y calculation
    // TH2F* hist_spline_turned = new TH2F("hist_turned", "E(nhits) vs. E(sim)", 200, 0, 20, nhit_bins, xbins_sp);

    // for (int i = 0; i < entries_d; i++){
    //     // cout << i << endl;
    //     // tpart->GetEntry(i);
    //     tree->GetEntry(i);
    //     if (nhit_len>36 && nhit_len < mean_hits[197] && e_i < input_energy[197]){
    //     hist_spline_turned->Fill(e_i, spline_hits->Eval(nhit_len));
    //     }
    //     // hist->Fill(nhit_len,e_i[0]);
    //     // hist->Fill(e_i[0], nhit_len);
    // }
    // hist_spline_turned->Draw();


    double energies[17] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4};
    
    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_yproj_toy.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss_energy_corrected_yproj_toy.txt");
    gauss_file << "energy,mean,sigma" << "\n";

    ofstream mean_file_x;
    mean_file_x.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_xproj_toy.txt");
    mean_file_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean_x;
    every_bin_mean_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_mean_xproj_toy.txt");
    every_bin_mean_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss_x;
    every_bin_gauss_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_gauss_xproj_toy.txt");
    every_bin_gauss_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean;
    every_bin_mean.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_mean_yproj_toy.txt");
    every_bin_mean << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss;
    every_bin_gauss.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_energy_every_bin_gauss_yproj_toy.txt");
    every_bin_gauss << "energy,mean,sigma" << "\n";

    TCanvas *c2 = new TCanvas("new", "new");


    for (int i = 1; i < 17; i++){
        TCanvas *c2 = new TCanvas(("new" + to_string(energies[i])).c_str(), to_string(energies[i]).c_str());
        int bin = hist_spline->GetXaxis()->FindBin(energies[i]);
        TH1D* histy = hist_spline->ProjectionY((to_string(energies[i]) + "_py").c_str(), bin, bin);
        histy->Draw();

        mean_file << energies[i] << "," << histy->GetMean() << "," << histy->GetStdDev() << "\n";

        histy->Fit("gaus");
        TF1 *fit = histy->GetFunction("gaus");
        double mean = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        gauss_file << energies[i] << "," << mean << "," << sigma << "\n";

        histy->Draw();

        int bin_y = hist_spline->GetYaxis()->FindBin(energies[i]);
        TH1D* histy_x = hist_spline->ProjectionX((to_string(energies[i]) + "px").c_str(), bin_y, bin_y);

        mean_file_x << energies[i] << "," << histy_x->GetMean() << "," << histy_x->GetStdDev() << "\n";
    }

    int num_xbins = int(ceil(hist_spline->GetNbinsX() / 2.));
    cout << "number of xbins: " << num_xbins << endl;
    for (int i = 2; i < num_xbins; i++){
        double energy_bin = hist_spline->GetXaxis()->GetBinCenter(i);
        cout << "number of hits: " << i << " bin energy: " << energy_bin << endl;
        // TCanvas *c2 = new TCanvas(("new" + to_string(energy_bin)).c_str(), to_string(energy_bin).c_str());
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

    int num_ybins = int(ceil(hist_spline->GetNbinsY()));
    cout << "number of ybins: " << num_ybins << endl;
    for (int i = 2; i < num_ybins; i++){
        double energy_bin = hist_spline->GetYaxis()->GetBinCenter(i);
        if (energy_bin <= 11.){
            // TCanvas *c2 = new TCanvas(("new" + to_string(energy_bin)).c_str(), to_string(energy_bin).c_str());
            TH1D* histy_bins = hist_spline->ProjectionX((to_string(energy_bin) + "_py").c_str(), i, i);
            // histy_bins->Draw();
            every_bin_mean_x << energy_bin << "," << histy_bins->GetMean() << "," << histy_bins->GetStdDev() << "\n";
            
            histy_bins->Fit("gaus");
            TF1* fit_gauss = histy_bins->GetFunction("gaus");
            // histy_bins->Draw();
            double mean = fit_gauss->GetParameter(1);
            double sigma = fit_gauss->GetParameter(2);
            // cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
            
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