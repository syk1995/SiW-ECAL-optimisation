void continuous_nhits_nhits_toy(){
    TFile *input = new TFile("/home/llr/ilc/ritzmann/ToySim.root");
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

    TH2F* hist = new TH2F("hist", "E(sim) vs. nhits", nbins, xbins, 300, 0, 20);

    double_t nhit_len = 0;
    tree->SetBranchAddress("nhit_len", &nhit_len);

    double_t e_i;
    tree->SetBranchAddress("e_i", &e_i);

    int entries_d = tree->GetEntries();
    cout << entries_d << endl;

    for (int i = 0; i < entries_d; i++){
        // cout << i << endl;
        tree->GetEntry(i);
        if (nhit_len!=0){
        hist->Fill(nhit_len,e_i);
        }
    }

    TCanvas *correction = new TCanvas("correction", "nhits");
    hist->Draw();

    double number_hits[350];
    double mean_energy[350];
    double std_energy[350];


    for (int i = 2; i < 350 + 2; i++){
        TH1D *slice = hist->ProjectionY((to_string(i) + "_py_nhits").c_str(), i, i);
        number_hits[i-2] = i-1;
        mean_energy[i-2] = slice->GetMean();
        std_energy[i-2] = slice->GetStdDev();
        // cout << number_hits[i-2] << " " << mean_energy[i-2]<< " " << number_hits[0] << endl;
    }


    TGraphErrors *enehits = new TGraphErrors(350, number_hits, mean_energy, 0, std_energy);
    TCanvas *data = new TCanvas("data", "nhits vs. energy");
    enehits->Draw();

    TSpline3 *spline_hits = new TSpline3("Spline Interpolation", number_hits, mean_energy, 350);
    spline_hits->SetLineWidth(2);
    spline_hits->SetLineColor(kRed);
    spline_hits->Draw("same");


    // I excluded the first 15 bins, because there the spline always gives back the same values 
    Double_t xbins_sp[350 + 1 - 15];
    for (int i = 15; i < 350+1; i++){
        double val = spline_hits->Eval((i-0.5));
        // if (val > xbins[i-1]){
        //     xbins[i] = val;
        //     cout << xbins[i] << endl;
        // }
        // else {
        //     xbins[i] = xbins[i-1];
        //     xbins[i-1] = val;
        //     cout << 0 << " " << xbins[i]<< endl;
        // }
        xbins_sp[i-15] = val;
        cout << i << " " << xbins_sp[i-15] << endl;
    }

    TCanvas *data2 = new TCanvas("data2", "nhits vs. energy");
    
    TH2F* hist_spline = new TH2F("hist", "E(sim) vs. E(nhits)", 350 - 15, xbins_sp, 300, 0, 20);


    for (int i = 0; i < entries_d; i++){
        tree->GetEntry(i);
        if (nhit_len!=0 && nhit_len < 350 ){
        hist_spline->Fill(spline_hits->Eval(nhit_len),e_i);
        }
    }
    hist_spline->Draw();


    double energies[26] = {0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5};
    
    ofstream mean_file;
    mean_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_yproj_toy.txt");
    mean_file << "energy,mean,sigma" << "\n";

    ofstream gauss_file;
    gauss_file.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_gauss_nhits_yproj_gauss_toy.txt");
    gauss_file << "energy,mean,sigma" << "\n";

    ofstream mean_file_x;
    mean_file_x.open ("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_xproj_toy.txt");
    mean_file_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean_x;
    every_bin_mean_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_every_bin_mean_xproj_toy.txt");
    every_bin_mean_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss_x;
    every_bin_gauss_x.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_every_bin_gauss_xproj_toy.txt");
    every_bin_gauss_x << "energy,mean,sigma" << "\n";

    ofstream every_bin_mean;
    every_bin_mean.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_every_bin_mean_yproj_toy.txt");
    every_bin_mean << "energy,mean,sigma" << "\n";

    ofstream every_bin_gauss;
    every_bin_gauss.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_hits_corr_mean_nhits_every_bin_gauss_yproj_toy.txt");
    every_bin_gauss << "energy,mean,sigma" << "\n";

    

    for (int i = 0; i < 26; i++){
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

        int bin_y = hist_spline->GetYaxis()->FindBin(energies[i]);
        TH1D* histy_x = hist_spline->ProjectionX((to_string(energies[i]) + "px").c_str(), bin_y, bin_y);

        mean_file_x << energies[i] << "," << histy_x->GetMean() << "," << histy_x->GetStdDev() << "\n";
    }

    int num_xbins = hist_spline->GetNbinsX();
    cout << "number of xbins: " << num_xbins << endl;
    for (int i = 15; i < num_xbins; i++){
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
        // // cout << "energy " << energies[i] << " GeV goodness-of-fit: " << fit->GetChisquare() / fit->GetNDF() << endl;
        every_bin_gauss << energy_bin << "," << mean << "," << sigma << "\n";
    }

    int num_ybins = hist_spline->GetNbinsY();
    cout << "number of ybins: " << num_ybins << endl;
    for (int i = 2; i < num_ybins; i++){
        double energy_bin = hist_spline->GetYaxis()->GetBinCenter(i);
        if (energy_bin <= 11.){
            cout << "number of hits: " << i << " bin energy: " << energy_bin << endl;
            // TCanvas *c2 = new TCanvas(("new" + to_string(energy_bin)).c_str(), to_string(energy_bin).c_str());
            TH1D* histy_bins = hist_spline->ProjectionX((to_string(energy_bin) + "_px").c_str(), i, i);
            // histy_bins->Draw();
            every_bin_mean_x << energy_bin << "," << histy_bins->GetMean() << "," << histy_bins->GetStdDev() << "\n";
        
            histy_bins->Fit("gaus");
            TF1* fit_gauss = histy_bins->GetFunction("gaus");
            // histy_bins->Draw();
            // fit_gauss->Draw();
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