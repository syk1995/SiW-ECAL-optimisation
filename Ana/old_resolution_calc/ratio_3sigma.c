
double calc_ratio(TH1F * hist, int entries){
    double_t mean = hist->GetMean();
    double_t sigma = hist->GetStdDev();
    
    int_t nbin_left = hist->GetXaxis()->FindBin(mean - 3*simga);
    int_t nbin_right = hist->GetXaxis()->FindBin(mean + 3*sigma) + 1;

    double_t ex_entries = 0;

    for(int i = 0; i<nbin_range_left; i++){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    return ex_entries / entries;
}


void ratio_3sigma.c(){

    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    ofstream myfile;
    myfile.open("/home/llr/ilc/ritzmann/work/res_csv_files/conf6_e-_GeV_5kevt_-42_-42_build_masked_ratio_3sigmas.txt")
    myfile_events << "energy,nhits,sumE,wsumE1,wsumE2" << "\n";

    for (int i = 0; i < 3; i++) {
        string argument1 = "ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        string argument2 = "keV.root";
        string beam_energy = to_string(1000*energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2;

        TFile *input = new TFile(TString::Format("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/%s", filename), "read");
        TTree *tree = (TTree*)input->Get("ecal");
        double_t entries = tree->GetEntries();

        TCanvas *c = new TCanvas(("c_" + to_string(input_energy)).c_str(), to_string(input_energy).c_str());

        tree->Draw("nhit_len"); 
        auto hist_nhits = (TH1F*)gPad->GetPrimitive("htemp");

        tree->Draw()
        
    }

}