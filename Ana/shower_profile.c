

double fitf(double *x, double *par){
        double val = par[0] * par[2] * pow(x[0] * par[2], par[1]-1) * exp(-x[0] * par[2]);
        return val;
     }

double f1(int layer){

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}     

TFile* open_file(float input_energy){
    // cout << "file with energy: " << input_energy << " openend!" << endl;
    string argument1, argument2, beam_energy;
    if(input_energy < 1){
        argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        argument2 = "keV.root"; 
        beam_energy = to_string(int(1000*input_energy));
    }

    else if((input_energy >= 1) && (input_energy <= 150)){
        argument1 = "/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_";
        argument2 = "GeV_5kevt_-42_-42_build_masked.root";
        beam_energy = to_string(int(input_energy));
    }

    else {
        argument1 = "/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_";
        argument2 = "GeV.root";
        beam_energy = to_string(int(input_energy));
    }

    string argument = argument1 + beam_energy + argument2;

    TFile *input = new TFile((char*)argument.c_str(), "read");
    // TTree *tree = (TTree*)input->Get("ecal");

    return input;
}

void shower_profile(){

    //write the shower leakage to a txt file
    ofstream eventfile;
    eventfile.open("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_shower_leakage.txt");
    eventfile << "energy,shower_leakage" << "\n";


    int energies[13] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    //iterate over all the different energies
    for (int i = 0; i < 13; i++) {
        TFile *input = open_file(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");

        float edges[16] = {0.6, 1.8, 3, 4.2, 5.4, 6.6, 7.8, 9., 10.4, 12., 13.6, 15.2, 16.8, 18.4, 20., 21.6};

        TH1F *hist = new TH1F("hist", "Histogram", 15, edges);

        hist->Draw();
        int entries = tree->GetEntries();
    
        cout << "Number of events:" << entries << endl;

        vector<int> *hit_energy = 0;
        vector<int> *hit_slab = 0;

        tree->SetBranchAddress("hit_energy", &hit_energy);
        tree->SetBranchAddress("hit_slab", &hit_slab);

        float centers[15] = {1.2,  2.4,  3.6,  4.8,  6. ,  7.2,  8.4,  9.6, 11.2, 12.8, 14.4, 16. , 17.6, 19.2, 20.8};

        for (int i = 0; i < entries; i++) {
            tree->GetEntry(i);
            ULong_t n_hit = hit_slab->size(); 
            int *energy = hit_energy->data(); 
            int *slice = hit_slab->data();
            double hit_energy = 0;
            int layer = 0;
            for (int j = 0; j < n_hit; j++) {
                hit_energy = 0;
                layer = *slice;
                hit_energy = *energy;
                hist->Fill(centers[layer], hit_energy);
                energy++;
                slice++;
         }
        }
        //hist->Draw();
        cout << "mean:" << hist->GetMean() << " maximum bin:" << hist->GetMaximumBin()<< endl;
     
        TF1 *func = new TF1("fit", fitf, 0, 10000, 3);
        func->SetParNames("E0", "a", "b");
        func->SetParameters(120016000,3.78,0.5);
        //hist->Fit("fit", "", "", 0, 21.6);

        TCanvas *c = new TCanvas();
        hist->SetTitle((to_string(energies[i]) + " GeV").c_str());
        hist->Fit("fit");
        hist->Draw();

        double full_integral = func->Integral(0, 10000);
        double part_integral = func->Integral(0, edges[15]);
        cout << full_integral << "  " << part_integral << endl;
        eventfile << energies[i] << "," << 1 - part_integral/full_integral << endl;
    
    }
    eventfile.close();

}