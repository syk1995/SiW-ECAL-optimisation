// double f(int layer){
//     // X0 stands for the radiation length in the corresponding material

// 	double X0W = 3.504;
// 	double X0Si = 93.70*1000; //microns!
// 	double Ecal_WThickness = double(0.7)/X0W;
// 	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
// 	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
// 	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
	
// 	return 1.0/f0;
// }



void sum_energy_weighted_masking(char* filename, int input_energy){

    ofstream myfile;
    myfile.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_masking.txt", std::ios_base::app);

    TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/%s.root", filename), "read");

    TTree *tree = (TTree*)input->Get("ecal");

    tree -> Draw("Sum$((hit_isMasked==0)*hit_energy*f(hit_slab))");
    auto hist = (TH1F*)gPad->GetPrimitive("htemp");
    // cout << htemp->GetNbinsX() << endl;

    //TH1F *hist = new TH1F("hist", "Histogram", 1000, 0, 3400000);

    //number of events
    int entries = tree->GetEntries();
    
    cout << "Number of event:" << entries << endl;

    vector<int> *hit_energy = 0;
    vector<int> *hit_slab = 0;

    tree->SetBranchAddress("hit_energy", &hit_energy);
    tree->SetBranchAddress("hit_slab", &hit_slab);

    //iterating over all events
    // for (int i = 0; i < entries; i++) {
    //     tree->GetEntry(i);
    //     ULong_t n_hit = hit_energy->size(); //number of hits per event
    //     int *energy = hit_energy->data(); //energy of each hit
    //     int *slice = hit_slab->data(); //layer number of each hit
    //     double tot_energy = 0;
    //     for (int j = 0; j < n_hit; j++) {
    //         int layer = *slice;
    //         // if (layer < 8)
    //         // {    
    //         //     tot_energy = tot_energy + *energy * 4.2;
    //         // }

    //         // else if (layer >= 8) 
    //         // {
    //         //     tot_energy = tot_energy + *energy * 5.6;
    //         // } 
    //         //cout << "layer:" << layer << "  energy:" << *energy << "  sampling func:" << f2(layer) << endl;
    //         tot_energy = tot_energy + *energy * f2(layer);
            
    //         //else cout << *slice << endl;
    //         //tot_energy = tot_energy + *energy;
            
    //         //cout << i << " : " << j << " : " << *energy << " " << *slice << endl;
    //         energy++;
    //         slice++;
    //     }
    //     hist->Fill(tot_energy);
    //  }
     gStyle->SetOptFit(kTRUE);
     hist->Draw();
     //hist->Fit("gaus", "", "", 2200, 3400);
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

     myfile << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "\n";
     myfile.close(); 
 }