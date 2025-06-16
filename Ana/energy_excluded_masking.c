double f5(int layer){

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}

float calculateSD(vector<float> data, int length, float *sum_value, float *mean_value, float *std_deviation) {
  float sum = 0.0, mean, standardDeviation = 0.0;
  int i;

  for(i = 0; i < length; ++i) {
    sum += data[i];
  }

  mean = sum / (float)length;

  for(i = 0; i < length; ++i) {
    standardDeviation += pow(data[i] - mean, 2);
  }

  cout << sum << " " << mean << endl;  
  *sum_value = sum;
  *mean_value = mean;
  *std_deviation = sqrt(standardDeviation / length);
  //return sum, mean, sqrt(standardDeviation / length);
  return 0;
}



void energy_excluded_masking(){
  // calculates the amount of energy which is excluded by applying masking

    //write into file the energy, the mean and the total sume of energy excluded by applying the masking 
    //the mean is taken over all events, so the mean energy excluded per event
    ofstream myfile_events_exmasking;
    myfile_events_exmasking.open ("/home/llr/ilc/ritzmann/work/energy_excluded_masking.txt");
    myfile_events_exmasking << "energy, exenergy_mask_sum, exenergy_mask_mean, exenergy_mask_std" << "\n";

    //different energies
    int energies[11] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150};

    for (int i = 0; i < 11; i++){
        //read in root file for each energy
        TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_%sGeV_5kevt_build_masked.root", to_string(energies[i]).c_str()), "read");

        TTree *tree = (TTree*)input->Get("ecal");

        //number of events
        int entries = tree->GetEntries();

        vector<int> *hit_isMasked = 0;
        vector<int> *hit_energy = 0;
        vector<int> *hit_slab = 0;

        tree->SetBranchAddress("hit_energy", &hit_energy);
        tree->SetBranchAddress("hit_slab", &hit_slab);
        tree->SetBranchAddress("hit_isMasked", &hit_isMasked);


        vector<float> sum_energy;
        for (int i = 0; i < entries; i++) {
            float tot_energy = 0;
            tree->GetEntry(i);
            ULong_t n_hit = hit_energy->size(); //number of hits per event
            int *energy = hit_energy->data(); //energy of each hit
            int *slice = hit_slab->data(); //layer number of each hit
            int *masking = hit_isMasked->data(); //if hit is masked or not
            for (int j = 0; j < n_hit; j++) {
                int layer = *slice;
                int masked_bool = *masking;

                if (masked_bool == 1){
                    tot_energy = tot_energy + *energy * f5(layer);
                }
    
                energy++;
                slice++;
                masking++;
            }
            sum_energy.push_back(tot_energy);
        }
        float sum1, mean1, std_dev;
        calculateSD(sum_energy, sum_energy.size(), &sum1, &mean1, &std_dev);
        cout << sum1 << " " << mean1 << endl;
        myfile_events_exmasking << (double)energies[i] << "," << sum1 << "," << mean1 << "," << std_dev << "\n";
    }
    myfile_events_exmasking.close();

    




}