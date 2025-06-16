
double f6(int layer){
    // sampling fraction first version -> sampling fraction 1

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}

double f7(int layer){
    // sampling fraction second version -> sampling fraction 2

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = 1/(d_SiX0[layer] + d_WX0[layer]);
	
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


void get_data_shower(){

    float energies[16] = {0.25, 0.5, 0.7, 1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150, 200, 250};

    //iterate over all the different energies
    for (int i = 0; i < 16; i++) {
       //write the number of hits and weighted energy sum to a txt file
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );

        ofstream eventfile;
        eventfile.open(TString::Format("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_data_shower_%sGeV.txt", beam_energy.c_str()));
        eventfile << "nhits,wnhits,wenergy1,wenergy2,nonwenergy" << "\n";

        // TFile *input = new TFile(TString::Format("/home/llr/ilc/ritzmann/work/TB2022-06/CONF6/build/ECAL_QGSP_BERT_conf6_e-_%sGeV_5kevt_-42_-42_build_masked.root", to_string(energies[i]).c_str()), "read");
        TFile *input = open_file(energies[i]);
        TTree *tree = (TTree*)input->Get("ecal");

        //number of events
        int entries = tree->GetEntries();

        vector<int> *hit_energy = 0;
        vector<int> *hit_slab = 0;

        tree->SetBranchAddress("hit_energy", &hit_energy);
        tree->SetBranchAddress("hit_slab", &hit_slab);


        for (int i = 0; i < entries; i++) {
            tree->GetEntry(i);
            ULong_t n_hit = hit_energy->size(); //number of hits per event
            int *energy = hit_energy->data(); //energy of each hit
            int *slice = hit_slab->data(); //layer number of each hit
            double tot_energy_nonw = 0;
            double tot_energy_w1 = 0;
            double tot_energy_w2 = 0;
            double wnhits = 0;
            for (int j = 0; j < n_hit; j++) {
                int layer = *slice;
                tot_energy_w1 = tot_energy_w1 + *energy * f6(layer);
                tot_energy_w2 = tot_energy_w2 + *energy * f7(layer);
                tot_energy_nonw = tot_energy_nonw + *energy;
                wnhits = wnhits + f7(layer);
                
                //else cout << *slice << endl;
                //tot_energy = tot_energy + *energy;
                
                //cout << i << " : " << j << " : " << *energy << " " << *slice << endl;
                energy++;
                slice++;
            }
            eventfile << n_hit << "," << wnhits << "," << tot_energy_w1 << "," << tot_energy_w2 << "," << tot_energy_nonw << "\n";
        }

        eventfile.close();  
    }

}
