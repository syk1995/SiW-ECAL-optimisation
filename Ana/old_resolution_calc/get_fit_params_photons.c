void get_fit_params_photons(){
    float energies[11] = {0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1, 4, 5, 10, 20};

    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events.txt");
    myfile_events << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_events.c((char*)\"ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_";
        string argument2 = "GeV\",";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //just summing the events, with applying masking
    // ofstream myfile_events_mask;
    // myfile_events_mask.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_events_masking.txt");
    // myfile_events_mask << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    // myfile_events_mask.close();

    // for (int i = 0; i < 4; i++) {
    //     string argument1 = ".x work/root_macros/sum_events_masking.c((char*)\"ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_";
    //     string argument2 = "GeV\",";
    //     string beam_energy = to_string(energies[i]);
    //     string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
    //     cout << argument << endl;
    //     gROOT->ProcessLine(argument.c_str());
    // }

    //just summing the energy, without weights
    ofstream myfile_events_energy;
    myfile_events_energy.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_energy.txt");
    myfile_events_energy << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_energy.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy.c((char*)\"ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_";
        string argument2 = "GeV\",";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights
    ofstream myfile_events_wenergy;
    myfile_events_wenergy.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_energy_weighted.txt");
    myfile_events_wenergy << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_wenergy.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy_weighted.c((char*)\"ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_";
        string argument2 = "GeV\",";
        string beam_energy = to_string(energies[i]);
        beam_energy.erase (beam_energy.find_last_not_of('0') + 1, std::string::npos);
        beam_energy.erase ( beam_energy.find_last_not_of('.') + 1, std::string::npos );
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights applying masking
    // ofstream myfile_events_wenergy_mask;
    // myfile_events_wenergy.open ("/home/llr/ilc/ritzmann/work/photon_analysis/ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_GeV_params_sum_energy_weighted_masking.txt");
    // myfile_events_wenergy << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    // myfile_events_wenergy.close();

    // for (int i = 0; i < 4; i++) {
    //     string argument1 = ".x work/root_macros/sum_energy_weighted_masking.c((char*)\"ECAL_QGSP_BERT_TB2022-06_CONF6_gamma_";
    //     string argument2 = "GeV\",";
    //     string beam_energy = to_string(energies[i]);
    //     string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
    //     cout << argument << endl;
    //     gROOT->ProcessLine(argument.c_str());
    // }

}