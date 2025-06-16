void get_fit_params(){
    int energies[11] = {1, 2, 5, 8, 10, 20, 40, 60 , 80, 100, 150};

    //just summing the events, without applying masking
    ofstream myfile_events;
    myfile_events.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_events.txt");
    myfile_events << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_events.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //just summing the events, with applying masking
    ofstream myfile_events_mask;
    myfile_events_mask.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_events_masking.txt");
    myfile_events_mask << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_mask.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_events_masking.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //just summing the energy, without weights
    ofstream myfile_events_energy;
    myfile_events_energy.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy.txt");
    myfile_events_energy << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_energy.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights
    ofstream myfile_events_wenergy;
    myfile_events_wenergy.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted.txt");
    myfile_events_wenergy << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_wenergy.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy_weighted.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights applying masking
    ofstream myfile_events_wenergy_mask;
    myfile_events_wenergy_mask.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted_masking.txt");
    myfile_events_wenergy_mask << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_wenergy_mask.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy_weighted_masking.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights applying masking
    ofstream myfile_events_wenergy2;
    myfile_events_wenergy2.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted2.txt");
    myfile_events_wenergy2 << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_wenergy2.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy_weighted2.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

    //summing energy with weights applying masking
    ofstream myfile_events_wenergy2_mask;
    myfile_events_wenergy2_mask.open ("/home/llr/ilc/ritzmann/work/ECAL_QGSP_BERT_conf6_e-_GeV_5kevt_-42_-42_build_masked_params_sum_energy_weighted2_masking.txt");
    myfile_events_wenergy2_mask << "energy,mean,mean_err,sigma,sigma_err" << "\n";
    myfile_events_wenergy2_mask.close();

    for (int i = 0; i < 11; i++) {
        string argument1 = ".x work/root_macros/sum_energy_weighted2_masking.c((char*)\"ECAL_QGSP_BERT_conf6_e-_";
        string argument2 = "GeV_5kevt_-42_-42_build_masked\",";
        string beam_energy = to_string(energies[i]);
        string argument = argument1 + beam_energy + argument2 + beam_energy + ")";
        cout << argument << endl;
        gROOT->ProcessLine(argument.c_str());
    }

}