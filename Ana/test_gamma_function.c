// this was a file to test the fit of the gamma function
void test_gamma_function(){

    TFile *input = new TFile("/data_ilc/flc/ritzmann/simulations/TB2022-06/CONF6/build/ECAL_QGSP_BERT_TB2022-06_CONF6_e-_700keV.root");
    TTree *tree = (TTree*)input->Get("ecal");
    double_t entries = tree->GetEntries();

    TCanvas *c = new TCanvas("250MeV", "250MeV");

    tree->Draw("Sum$(hit_energy)");

    TH1F *hist = (TH1F*)gPad->GetPrimitive("htemp");

    TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-TMath::Power([1]/[2],2)*((x-[1])/[1]-TMath::Log(x/[1])))", 0, 100); // "xmin" = 0, "xmax" = 100
    f1->SetParameters(250, 20, 20); // you MUST set non-zero initial values for parameters
    f1->SetParNames("A", "x0", "s0");
    hist->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
    double A = f1->GetParameter("A");
    double x0 = f1->GetParameter("x0");
    double s0 = f1->GetParameter("s0");
    gStyle->SetOptFit(11111);

    cout << x0 + TMath::Power(s0,2)/x0 << " " << hist->GetMean() << endl;
    cout << s0*TMath::Sqrt(1 + TMath::Power(s0,2)/TMath::Power(x0,2)) << " " << hist->GetStdDev() << endl;



    hist->Draw();
}