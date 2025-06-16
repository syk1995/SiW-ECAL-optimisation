// all the functions used for the energy resolution caculation

double calc_ratio(TH1F * hist, int entries);
double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries);

void gamma_fit(TH1F *h, float input_energy, ofstream &gamma_file);
void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries);
void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy);

void rms90(TH1F *h, double &mean, double &rms_val);
void rms68(TH1F *h, double &mean, double &rms_val);

// functions for the weighting
double f(int layer);
double f2(int layer);



double f(int layer){
    // X0 stands for the radiation length in the corresponding material

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
    // double f0 = 1/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}

double f2(int layer){
    // X0 stands for the radiation length in the corresponding material

	double X0W = 3.504;
	double X0Si = 93.70*1000; //microns!
	double Ecal_WThickness = double(0.7)/X0W;
	double d_SiX0[15] = {650.0/X0Si, 650.0/X0Si, 650.0/X0Si, 650.0/X0Si,500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si, 500.0/X0Si,320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si, 320.0/X0Si};
	double d_WX0[15] = {6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,6*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness,8*Ecal_WThickness};
	
	// double f0 = (d_SiX0[layer])/(d_SiX0[layer] + d_WX0[layer]);
    double f0 = 1/(d_SiX0[layer] + d_WX0[layer]);
	
	return 1.0/f0;
}


void gamma_fit(TH1F *h, float input_energy, ofstream &gamma_file){
    // fit gamma function to histogram

    TCanvas *c = new TCanvas(("gamma" + to_string(input_energy)).c_str(), ("gamma" + to_string(input_energy)).c_str());

    TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-TMath::Power([1]/[2],2)*((x-[1])/[1]-TMath::Log(x/[1])))", 0, h->GetBinLowEdge(h->GetNbinsX()) + 20); // "xmin" = 0, "xmax" = 100
    f1->SetParameters(h->GetMaximum(), h->GetMean(), h->GetStdDev()); // you MUST set non-zero initial values for parameters
    f1->SetParNames("A", "x0", "s0");
    h->Fit("f1", "R"); // "R" = fit between "xmin" and "xmax" of the "f1"
    double A = f1->GetParameter("A");
    double x0 = f1->GetParameter("x0");
    double s0 = f1->GetParameter("s0");
    gStyle->SetOptFit(11111);
    double chisq = f1->GetChisquare();
    double ndf = f1->GetNDF();

    gamma_file << input_energy << "," << x0 + TMath::Power(s0,2)/x0 << "," << s0*TMath::Sqrt(1 + TMath::Power(s0,2)/TMath::Power(x0,2)) << "," << chisq << "," << ndf << "\n";
    // cout << x0 + TMath::Power(s0,2)/x0 << " " << h->GetMean() << endl;
    // cout << s0*TMath::Sqrt(1 + TMath::Power(s0,2)/TMath::Power(x0,2)) << " " << h->GetStdDev() << endl;
    h->Draw();   
}

void fit_gauss(TH1F *h, ofstream &gauss_file, float input_energy, ofstream &myfile_3sig_fit, double entries){
    // fit gaussian to histogram 

    h->Fit("gaus");
    TF1 *fit = h->GetFunction("gaus");
    double_t mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double mean_err = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    float res = sigma/mean;
    double goodness = fit->GetChisquare() / fit->GetNDF();
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    h->Draw();
    gStyle->SetOptFit(11111);
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    gauss_file << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";

    double ratio = calc_ratio_fitted(h, mean, sigma, entries);
    myfile_3sig_fit << input_energy << "," << ratio << "\n";
}

void fit_gauss_2sig(TH1F *h, ofstream &gauss_file2s, double mean_hist, double std_hist, float input_energy){
    // fit gaussian within two standard deviations of the mean of the histogram

    h->Fit("gaus", "", "", mean_hist - 2 * std_hist, mean_hist + 2 * std_hist);
    TF1 *fit = h->GetFunction("gaus");
    double_t mean = fit->GetParameter(1);
    double sigma = fit->GetParameter(2);
    double mean_err = fit->GetParError(1);
    double sigma_err = fit->GetParError(2);
    float res = sigma/mean;
    double chisq = fit->GetChisquare();
    double ndf = fit->GetNDF();
    double goodness = chisq / ndf;
    cout << "resolution:" << res << endl;
    cout << "goodness-of-fit:" << goodness << endl;
    h->Draw();
    gStyle->SetOptFit(11111);
    gSystem->ProcessEvents();
    gSystem->ProcessEvents();
    gauss_file2s << input_energy << "," << mean << "," << mean_err << "," << sigma << "," << sigma_err << "," << goodness << "\n";
}

double calc_ratio(TH1F * hist, int entries){
    // calculate the number of events outside of three sigmas divided by the total number of events
    // sigma just the standard deviation of the mean
    double_t mean = hist->GetMean();
    double_t sigma = hist->GetStdDev();
    
    int nbin_left = hist->GetXaxis()->FindBin(mean - 3*sigma);
    int nbin_right = hist->GetXaxis()->FindBin(mean + 3*sigma) + 1;

    double_t ex_entries = 0;

    for(int i = 0; i<nbin_left; i++){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    return ex_entries / entries;
}

double calc_ratio_fitted(TH1F *hist, double mean_fit, double sigma_fit, double entries){
    // calculate the number of events outside of three sigmas divided by the total number of events
    // sigma determined by gaussian fit in function fit_gauss
    int nbin_left = hist->GetXaxis()->FindBin(mean_fit - 3*sigma_fit);
    int nbin_right = hist->GetXaxis()->FindBin(mean_fit + 3*sigma_fit) + 1;

    double_t ex_entries = 0;

    for(int i = 0; i<nbin_left; i++){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    for(int i = nbin_right; i <= hist->GetNbinsX(); i++ ){
        ex_entries = ex_entries + hist->GetBinContent(i);
    }

    return ex_entries / entries;
}

void rms90(TH1F *h, double &mean, double &rms_val) {
    // calculate the mean and the standard deviation of the central 90 percent

    double median, percentile;
    percentile = 0.5;
    h->ComputeIntegral();
    h->GetQuantiles(1, &median, &percentile);
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t imedian = axis->FindBin(median);
    Double_t entries = 0.9*h->GetEntries();
    Double_t w = h->GetBinContent(imedian); //content of central bin
    Double_t x = h->GetBinCenter(imedian); //center of central bin
    Double_t sumw = w; //sum of all the entries
    // Double_t sumwx = w*x; //sum of all the entry values
    // Double_t sumwx2 = w*x*x; 
    Int_t range;
    // go out from the central bin until 90% of the entries are included
    for (Int_t i = 1; i<nbins;i++) {
        if (i> 0) {
            w = h->GetBinContent(imedian-i);
            x = h->GetBinCenter(imedian-i);
            sumw += w;
            // sumwx += w*x;
            // sumwx2 += w*x*x;
        }
        if (i<= nbins) {
            w = h->GetBinContent(imedian+i);
            x = h->GetBinCenter(imedian+i);
            sumw += w;
            // sumwx += w*x;
            // sumwx2 += w*x*x;
        }
        if (sumw > entries) {
            range = i;
            break;
        }
    }
    // cout << "median bin: " << imedian << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(imedian-range);
    Double_t x2 = h->GetBinLowEdge(imedian+range+1);

    // make new histogram containing central 90%
    TH1F* hist90 = new TH1F("h90", "h90", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist90->GetNbinsX(); j++){
        hist90->SetBinContent(j, h->GetBinContent(imedian-range + j));
    }
    // hist90->Draw();
    mean = hist90->GetMean();
    rms_val = hist90->GetStdDev(); 
    // printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}

void rms68(TH1F *h, double &mean, double &rms_val) {
    // calculate the mean and the standard deviation of the central 68 percent
    double median, percentile;
    percentile = 0.5;
    h->ComputeIntegral();
    h->GetQuantiles(1, &median, &percentile);
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t imedian = axis->FindBin(median);
    Double_t entries = 0.68*h->GetEntries();
    Double_t w = h->GetBinContent(imedian);
    Double_t x = h->GetBinCenter(imedian);
    Double_t sumw = w;
    // go out from the central bin until 68% of the entries are included
    Int_t range;
    for (Int_t i = 1; i<nbins;i++) {
        if (i> 0) {
            w = h->GetBinContent(imedian-i);
            x = h->GetBinCenter(imedian-i);
            sumw += w;
        }
        if (i<= nbins) {
            w = h->GetBinContent(imedian+i);
            x = h->GetBinCenter(imedian+i);
            sumw += w;
        }
        if (sumw > entries) {
            range = i;
            break;
        }
    }
    // cout << "median bin: " << imedian << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(imedian-range);
    Double_t x2 = h->GetBinLowEdge(imedian+range+1);

    // new histogram with central 68 percent
    TH1F* hist68 = new TH1F("h68", "h68", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist68->GetNbinsX(); j++){
        hist68->SetBinContent(j, h->GetBinContent(imedian-range + j));
    }

    // hist68->Draw();
    mean = hist68->GetMean();
    rms_val = hist68->GetStdDev(); 
    // printf("RMS of central 68 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}
