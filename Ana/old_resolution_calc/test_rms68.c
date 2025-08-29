// this was a file just to test the implementation of the rms method and to calculate the correction factor
// for the standard deviation

void rms68_(TH1F *h, double &mean, double &rms_val) {
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t maxbin = h->GetMaximumBin();
    Int_t imean = axis->FindBin(h->GetMean());
    Double_t entries = 0.68*h->GetEntries();
    Double_t w = h->GetBinContent(maxbin);
    Double_t x = h->GetBinCenter(maxbin);
    Double_t sumw = w;
    Double_t sumwx = w*x;
    Double_t sumwx2 = w*x*x;
    Int_t range;
    for (Int_t i = 1; i<nbins;i++) {
    if (i> 0) {
        w = h->GetBinContent(maxbin-i);
        x = h->GetBinCenter(maxbin-i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (i<= nbins) {
        w = h->GetBinContent(maxbin+i);
        x = h->GetBinCenter(maxbin+i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (sumw > entries) {
        range = i;
        break;
    }
    }
    cout << "maxbin: " << maxbin << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(maxbin-range);
    Double_t x2 = h->GetBinLowEdge(maxbin+range+1);

    TH1F* hist68 = new TH1F("h68", "h68", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist68->GetNbinsX(); j++){
        hist68->SetBinContent(j, h->GetBinContent(maxbin-range + j));
    }

    hist68->Draw();
    mean = hist68->GetMean();
    rms_val = hist68->GetStdDev(); //this would be the mean value of the central 90 percent
    // 
    // Double_t rms2 = TMath::Abs(sumwx2/sumw);
    // rms_val = TMath::Sqrt(rms2);
    printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}

void rms90_(TH1F *h, double &mean, double &rms_val);

void rms90_(TH1F *h, double &mean, double &rms_val) {
    double median, percentile;
    percentile = 0.5;
    h->ComputeIntegral();
    h->GetQuantiles(1, &median, &percentile);
    cout << "median: " << median << endl;
    TAxis *axis = h->GetXaxis();
    Int_t nbins = axis->GetNbins();
    Int_t imedian = axis->FindBin(median);
    Double_t entries = 0.68*h->GetEntries();
    Double_t w = h->GetBinContent(imedian);
    Double_t x = h->GetBinCenter(imedian);
    Double_t sumw = w;
    Double_t sumwx = w*x;
    Double_t sumwx2 = w*x*x;
    Int_t range;
    for (Int_t i = 1; i<nbins;i++) {
    if (i > 0) {
        w = h->GetBinContent(imedian-i);
        x = h->GetBinCenter(imedian-i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (i<= nbins) {
        w = h->GetBinContent(imedian+i);
        x = h->GetBinCenter(imedian+i);
        sumw += w;
        sumwx += w*x;
        sumwx2 += w*x*x;
    }
    if (sumw > entries) {
        range = i;
        break;
    }
    }
    cout << "median bin: " << imedian << " " << "range: " << range << endl;
    Double_t x1 = h->GetBinLowEdge(imedian-range);
    Double_t x2 = h->GetBinLowEdge(imedian+range+1);

    TH1F* hist90 = new TH1F("h90", "h90", 2*range, x1, x2);
    for(Int_t j = 1; j <= hist90->GetNbinsX(); j++){
        hist90->SetBinContent(j, h->GetBinContent(imedian-range + j));
    }

        TCanvas *c = new TCanvas("40GeV", "40GeV");


    hist90->Draw();
    mean = hist90->GetMean();
    rms_val = hist90->GetStdDev(); //this would be the mean value of the central 90 percent
    // 
    // Double_t rms2 = TMath::Abs(sumwx2/sumw);
    // rms_val = TMath::Sqrt(rms2);
    printf("RMS of central 90 percent = %g, RMS total = %g\n",rms_val,h->GetRMS());
}



void test_rms68(){
    double sum = 0;
    for(int i = 0; i < 10000; i++){
        TH1F* h_rand = new TH1F("h1", "Histogram from a Gaussian", 100, -5, 5);
        h_rand->FillRandom("gaus", 5000);

        h_rand->Draw();

        double_t mean = h_rand->GetMean();
        double_t std_dev = h_rand->GetStdDev();

        double mean_rms, rms_val;
        rms90_(h_rand, mean_rms, rms_val);

        sum += std_dev / rms_val;

        // cout << "mean: " << mean << " " << "stddev: " << std_dev << endl;
        // cout << "rms mean: " << mean_rms << " " << "rms_stddev: " << rms_val << endl;
        // cout << std_dev / rms_val << endl;
    }

    cout << "mean correction: " << sum / 10000 << endl;
}
