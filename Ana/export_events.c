void export_events(char* filename){

	TFile *f = new TFile(TString::Format("/home/svadhera/Documents/EP/ECAL_RES/rootfiles/%s.root", filename));
	TTree *ecal = (TTree*)f->Get("ecal");
	TCanvas *c1 = new TCanvas();
	ecal -> Draw("sum_energy:nhit_len");
	auto graph = (TGraph*)gPad->GetPrimitive("Graph");
	streambuf *coutbuf = cout.rdbuf();
	ofstream out(TString::Format("/home/svadhera/Documents/EP/ECAL_RES/csv_files/events_data_%s.txt", filename));
	cout.rdbuf(out.rdbuf());
	
	int nevents = 4990;
	
	cout << "nhit_len" << "," << "sum_energy" <<endl;
	for(int event=0; event < nevents; event++) {
		
		double x= 0.0;
		double y= 0.0;
		
		graph -> GetPoint(event, x, y);
		cout << x << "," << y << endl;
		
    			
    	}
    	
    	cout.rdbuf(coutbuf);
}
