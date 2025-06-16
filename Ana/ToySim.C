/// \file
/// \ingroup tutorial_tree
/// \notebook -nodraw
/// Simple Detector Toy model
///
/// Code adapted from tree0.C Event class example
///
/// execute as: .x ToySim.C++
///
///
/// ### Effect of ClassDef() and ClassImp() macros
///
/// After running this macro create an instance of Det and Event
///
/// ~~~
///   Det d;
///   Event e;
/// ~~~
///
/// now you can see the effect of the  ClassDef() and ClassImp() macros.
/// (for the Det class these commands are commented!)
/// For instance 'e' now knows who it is:
///
/// ~~~
///   cout<<e.Class_Name()<<endl;
/// ~~~
///
/// whereas d does not.
///
/// The methods that are added by the ClassDef()/Imp() macro can be listed with
///
/// ~~~
/// .class
///   .class Event
///   .class Det
/// ~~~
///
/// \macro_code
///
/// \author Heiko.Scheit@mpi-hd.mpg.de

#include <TRandom.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>

#include <Riostream.h>

// class Event { //TObject is not required by this example
class Event : public TObject
{
public:
    Double_t e_i; // input energy
    Double_t sum_energy; // output energy
    Double_t nhit_len; // output hits
    ClassDef(Event, 1)
};

double eLin(double x){
    double eLin = 0.2 + 1.01*x + 0.0002*x*x; // small quadratic non-linarity
    return eLin;
}

double eRes(double x){
    if( x <0 ) return 0.0;
    double sigma = 0.2/sqrt(x)+0.01;
    double eRes= x*(1 + gRandom->Gaus(0., sigma)); 
    eRes = max(0.0, eRes);
    return eRes;
}

double hLin(double x){
    double hLin = max(0.0, 100*sqrt(x)); // strong quadratic non-linarity
    return hLin;
}

double hRes(double x){
    if( x <0 ) return 0.0;
    double sigma = 1.0/sqrt(x); // Mettre une poisson
    double hRes= x*(1 + gRandom->Gaus(0., sigma)); 
    hRes = max(0.0, hRes);
    return hRes;
}

ClassImp(Event)

void ToySim()
{
    // create a file
    TFile * f = new TFile("ToySim.root", "RECREATE");
    // create a TTree
    TTree *tree = new TTree("tree", "treelibrated tree");
    Event *e = new Event;

    // create a branch with energy
    tree->Branch("event", &e);

    // fill some events with random numbers
    Int_t nevent = 1000000;

    double_t Emin =  0.0;
    double_t Emax = 20.0;

    for (Int_t iev = 0; iev < nevent; iev++)
    {
        if (iev % 1000 == 0)
            cout << "Processing event " << iev << "..." << endl;

        Float_t e_i;
        e_i = gRandom->Rndm()*(Emax-Emin)+Emin;
        e->e_i = e_i;
        e->sum_energy = eRes(eLin(e_i));
        e->nhit_len = hRes(hLin(e_i));

        tree->Fill(); // fill the tree with the current event
    }

    // start the viewer
    // here you can investigate the structure of your Event class
    tree->StartViewer();

    // gROOT->SetStyle("Plain");   // uncomment to set a different style

    // now draw some tree variables
    TCanvas *c1 = new TCanvas();
    c1->Divide(1, 3);
    c1->cd(1);
    tree->Draw("e_i");                                  // input energy
    c1->cd(2);
    tree->Draw("sum_energy:e_i>>hE(200,0,20, 200, 0, 20)", "", "colz"); // one energy against the other
    c1->cd(3);
    tree->Draw("nhit_len:e_i>>hH(200,0,20, 200, 0, 1000)", "", "colz");    // time of b with errorbars

    cout << endl;
    cout << "You can now examine the structure of your tree in the TreeViewer" << endl;
    cout << endl;


    
    tree->Write();
    f->Write();

}
