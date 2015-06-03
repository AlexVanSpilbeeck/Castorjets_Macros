{
 gSystem->Load("libFWCoreFWLite.so");
 AutoLibraryLoader::enable();
 gSystem->Load("libDataFormatsFWLite.so");
 gROOT->SetStyle ("Plain");
 gStyle->SetPalette(1);
 gStyle->SetOptTitle(0);
 gStyle->SetLegendBorderSize(0);
 gStyle->SetStatStyle(0);
 gStyle->SetTitleStyle(0);

gStyle->SetCanvasColor(-1);
gStyle->SetPadColor(-1);
gStyle->SetFrameFillColor(-1);
gStyle->SetHistFillColor(-1);
gStyle->SetTitleFillColor(-1);
gStyle->SetFillColor(-1);
gStyle->SetStatStyle(0);
gStyle->SetTitleStyle(0);



 gSystem->Load("libRooFit") ;
 using namespace RooFit ;
 gSystem->Load("~/RooUnfold-1.1.1/libRooUnfold");
 cout << "RooUnfold loaded" << endl;
}

