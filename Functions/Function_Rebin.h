using namespace std;


void Rebin_to(TH1D* &original, int maxBins){

  int nbins_default = original->GetNbinsX();
  int nbins_current = nbins_default;
  int rebinner = 1;
  while( ! (nbins_current < maxBins &&  nbins_default%rebinner == 0) ){
    rebinner++;

    // Is the number of bins divisible by the current rebinning factor?
    if( nbins_default%rebinner == 0){
      nbins_current = nbins_default/rebinner;

      cout << "\tnbins_default%rebinner\t" << endl;

      // Did we reduce the number of bins to the desired number?
      if(  nbins_current < maxBins ){
	cout << "\tRebinning!" << endl;
        original->Rebin( rebinner );
      }
    }  
  }
}





