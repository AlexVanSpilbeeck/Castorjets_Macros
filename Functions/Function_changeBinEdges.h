void ChangeBinEdges( TH1D* &hResult, TH1D* hData, TH1D* hTemplate){



}




void ChangeBinEdges( TH2D* &hResult, TH2D* hData, TH2D* hTemplate){

  if( hData->GetNbinsX() != hTemplate->GetNbinsX() || hData->GetNbinsY() != hTemplate->GetNbinsY() ){
    std::cout << "Changing bin edges not possible." << std::endl;
    std::cout << "X:\t" << hData->GetNbinsX() << "\t" << hTemplate->GetNbinsX() << std::endl;
    std::cout << "Y:\t" << hData->GetNbinsY() << "\t" << hTemplate->GetNbinsY() << std::endl;
  }

  else{ 
    hResult = (TH2D*)hTemplate->Clone( hData->GetName() );
    hResult->GetXaxis()->SetTitle( hData->GetXaxis()->GetTitle() );
    hResult->GetYaxis()->SetTitle( hData->GetYaxis()->GetTitle() );

    for(int binx = 0; binx <= hData->GetNbinsX()+1; binx++){
      for(int biny = 0; biny <= hData->GetNbinsY()+1; biny++){

        double bincontent = hData->GetBinContent( binx, biny );
        hResult->SetBinContent( binx, biny, bincontent );

      } // Loop over y-bins
    } // Loop over x-bins
  } // else.
}



void ChangeBinEdges( TH2D* &hResult, TMatrixD& mData, TH2D* hTemplate){

  if( mData.GetNcols() != hTemplate->GetNbinsX() || mData.GetNrows() != hTemplate->GetNbinsY() ){
    std::cout << "Changing bin edges not possible." << std::endl;
    std::cout << "X:\t" << mData.GetNcols() << "\t" << hTemplate->GetNbinsX() << std::endl;
    std::cout << "Y:\t" << mData.GetNrows() << "\t" << hTemplate->GetNbinsY() << std::endl;
  }

  else{ 
    hResult = (TH2D*)hTemplate->Clone( mData.GetName() );
    //hResult->GetXaxis()->SetTitle( hData->GetXaxis()->GetTitle() );
    //hResult->GetYaxis()->SetTitle( hData->GetYaxis()->GetTitle() );

    for(int binx = 0; binx < mData.GetNcols(); binx++){
      for(int biny = 0; biny < mData.GetNrows(); biny++){

        double bincontent = (mData[binx])[biny];
        std::cout << "(x,y)\t" << binx << "\t" << biny << "\t" << bincontent << std::endl;

        hResult->SetBinContent( binx+1, biny+1, bincontent );

      } // Loop over y-bins
    } // Loop over x-bins
  } // else.
}
