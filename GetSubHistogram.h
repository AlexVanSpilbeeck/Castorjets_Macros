void GetSubHistogram(TH1D* original, TH1D* & newhist, double lowedge, double highedge){

  TH1D* originalhist = (TH1D*)original->Clone("Original");
  int bin = 0;
  int lowbin;

  while( originalhist->GetBinLowEdge( bin ) < lowedge){ bin++; }
  lowbin = bin;

  newhist = new TH1D(	originalhist->GetName(), 
			originalhist->GetTitle(), 
			originalhist->GetNbinsX() - lowbin + 1,
			originalhist->GetBinLowEdge(lowbin),
			originalhist->GetBinLowEdge( original->GetNbinsX()+1  )
  );
  for(; bin <= originalhist->GetNbinsX()+1; bin++){
    newhist->SetBinContent( bin - lowbin + 1, originalhist->GetBinContent( bin ) );
    newhist->SetBinError( bin - lowbin + 1, originalhist->GetBinError( bin ) );
  }

}




void GetSubHistogram(TH2D* original, TH2D* & newhist, double lowedge, double highedge){

  TH2D* originalhist = (TH2D*)original->Clone("Original");
  int bin = 0;
  int lowbin;

  while( originalhist->GetXaxis()->GetBinLowEdge( bin ) < lowedge){ bin++; }
  lowbin = bin;

  newhist = new TH2D(	originalhist->GetName(), 
			originalhist->GetTitle(), 
			originalhist->GetNbinsX() - lowbin + 1,
			originalhist->GetXaxis()->GetBinLowEdge(lowbin),
			originalhist->GetXaxis()->GetBinUpEdge( original->GetNbinsX()  ),
			originalhist->GetNbinsY() - lowbin + 1,
			originalhist->GetYaxis()->GetBinLowEdge(lowbin),
			originalhist->GetYaxis()->GetBinUpEdge( original->GetNbinsY()  )
  );
  for(int binx = lowbin; binx <= originalhist->GetNbinsX()+1; binx++){
    for(int biny = lowbin; biny <= originalhist->GetNbinsY()+1; biny++){
      newhist->SetBinContent( binx - lowbin + 1, biny - lowbin + 1, originalhist->GetBinContent( binx, biny ) );
      newhist->SetBinError( binx - lowbin + 1, biny - lowbin + 1, originalhist->GetBinError( binx, biny ) );
    }
  }
}
