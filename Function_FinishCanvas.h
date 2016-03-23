#include <TCanvas.h>
#include <TLatex.h>

void Finish_canvas(TCanvas* &can_){
   can_->cd();

  TLatex *   tex = new TLatex(0.96,0.936,"7 TeV (pp)");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(0.96,0.936,"#mathcal{L}=0.12 nb^{-1}");
   tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(0.96,0.936,"L=0.12 nb^{-1}");

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - 0.015,
	"CMS");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - can_->GetTopMargin()/1.5,
	"Preliminary");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(
	can_->GetLeftMargin() + 0.06*4,
	0.936,
	"anti-k_{t} (R=0.5)");
   tex->SetNDC();
   tex->SetTextAlign(11);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(can_->GetLeftMargin() + 0.06*4,
	0.936,
	"-6.6<#eta<-5.2");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(can_->GetLeftMargin() + 0.06*4,
	1. - can_->GetTopMargin()/1.5,
	"-6.6<#eta<-5.2");

   can_->Modified();
   can_->cd();
   can_->SetSelected(can_);
}

//== Altered eta range.

void Finish_canvas_narrow(TCanvas* &can_){
   can_->cd();

  TLatex *   tex = new TLatex(0.96,0.936,"7 TeV (pp)");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(0.96,0.936,"#mathcal{L}=0.12 nb^{-1}");
   tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(0.96,0.936,"L=0.12 nb^{-1}");

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - 0.015,
	"CMS");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - can_->GetTopMargin()/1.5,
	"Preliminary");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(
	can_->GetLeftMargin() + 0.06*4,
	0.936,
	"anti-k_{t} (R=0.5)");
   tex->SetNDC();
   tex->SetTextAlign(11);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->Draw();

   tex = new TLatex(can_->GetLeftMargin() + 0.06*4,
	0.936,
	"-6.1<#eta<-5.7");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(can_->GetLeftMargin() + 0.06*4,
	1. - can_->GetTopMargin()/1.5,
	"-6.6<#eta<-5.2");

   can_->Modified();
   can_->cd();
   can_->SetSelected(can_);
}
