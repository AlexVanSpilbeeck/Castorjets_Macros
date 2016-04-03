#include <TCanvas.h>
#include <TLatex.h>

#include <iostream>
void Finish_canvas(TCanvas* &can_){

   cout << "\tfinish canvas\t" << endl;
   can_->cd();

   int iPosX = 22;

   float lumiTextSize     = 0.6;
   float lumiTextOffset   = 0.2;
   float cmsTextSize      = 0.75;
   float cmsTextOffset    = 0.1;  

   float H = can_->GetWh();
   float W = can_->GetWw();
   float l = can_->GetLeftMargin();
   float t = can_->GetTopMargin();
   float r = can_->GetRightMargin();
   float b = can_->GetBottomMargin();

   float relPosX    = 0.045;
   float relPosY    = 0.035;
   float relExtraDY = 1.2;

   float posX_=0;
   if( iPosX%10<=1 )
     {
       posX_ =   l + relPosX*(1-l-r);
     }
   else if( iPosX%10==2 )
     {
       posX_ =  l + 0.5*(1-l-r);
     }
   else if( iPosX%10==3 )
     {
       posX_ =  1-r - relPosX*(1-l-r);
     }
   float posY_ = 1-t - relPosY*(1-t-b);

   posX_ -= 0.1;

//   posX_ =   l + 0.045*(1-l-r)*W/H;
//   posY_ = 1-t - 0.045*(1-t-b);

   float xl_0 = posX_;
   float yl_0 = posY_ - 0.15;
   float xl_1 = posX_ + 0.15*H/W;
   float yl_1 = posY_;



   TLatex *   tex = new TLatex(1-r,1-t+lumiTextOffset*t,"7 TeV (pp)");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(1-r,1-t+lumiTextOffset*t,"7 TeV (pp)");

   tex = new TLatex(1-r,1-t+lumiTextOffset*t,"#mathcal{L}=0.12 nb^{-1}");
   tex->SetNDC();
   tex->SetTextAlign(33);
   tex->SetTextFont(42);
   tex->SetTextSize( lumiTextSize*t );
   tex->SetLineWidth(2);
   tex->DrawLatex(1-r,1-t+lumiTextOffset*t,"L=0.12 nb^{-1}");

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - 0.015,
	"CMS");
   tex->SetNDC();
   tex->SetTextAlign(11);
   tex->SetTextFont(61);
   tex->SetTextSize(0.06);
   tex->SetLineWidth(2);
   tex->DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, "CMS");

   tex = new TLatex( can_->GetLeftMargin() ,
	1. - can_->GetTopMargin()/1.5,
	"Preliminary");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(52);
   tex->SetTextSize(0.0456);
   tex->SetLineWidth(2);
   tex->DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, "Preliminary");

   tex = new TLatex(
	can_->GetLeftMargin() + 0.06*4,
	0.936,
	"anti-k_{t} (R=0.5)");
   tex->SetNDC();
   tex->SetTextAlign(11);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(can_->GetLeftMargin() + 0.06*4, 1-t+lumiTextOffset*t, "anti-k_{t} (R=0.5)");

   tex = new TLatex(can_->GetLeftMargin() + 0.06*4,
	0.936,
	"-6.6<#eta<-5.2");
   tex->SetNDC();
   tex->SetTextAlign(13);
   tex->SetTextFont(42);
   tex->SetTextSize(0.048);
   tex->SetLineWidth(2);
   tex->DrawLatex(can_->GetLeftMargin() + 0.06*4, 1-t+lumiTextOffset*t, "-6.6<#eta<-5.2");

   can_->Modified();
   can_->cd();
   can_->SetSelected(can_);
}

