/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Yasmine Demane     <ydemane@ikp.tu-darmstadt.de>
//*-- Date: 09/03/2020
//*-- Last Update: 30/07/2020
// --------------------------------------------------------------
// Comments: Plots
//
// --------------------------------------------------------------

#include "Rtypes.h"
#include "TGraph.h"
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdio.h>

void Plots() {

  // zsref vs q

  auto c1 = new TCanvas("c1", "c1", 200, 10, 800, 600);

  auto mg = new TMultiGraph();
  mg->SetTitle("Vaporization height vs Heat flux (Lhori=30cm, Lverti=70cm, "
               "D=6mm, 15cm above target)");

  auto plot2 = new TGraph("./data/thermosiphon10cm.txt", "%lg%*lg%*lg%*lg%lg%*lg");
  plot2->SetTitle("Heating height = 10cm");
  plot2->SetLineColor(kBlue);
  plot2->SetMarkerColor(kBlue);
  plot2->SetMarkerSize(0.8);
  plot2->Draw("APL*");

  /*   auto plot4 = new TGraph("thermosiphon20cm.txt","%lg%*lg%*lg%*lg%lg%*lg");
     plot4->SetTitle("Heating height = 20cm");
     plot4->SetLineColor(kGreen);
     plot4->SetMarkerColor(kGreen);
     plot4->SetMarkerSize(0.8);
     plot4->Draw("APL*");

     auto plot5 = new TGraph("thermosiphon30cm.txt","%lg%*lg%*lg%*lg%lg%*lg");
     plot5->SetTitle("Heating height = 30cm");
     plot5->SetLineColor(6);
     plot5->SetMarkerColor(6);
     plot5->SetMarkerSize(0.8);
     plot5->Draw("APL*");

     auto plot6 = new TGraph("thermosiphon40cm.txt","%lg%*lg%*lg%*lg%lg%*lg");
     plot6->SetTitle("Heating height = 40cm");
     plot6->SetLineColor(kRed);
     plot6->SetMarkerColor(kRed);
     plot6->SetMarkerSize(0.8);
     plot6->Draw("APL*");

     auto plot1 = new TGraph("thermosiphon50cm.txt","%lg%*lg%*lg%*lg%lg%*lg");
     plot1->SetTitle("Heating height = 50cm");
     plot1->SetLineColor(kYellow);
     plot1->SetMarkerColor(kYellow);
     plot1->SetMarkerSize(0.8);
     plot1->Draw("APL*");

     auto plot3 = new TGraph("thermosiphon60cm.txt","%lg%*lg%*lg%*lg%lg%*lg");
     plot3->SetTitle("Heating height = 60cm");
     plot3->SetLineColor(kBlack);
     plot3->SetMarkerColor(kBlack);
     plot3->SetMarkerSize(0.8);
     plot3->Draw("APL*"); */

  mg->Add(plot2, "APL*");
  /*  mg->Add(plot4,"APL*");
    mg->Add(plot5,"APL*");
    mg->Add(plot6,"APL*");
    mg->Add(plot1,"APL*");
    mg->Add(plot3,"APL*"); */
  mg->GetYaxis()->SetTitle("zsref (m)");
  mg->GetYaxis()->SetTitleOffset(1.3);
  plot2->GetYaxis()->CenterTitle(true);
  mg->GetXaxis()->SetTitle("q (W/m^2)");
  mg->GetXaxis()->SetTitleOffset(1.3);
  plot2->GetXaxis()->CenterTitle(true);
  mg->GetYaxis()->SetRangeUser(0.2, 0.25);
  mg->GetXaxis()->SetRange(0.0, 160.0);
  mg->Draw("APL*");

  c1->BuildLegend();
  gPad->SetGrid(1, 1);

  //------14----------------------------------------------------------------------------------------------------------------------------

  // x vs q (14mm)

  auto c2 = new TCanvas("c2", "c2", 200, 10, 800, 600);

  auto ms = new TMultiGraph();
  ms->SetTitle("Vapor quality vs Heat flux (Lhori=30cm, Lverti=70cm, D=6mm, "
               "15cm above target)");

  auto plot11 = new TGraph("./data/thermosiphon10cm.txt", "%lg%*lg%*lg%lg%*lg%*lg");
  plot11->SetTitle("Heating height = 10cm");
  plot11->SetLineColor(kBlue);
  plot11->SetMarkerColor(kBlue);
  plot11->SetMarkerSize(0.8);
  plot11->Draw("APL*");

  /*    auto plot22 = new
     TGraph("thermosiphon20cm.txt","%lg%*lg%*lg%lg%*lg%*lg");
      plot22->SetTitle("Heating height = 20cm");
      plot22->SetLineColor(kGreen);
      plot22->SetMarkerColor(kGreen);
      plot22->SetMarkerSize(0.8);
      plot22->Draw("APL*");

      auto plot33 = new TGraph("thermosiphon30cm.txt","%lg%*lg%*lg%lg%*lg%*lg");
      plot33->SetTitle("Heating height = 30cm");
      plot33->SetLineColor(6);
      plot33->SetMarkerColor(6);
      plot33->SetMarkerSize(0.8);
      plot33->Draw("APL*");

      auto plot44 = new TGraph("thermosiphon40cm.txt","%lg%*lg%*lg%lg%*lg%*lg");
      plot44->SetTitle("Heating height = 40cm");
      plot44->SetLineColor(kRed);
      plot44->SetMarkerColor(kRed);
      plot44->SetMarkerSize(0.8);
      plot44->Draw("APL*");

      auto plot55 = new TGraph("thermosiphon50cm.txt","%lg%*lg%*lg%lg%*lg%*lg");
      plot55->SetTitle("Heating height = 50cm");
      plot55->SetLineColor(kYellow);
      plot55->SetMarkerColor(kYellow);
      plot55->SetMarkerSize(0.8);
      plot55->Draw("APL*");

      auto plot66 = new TGraph("thermosiphon60cm.txt","%lg%*lg%*lg%lg%*lg%*lg");
      plot66->SetTitle("Heating height = 60cm");
      plot66->SetLineColor(kBlack);
      plot66->SetMarkerColor(kBlack);
      plot66->SetMarkerSize(0.8);
      plot66->Draw("APL*"); */

  ms->Add(plot11, "APL*");
  /*  ms->Add(plot22,"APL*");
    ms->Add(plot33,"APL*");
    ms->Add(plot44,"APL*");
    ms->Add(plot55,"APL*");
    ms->Add(plot66,"APL*");*/
  ms->GetYaxis()->SetTitle("x");
  ms->GetYaxis()->SetTitleOffset(1.3);
  plot11->GetYaxis()->CenterTitle(true);
  ms->GetXaxis()->SetTitle("q (W/m^2)");
  ms->GetXaxis()->SetTitleOffset(1.3);
  plot11->GetXaxis()->CenterTitle(true);
  ms->GetYaxis()->SetRangeUser(0.0, 0.);
  ms->GetXaxis()->SetLimits(0.0, 160.0);
  ms->Draw("APL*");

  c2->BuildLegend();
  gPad->SetGrid(1, 1);

  //-----14--------------------------------------------------------------------------------------------------------------------------

  // mt vs q (14mm)

  auto c3 = new TCanvas("c3", "c3", 200, 10, 800, 600);

  auto mb = new TMultiGraph();
  mb->SetTitle("Total mass flux vs Heat flux (Lhori=30cm, Lverti=70cm, D=6mm, "
               "15cm above target)");

  auto plot111 = new TGraph("./data/thermosiphon10cm.txt", "%lg%*lg%*lg%*lg%*lg%lg");
  plot111->SetTitle("Heating height = 10cm");
  plot111->SetLineColor(kBlue);
  plot111->SetMarkerColor(kBlue);
  plot111->SetMarkerSize(0.8);
  plot111->Draw("APL*");

  /* auto plot222 = new TGraph("thermosiphon20cm.txt","%lg%*lg%*lg%*lg%*lg%lg");
     plot222->SetTitle("Heating height = 20cm");
     plot222->SetLineColor(kGreen);
     plot222->SetMarkerColor(kGreen);
     plot222->SetMarkerSize(0.8);
     plot222->Draw("APL*");

     auto plot333 = new TGraph("thermosiphon30cm.txt","%lg%*lg%*lg%*lg%*lg%lg");
     plot333->SetTitle("Heating height = 30cm");
     plot333->SetLineColor(6);
     plot333->SetMarkerColor(6);
     plot333->SetMarkerSize(0.8);
     plot333->Draw("APL*");

     auto plot444 = new TGraph("thermosiphon40cm.txt","%lg%*lg%*lg%*lg%*lg%lg");
     plot444->SetTitle("Heating height = 40cm");
     plot444->SetLineColor(kRed);
     plot444->SetMarkerColor(kRed);
     plot444->SetMarkerSize(0.8);
     plot444->Draw("APL*");

     auto plot555 = new TGraph("thermosiphon50cm.txt","%lg%*lg%*lg%*lg%*lg%lg");
     plot555->SetTitle("Heating height = 50cm");
     plot555->SetLineColor(kYellow);
     plot555->SetMarkerColor(kYellow);
     plot555->SetMarkerSize(0.8);
     plot555->Draw("APL*");

     auto plot666 = new TGraph("thermosiphon60cm.txt","%lg%*lg%*lg%*lg%*lg%lg");
     plot666->SetTitle("Heating height = 60cm");
     plot666->SetLineColor(kBlack);
     plot666->SetMarkerColor(kBlack);
     plot666->SetMarkerSize(0.8);
     plot666->Draw("APL*");   */

  mb->Add(plot111, "APL*");
  /*  mb->Add(plot222,"APL*");
    mb->Add(plot333,"APL*");
    mb->Add(plot444,"APL*");
    mb->Add(plot555,"APL*");
    mb->Add(plot666,"APL*");*/
  mb->GetYaxis()->SetTitle("mt (kg/s)");
  mb->GetYaxis()->SetTitleOffset(1.3);
  plot111->GetYaxis()->CenterTitle(true);
  mb->GetXaxis()->SetTitle("q (W/m^2)");
  mb->GetXaxis()->SetTitleOffset(1.3);
  plot111->GetXaxis()->CenterTitle(true);
  mb->GetYaxis()->SetRangeUser(0.0, 0.00055);
  mb->GetXaxis()->SetLimits(0.0, 160.0);
  mb->Draw("APL*");

  c3->BuildLegend();
  gPad->SetGrid(1, 1);
}
