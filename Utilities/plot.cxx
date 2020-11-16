//Macro for plotting the data from Cryosim

void plot(){


//__________________Plot of mt_________________

TCanvas *c1 = new TCanvas("c1","c1",800,600);

TGraph *g1 = new TGraph("./data/helium/data_tube10_m.txt");

g1->SetLineColor(kBlue);
g1->SetMarkerStyle(22);
g1->SetTitle("graph1");
//g1->Draw("AL*");



TGraph *g2 = new TGraph("./data/helium/thermosiphon10cm.txt", "%lg %*lg %*lg %*lg %*lg %lg");

g2->SetLineColor(kRed);
g2->SetMarkerStyle(23);
g2->SetTitle("graph2");

TMultiGraph *mg1 = new TMultiGraph();



mg1->Add(g1);
mg1->Add(g2);
mg1->SetTitle("; q [W/m²]; mt [kg/s]");

mg1->Draw("ALP");

TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);


legend->AddEntry(g1,"Data");
legend->AddEntry(g2,"Simulation");
legend->Draw(); 

//_____________________Plot of z________________________

TCanvas *c2 = new TCanvas("c2","c2",800,600);

TGraph *g3 = new TGraph("./data/helium/data_tube10_z.txt");

g3->SetLineColor(kBlue);
g3->SetMarkerStyle(22);
g3->SetTitle("graph1");



TGraph *g4 = new TGraph("./data/helium/thermosiphon10cm.txt", "%lg %*lg %*lg %*lg %lg");

g4->SetLineColor(kRed);
g4->SetMarkerStyle(23);
g4->SetTitle("graph2");

TMultiGraph *mg2 = new TMultiGraph();



mg2->Add(g3);
mg2->Add(g4);
mg2->SetTitle("; q [W/m²]; z [cm]");

mg2->Draw("ALP");

TLegend *legend2 = new TLegend(0.1,0.7,0.48,0.9);


legend2->AddEntry(g3,"Data");
legend2->AddEntry(g4,"Simulation");
legend2->Draw(); 


//_______________________Plot of x____________________________

TCanvas *c3 = new TCanvas("c3","c3",800,600);

TGraph *g5 = new TGraph("./data/helium/data_tube10_x.txt");

g5->SetLineColor(kBlue);
g5->SetMarkerStyle(22);
g5->SetTitle("graph1");



TGraph *g6 = new TGraph("./data/helium/thermosiphon10cm.txt", "%lg %*lg %*lg %lg");

g6->SetLineColor(kRed);
g6->SetMarkerStyle(23);
g6->SetTitle("graph2");

TMultiGraph *mg3 = new TMultiGraph();



mg3->Add(g5);
mg3->Add(g6);
mg3->SetTitle("; q [W/m²]; x");

mg3->Draw("ALP");

TLegend *legend3 = new TLegend(0.1,0.7,0.48,0.9);


legend3->AddEntry(g5,"Data");
legend3->AddEntry(g6,"Simulation");
legend3->Draw(); 

}
