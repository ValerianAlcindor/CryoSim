//Macro for plotting the data from Cryosim

void plot(){


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
mg1->SetTitle("; q [cm]; z [cm]");

mg1->Draw("ALP");

TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);


legend->AddEntry(g1,"Data");
legend->AddEntry(g2,"Simulation");
legend->Draw(); 

}
