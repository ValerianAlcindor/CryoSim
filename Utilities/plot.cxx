//Macro for plotting the data from Cryosim

void plot(){


TGraph *g1 = new TGraph("./data/helium/data_tube14_z.txt");

g1->Draw("AL*");


}
