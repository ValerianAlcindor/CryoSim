// Standard Lib
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "cryostat.h"

int main(int argc, char* argv[]) {

  if (std::string(argv[1]) != "He" && std::string(argv[1]) != "H2") {
    std::cout
        << "Error: Wrong argument for Material. Add \"H2\" for hydrogen or "
           "\"He\" for Helium as an argument, when starting the program."
        << std::endl;
    return 0;
  }

  CryoSim::cryostat cryo(argv[1]);

  double q          = 0;
  double mt         = 0;
  double x          = 0;
  double zsrefInput = 0;
  double HL         = 0;

  std::vector<double> tabq;
  std::vector<double> tabPres;
  std::vector<double> tabHL1;
  std::vector<double> tabHL2;
  std::vector<double> tabHL3;

  //_________Input for Simulation_________
  if (std::string(argv[1]) == "He") {
    std::ifstream inputfile;
    inputfile.open("./data/helium/input10.txt");
    double a, b, c;
    while (inputfile >> a >> b >> c) {
      tabq.push_back(a);
      tabPres.push_back(b);
      tabHL1.push_back(c);
    }
  } else if (std::string(argv[1]) == "H2") {
    std::ifstream inputfile;
    inputfile.open("./data/hydrogen/input10.txt");
    double a, b, c;
    while (inputfile >> a >> b >> c) {
      tabq.push_back(a);
      tabPres.push_back(b);
      tabHL1.push_back(c);
    }
  }

  //________Output File__________
  FILE* fout;
  if (std::string(argv[1]) == "He") {
    fout = fopen("./data/helium/thermosiphon10cm.txt", "w");
  } else if (std::string(argv[1]) == "H2") {
    fout = fopen("./data/hydrogen/thermosiphon10cm.txt", "w");
  }

  fprintf(fout,
          "    %s                %s               %s               %s "
          " %s         %s\n",
          "q", "Pres", "HL", "x", "zsrefInput", "mt");

  //________Starting the computation_______
  for (int i = 0; i < tabq.size(); i++) {
    cryo.compute(tabq[i], tabPres[i], tabHL1[i], argv[1], x, zsrefInput, mt);

    if (x > 0) {
      if (mt > 0) {
        fprintf(fout, "%f       %f       %f         %f        %f          %f\n",
                tabq[i], tabPres[i], tabHL1[i], x, zsrefInput, mt);
      }
    }
  }
  fclose(fout);

  return 0;
}
