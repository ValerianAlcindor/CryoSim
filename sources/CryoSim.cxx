// Standard Lib
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "cryostat.h"

int main() {

  CryoSim::cryostat cryo;

  double q          = 0;
  double mt         = 0;
  double x          = 0;
  double zsrefInput = 0;
  double HL         = 0;

  const int size = 15;

  double tabq[size];
  double tabPres[size];
  double tabHL1[size];
  double tabHL2[size];
  double tabHL3[size];

  for (int i = 0; i < size; i++) {
    tabq[i]    = 10 + i * 10;
    tabPres[i] = 21462.9;
    tabHL1[i]  = 0.1;
    tabHL2[i]  = 0.2;
    tabHL3[i]  = 0.3;
  }

  FILE* fout = fopen("./data/thermosiphon10cm.txt", "w");

  fprintf(fout,
          "    %s                %s               %s               %s "
          " %s         %s\n",
          "q", "Pres", "HL", "x", "zsrefInput", "mt");

  for (int i = 0; i < 15; i++) {
    cryo.compute(tabq[i], tabPres[i], tabHL1[i], x, zsrefInput, mt);
    std::cout << "i :" << i << "  tabq[i] = " << tabq[i] << std::endl;

    if (x > 0) {
      // if (zsrefInput > 0) {
      if (mt > 0) {
        fprintf(fout, "%f       %f       %f         %f        %f          %f\n",
                tabq[i], tabPres[i], tabHL1[i], x, zsrefInput + 0.15, mt);
      }
      // }
    }
  }

  fclose(fout);

  return 0;
}
