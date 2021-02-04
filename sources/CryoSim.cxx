// Standard Lib
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

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

  const int size = 15;
  // const int size = 30;

  double tabq[size];
  double tabPres[size];
  double tabHL1[size];
  double tabHL2[size];
  double tabHL3[size];

  /*
      tabHL1[0]=1730.99;
      tabHL1[1]=1650.95;
      tabHL1[2]=1600.98;
      tabHL1[3]=1550.92;
      tabHL1[4]=1500.26;
      tabHL1[5]=1480.64;
      tabHL1[6]=1450.61;
      tabHL1[7]=1430.22;
      tabHL1[8]=1415.08;
      tabHL1[9]=1400.14;
      tabHL1[10]=1395.53;
      tabHL1[11]=1385.55;
      tabHL1[12]=1370.44;
      tabHL1[13]=1360.54;
      tabHL1[14]=1350.48;
  */


  //_________Input for Simulation_________

  if (std::string(argv[1]) == "He") {
    std::ifstream inputfile;
    inputfile.open("./data/helium/input10.txt");
    int i_row = 0;

    std::string line;
    while (std::getline(inputfile, line))
    {
	std::stringstream stream(line);
	double a, b, c;
	if(stream >> a >> b >> c)
	{
	tabq[i_row] = a;
	tabPres[i_row] = b;
	tabHL1[i_row] = c;
	i_row++;
	std::cout<<a<<std::endl;
	std::cout<<b<<std::endl;
	std::cout<<c<<std::endl;
	}
    }

/*    int    n_data = 0;
    int    i_row  = 0;
    double data;

    char temp;
    inputfile >> temp;
    inputfile >> temp;

    while (inputfile >> data) {

      if (n_data % 2 == 0) {

        tabq[i_row] = data;

      } else {

        tabPres[i_row] = data;
        i_row++;
      }

      n_data++;
    }
*/
  } else if (std::string(argv[1]) == "H2") {
    for (int i = 0; i < size; i++) {
 //     tabq[i]    = 10 + i * 10;
 //     tabPres[i] = 102000.; // before 21462.9; ??

    std::ifstream inputfile;
    inputfile.open("./data/hydrogen/input10.txt");
    int i_row = 0;

    std::string line;
    while (std::getline(inputfile, line))
    {
	std::stringstream stream(line);
	double a, b, c;
	if(stream >> a >> b >> c)
	{
	tabq[i_row] = a;
	tabPres[i_row] = b;
	tabHL1[i_row] = c;
	i_row++;
	std::cout<<a<<std::endl;
	std::cout<<b<<std::endl;
	std::cout<<c<<std::endl;
	}
    }



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

  for (int i = 0; i < size; i++) {
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
