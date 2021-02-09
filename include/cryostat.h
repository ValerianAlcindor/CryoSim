#ifndef cryostat_h
#define cryostat_h

#include <cerrno>
#include <cfenv>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <math.h>

#include <gsl/gsl_poly.h>

namespace CryoSim {

class cryostat {

public:
  cryostat(char* material);
  ~cryostat();

public:
  // Inputs
  double beta; // Thermal expansion coefficient [K^-1]
  double rhoL; // Liquid density[kg/m^3]
  double rhoG; // Gas density [kg/m^3]
  double Lv; // Latent heat [J/Kg]
  double muL; // LViscosity, [Pa*s] = [kg/m*s]
  double muG; // GViscosity, [Pa*s] = [kg/m*s]
  double HC; // specific Heat capacity [J/K*kg]
  double g; // gravity [m/s^2]
  double D; // diameter of the return line [m]
  //  double l; // Length of the supply line [m] (careful : length between the
  // heating and the top of the return line)

  // Thermosiphon properties
  double TubeArea; // tube cross section area
  double zmax; // l = zmax
  double HL; // heated part of the return line [m]  //(0.1m, 0.2m or 0.3m)

  // Others parameters
  double Tsat; // Saturation temperature [K]
  double Tsref; // Saturation temperature at zsref [K]
  double z;
  double VF1; // void fraction from homogeneous model
  double delta13input; //
  double zsref; // Vaporization height [m]
  double zch; // supposed heat return line [m]
  double FinalVQ; // final vapor quality
  double phisquare1; // Two phase multiplier for homogen model

public:
  double FrictionCoefficient(double mt);
  double EC(double mt, double x0, double q, double zsrefInput);
  double PressureDrop(double mtot0, double q, double zsrefInput);
  double Temperature(double P, char* material);
  void   compute(double q, double Pres, double HL, char* material, double& x,
                 double& zsrefInput, double& mt);
};
} // namespace CryoSim
#endif
