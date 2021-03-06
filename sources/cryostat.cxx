#include "cryostat.h"
#include <fstream>
#include <iostream>
#include <sstream>

CryoSim::cryostat::cryostat(char* material) {

  if (std::string(material) == "H2") {
    // Properties of H2 at 1atm and 20.238K
    beta = 0.0164; // Thermal expansion coefficient [K^-1]
    rhoL = 70.79; // Liquid density[kg/m^3]
    rhoG = 1.34; // Gas density [kg/m^3]
    Lv   = 445400; // Latent heat [J/Kg]
    muL  = 13.92E-6; // LViscosity, [Pa*s] = [kg/m*s]
    muG  = 10.93E-6; // GViscosity, [Pa*s] = [kg/m*s]
    HC   = 10300; // specific Heat capacity [J/K*kg]
  } else if (std::string(material) == "He") {
    // Properties of LHe4 at 1atm
    beta = 0.23; // Thermal expansion coefficient at 4,2K [K^-1]
    rhoL = 125; // density, LHe [kg/m^3]
    rhoG = 17; // density, He4 Gas [kg/m^3]
    Lv   = 20750; // Latent heat [J/Kg] (p64)
    //  Tres = 4.2; //Reservoir temperature [K]
    muL = 3.60E-6; // LViscosity, [Pa*s] = [kg/m*s]
    muG = 1.04E-6; // GViscosity, [Pa*s] = [kg/m*s]
    HC  = 4480; // specific Heat capacity, LHe [J/K*kg]
  }
  // general properties and constants
  g = 9.80665; // gravity [m/s^2]
  //  l    = 0.55; // Length of the supply line [m] (careful : length between
  //  the
  // heating and the top of the return line)
}

CryoSim::cryostat::~cryostat() {}

double CryoSim::cryostat::FrictionCoefficient(double mt) {
  double val = 0.079 / pow(abs(mt) * D / (TubeArea * muL), 0.25);
  return val;
}

// Energy Conservation____
double CryoSim::cryostat::EC(double mt, double VQ, double q,
                             double zsrefInput) {
  // homogene model
  // void fraction from homogeneous model
  VF1 = (VQ * rhoL) / (VQ * rhoL + (1 - VQ) * rhoG);

  double val = mt * HC * (Tsat - Tsref);
  val += mt * Lv * VQ; // Lv(z) ???
  val += mt * g * (zch - zsref) - q * M_PI * D * (zch - zsref);
  val += ((8. * pow(mt, 3)) / (pow(M_PI, 2.) * pow(D, 4.) * pow(rhoL, 2.)))
         * ((pow(VQ, 3.) * pow(rhoL, 2.)) / (pow(VF1, 2.) * pow(rhoG, 2.))
            + (pow((1 - VQ), 3.)) / (pow((1 - VF1), 2.)) - 1);

  return val;
}

// Pressure drop____
double CryoSim::cryostat::PressureDrop(double mtot0, double q,
                                       double zsrefInput) {

  // Liquid phase friction coefficient
  double Cflo = FrictionCoefficient(mtot0);
  double gap  = 0.15;

  // Pressure variation for the heated liquid
  double deltaP1phase
      = 2 * Cflo * pow(mtot0, 2) * zsref / (D * rhoL * pow(TubeArea, 2))
        + rhoL * g * zsref
              * (1 - beta * q * M_PI * D / (2 * mtot0 * HC) * zsref);

  // Pressure variation for the heated 2phases fluid
  double deltaP2phase = 2 * Cflo * pow(mtot0, 2) / (D * rhoL * pow(TubeArea, 2))
                            * (zch - zsref) * phisquare1
                        + pow(mtot0, 2) / ((pow(TubeArea, 2)) * rhoL)
                              * (pow(FinalVQ, 2) * rhoL / (VF1 * rhoG)
                                 + pow(1 - FinalVQ, 2) / (1 - VF1) - 1)
                        + g * (zch - zsref) * (VF1 * rhoG + (1 - VF1) * rhoL);

  // Pressure variation in the no heated line
  double deltaPriser = 2 * Cflo * pow(mtot0, 2) * (zmax - zch) * phisquare1
                           / (D * rhoL * pow(TubeArea, 2))
                       + (VF1 * rhoG + (1 - VF1) * rhoL) * g * (zmax - zch);

  double deltaPS = pow(mtot0, 2)
                   / (2 * pow(TubeArea, 2) * (VF1 * rhoG + (1 - VF1) * rhoL));

  // Return line pressure drop [Pa]
  double delta13 = deltaP1phase + deltaP2phase + deltaPriser + deltaPS
                   - delta13input; // + abovetarget;
  double val = delta13;

  return val;
}

// Temperature(Reservoir-Pres)
// Temperature(saturation-Pz)
double CryoSim::cryostat::Temperature(double P, char* material) {

  double x[3];
  double val;

  // Calculation of Temperature for helium
  if (std::string(material) == "He") {

    val = 3.146631;
    val += 1.357655 * pow(((log(P) - 10.3) / 1.9), 1);
    val += 0.413923 * (pow(((log(P) - 10.3) / 1.9), 2));
    val += 0.091159 * (pow(((log(P) - 10.3) / 1.9), 3));
    val += 0.016349 * (pow(((log(P) - 10.3) / 1.9), 4));
    val += 0.001826 * (pow(((log(P) - 10.3) / 1.9), 5));
    val -= 0.004325 * (pow(((log(P) - 10.3) / 1.9), 6));
    val -= 0.004973 * (pow(((log(P) - 10.3) / 1.9), 7));
  }

  // Calculation of Temperature for hydrogen
  else if (std::string(material) == "H2") {

    double AA = 15.46688;
    double BB = -1.013378E2;
    double C  = 5.432005E-2;
    double D  = -1.105632E-4;

    double a = C / D;
    double b = (AA - log(P)) / D;
    double c = BB / D;

    gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);

    val = x[1];
  }

  return val;
}

void CryoSim::cryostat::compute(double q, double Pres, double HL,
                                char* material, double& x, double& zsrefInput,
                                double& mt) {

  std::cout << "_________________Calculate q = " << q << "____________________"
            << std::endl;

  // What can be modified
  std::ifstream inputfile;
  inputfile.open("./data/inputcryostat.txt");
  int    i_row = 0;
  double Cryostatvariables[5];

  std::string line;
  while (std::getline(inputfile, line)) {
    std::stringstream stream(line);
    std::string       a;
    double            b;
    if (stream >> a >> b) {
      Cryostatvariables[i_row] = b;
      i_row++;
    }
  }

  D              = Cryostatvariables[0];
  double lsupply = Cryostatvariables[1];
  zmax = Cryostatvariables[2]; // maxium height (0 is at target height)
  zch  = Cryostatvariables[3]; // from 0 to maxium heated height (0 is at target
                              // height)
  mt = Cryostatvariables[4];

  double mtot0 = 10;
  double L     = 0.3; // Horizontal tube lenght [m]
  double v     = 0;

  // Pressure drop equation
  double eM;

  //_______________________Main loop for zref, mt and x______________________
  while (abs(mt - mtot0) >= 0.00001) {

    mtot0 = mt;

    // zsref = Vaporization height
    TubeArea = M_PI * pow(D / 2., 2);

    // Liquid phase friction coefficient
    double Cflo = FrictionCoefficient(mt);

    delta13input = HL;

    double frictionforce
        = 2 * Cflo * pow(mt, 2) * lsupply / (D * rhoL * pow(TubeArea, 2));

    double Pe = Pres + delta13input - frictionforce;

    z         = 0;
    double dz = 0.00001;
    double Pz;
    double Tz;

    while (abs(Temperature(Pz, material) - Tz) >= 0.0001) {
      z += dz;
      if (z > zmax) {
        break;
      }

      // Temperature and pressure calculations at z
      Pz = Pe + 2 * Cflo * pow(mt, 2) * z / (D * rhoL * pow(TubeArea, 2))
           + rhoL * g * z * (1 - beta * (q * M_PI * D / (2 * mt * HC)) * z);

      // Temperature at z [K]
      Tz = Temperature(Pres, material)
           + (q * M_PI * D - mt * g) / (mt * HC) * z;
    }

    zsref = z;
    Tsref = Tz;
    //__________________________Loop for vapor quality___________

    x = 0.98; // Start value must be between 0 and 1. for 1 homgeneous model
              // does not work, because VF1 becomes 1 and the energy equation is
              // divided by 0

    double dx = 0.000001;
    double VQ;
    Pz = Pe + 2 * Cflo * pow(mt, 2) * zch / (D * rhoL * pow(TubeArea, 2))
         + rhoL * g * zch * (1 - beta * (q * M_PI * D / (2 * mt * HC)) * zch);

    Tsat = Temperature(Pz, material);

    while (abs(x - VQ) >= 0.00001) {

      VQ = x;

      // Energie conservation equation (thesis p111 IV-23)
      double eX = EC(mt, VQ, q, zsref) * dx
                  / (EC(mt, VQ + dx, q, zsref) - EC(mt, VQ, q, zsref));

      x = VQ - eX;
    }

    VF1        = abs(x) * rhoL / (abs(x) * rhoL + (1 - abs(x)) * rhoG);
    phisquare1 = (1 + abs(x) * (rhoL - rhoG) / rhoG)
                 * pow(1 + abs(x) * (muL - muG) / muG, -0.25);

    //_________________________Calculation of mt__________________________

    // Total mass flux
    double dm = 0.0001;
    FinalVQ   = x;

    // Pressure drop equation (thesis p111 IV-22)
    // 0.1 is the minimum step size required to converge on the solution via
    // the differential equation solving method of Newton
    eM = 0.1 * PressureDrop(mtot0, q, zsref) * dm
         / (PressureDrop(mtot0 + dm, q, zsref) - PressureDrop(mtot0, q, zsref));
    mt = mtot0 - eM;

    v += 1;
    if (v > 100) {
      std::cout << "Max iteration reached it= " << v << " No convergence"
                << std::endl;
      break;
    }
  }

  //________________End main loop____________________

  std::cout << "z: " << zsref << std::endl;
  std::cout << "x: " << FinalVQ << std::endl;
  std::cout << "mt: " << mt << std::endl;
  zsrefInput = zsref;
}
