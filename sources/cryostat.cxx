#include "cryostat.h"

CryoSim::cryostat::cryostat() {
  // Properties of H2 at 1atm and 20.238K
  beta = 0.0164; // Thermal expansion coefficient [K^-1]
  rhoL = 70.79; // Liquid density[kg/m^3]
  rhoG = 1.34; // Gas density [kg/m^3]
  Lv   = 445400; // Latent heat [J/Kg]
  muL  = 13.92E-6; // LViscosity, [Pa*s] = [kg/m*s]
  muG  = 10.93E-6; // GViscosity, [Pa*s] = [kg/m*s]
  HC   = 10300; // specific Heat capacity [J/K*kg]
  g    = 9.80665; // gravity [m/s^2]
  D    = 0.006; // diameter of the return line [m]
  l    = 0.55; // Length of the supply line [m] (careful : length between the
            // heating and the top of the return line)
}
CryoSim::cryostat::~cryostat() {}

// Energy Conservation____
double CryoSim::cryostat::EC(double mt, double VQ, double q,
                             double zsrefInput) {

  // homogene model
  // void fraction from homogeneous model
  VF1 = (VQ * rhoL) / (VQ * rhoL + (1 - VQ) * rhoG);

  double val = mt * HC * (Tsat - Tsref);
  val += mt * Lv * VQ;
  val += mt * g * (z - zsref) - q * M_PI * D * (z - zsref);
  val += ((8. * pow(mt, 3)) / (pow(M_PI, 2.) * pow(D, 4.) * pow(rhoL, 2.)))
         * ((pow(VQ, 3.) * pow(rhoL, 2.)) / (pow(VF1, 2.) * pow(rhoG, 2.))
            + (pow((1 - VQ), 3.)) / (pow((1 - VF1), 2.)) - 1);

  return val;
}

// Pressure drop____
double CryoSim::cryostat::PressureDrop(double mtot0, double q,
                                       double zsrefInput) {
  zsref = zsrefInput;

  // Liquid phase friction coefficient
  double Cflo = 0.079 / pow(mtot0 * D / (TubeArea * muL), 0.25);
  double gap  = 0.15;

  // Pressure variation for the heated liquid
  double deltaP1phase
      = ((2 * Cflo * pow(mtot0, 2)) * (zsref)) / (D * rhoL * pow(TubeArea, 2))
        + ((rhoL * g * (zsref))
           * (1 - beta * ((q * M_PI * D) / (2 * mtot0 * HC)) * (zsref)));

  double abovetarget
      = ((2 * Cflo * pow(mtot0, 2)) * (gap)) / (D * rhoL * pow(TubeArea, 2))
        + ((rhoL * g * (gap))
           * (1 - beta * ((q * M_PI * D) / (2 * mtot0 * HC)) * (gap)));

  // Pressure variation for the heated 2phases fluid
  double deltaP2phase
      = ((2 * Cflo * pow(mtot0, 2)) / (D * rhoL * pow(TubeArea, 2)))
            * ((zch - zsref) * phisquare1)
        + (pow(mtot0, 2) / ((pow(TubeArea, 2)) * rhoL))
              * (((pow(FinalVQ, 2)) * rhoL) / ((VF1 * rhoG))
                 + ((1 - FinalVQ) * (1 - FinalVQ)) / (1 - VF1) - 1)
        + g * (zch - zsref) * (VF1 * rhoG + (1 - VF1) * rhoL);

  // Pressure variation in the no heated line
  double deltaPriser = ((2 * Cflo * pow(mtot0, 2)) * (zmax - zch) * phisquare1)
                           / (D * rhoL * pow(TubeArea, 2))
                       + ((VF1 * rhoG + (1 - VF1) * rhoL) * g * (zmax - zch));

  double deltaPS = (pow(mtot0, 2))
                   / (2 * pow(TubeArea, 2) * (VF1 * rhoG + (1 - VF1) * rhoL));
  // Return line pressure drop [Pa]
  double delta13 = deltaP1phase + deltaP2phase + deltaPriser + deltaPS
                   - delta13input + abovetarget;
  double val = delta13;

  // cout<<"val: "<<val<<endl;
  return val;
}

// Temperature(Reservoir-Pres)
// Temperature(saturation-Pz)
double CryoSim::cryostat::Temperature(double P) {

  double AA = 15.46688;
  double BB = -1.013378E2;
  double C  = 5.432005E-2;
  double D  = -1.105632E-4;

  double a = C / D;
  double b = (AA - log(P)) / D;
  double c = BB / D;
  double x[3];
  double xT[3];

  gsl_poly_solve_cubic(a, b, c, &x[0], &x[1], &x[2]);
  return x[1];
}

void CryoSim::cryostat::compute(double q, double Pres, double HL, double& x,
                                double& zsrefInput, double& mt) {

  // What can be modified
  zmax           = 0.55; // (l = zmax)
  double mtot0   = 10;
  mt             = 0.0001;
  double L       = 0.3; // Horizontal tube lenght [m]
  double lsupply = 0.7; // Supply line length [m]

  double previousz = 10;
  double lastz     = 15;

  // Pressure drop equation
  double eM = 0.5;

  for (int y = 0; y < 1000; y++) {
    if (((abs(eM)) >= 0.0001) || ((abs(previousz - lastz)) >= 0.001)
        || ((x) < 0)) {
      mtot0 = mt;

      // zsref = Vaporization height
      TubeArea = M_PI * pow((D / 2.), 2);

      // Liquid phase friction coefficient
      double Cflo = 0.079 / pow(((abs(mt)) * D) / (TubeArea * muL), 0.25);

      delta13input         = rhoL * g * lsupply;
      double frictionforce = (2 * Cflo * pow(mt, 2) * (lsupply + L * 2 + 0.15))
                             / (D * rhoL * pow(TubeArea, 2));

      double Pe = Pres + delta13input - frictionforce - rhoL * g * 0.15;

      z             = 0;
      int    zsteps = 10000;
      double dz     = zmax / zsteps;

      previousz = zsref;
      for (int i = 0; i < zsteps; i++) {
        // Temperature and pressure calculations at z
        double Pz
            = Pe + (2 * Cflo * pow(mt, 2) * z) / (D * rhoL * pow(TubeArea, 2))
              + (rhoL * g * z)
                    * (1 - beta * ((q * M_PI * D) / (2 * mt * HC)) * z);
        Temperature(Pz);

        // additional radiation heat from the detecteur [K]
        double Thorizontal = ((2 * M_PI * D * L * 2) / (mt * HC)); // q=4W/m^2

        // Temperature at z [K]
        double Tz = Temperature(Pres)
                    + ((q * M_PI * D - mt * g) / (mt * HC)) * (z) + Thorizontal;

        if ((abs(Temperature(Pz) - Tz)) <= 0.0001) {
          zsref      = z;
          zsrefInput = zsref;
          Tsref      = Tz;
          lastz      = zsref;
          break;
        }
        if (Tz > (Temperature(Pz))) {
          zsref      = z;
          zsrefInput = zsref;
          Tsref      = Tz;
          break;
        } else {
          z = z + dz;
        }
      }

      // Vapor quality
      // vapor mass flux [kg/s]
      double mV = (q * M_PI * D * HL / Lv);
      double x1 = (mV / mt);

      // VQ and x have to be close to the expected value otherwise the
      // algorithm diverge
      double VQ = 100; // 0.09
      x         = 0.1; // 0.01
      double dx = 0.1;

      z          = zsref;
      int hsteps = 10000;
      dz         = zmax / hsteps;

      for (int i = 0; i < hsteps; i++) {
        z = HL;

        double Pz = Pe;
        Pz += (2 * Cflo * pow(mt, 2) * z) / (D * rhoL * pow(TubeArea, 2));
        Pz += (rhoL * g * z)
              * (1 - beta * ((q * M_PI * D) / (2 * mt * HC)) * z);
        Tsat = Temperature(Pz);

        // Energie conservation equation (thesis p111 IV-23)
        double eX = (EC(mt, VQ, q, zsref) * dx)
                    / (EC(mt, VQ + dx, q, zsref) - EC(mt, VQ, q, zsref));
        x = VQ - eX;

        if ((abs(x - VQ)) <= 0.00001) {
          zch = HL;
          // void fraction from homogeneous model
          VF1 = ((abs(x)) * rhoL) / ((abs(x)) * rhoL + (1 - (abs(x))) * rhoG);
          phisquare1 = (1 + (abs(x)) * (rhoL - rhoG) / rhoG)
                       * pow((1 + (abs(x)) * (muL - muG) / muG), -0.25);
          break;
        } else {
          VQ = x;
        }
      }

      // Total mass flux
      double dm = 0.00001;
      FinalVQ   = x;

      // Pressure drop equation (thesis p111 IV-22)
      eM = ((PressureDrop(mtot0, q, zsref)) * dm)
           / ((PressureDrop(mtot0 + dm, q, zsref))
              - (PressureDrop(mtot0, q, zsref)));
      mt = mtot0 - eM;
    } else {
      std::cout << "Total mass flow  = " << mt << " [kg/s], "
                << "  "
                << " Vapor quality  = " << x << ",  "
                << "  Vaporization height = " << zsref << "[m] " << std::endl;
      break;
    }
  }
}
