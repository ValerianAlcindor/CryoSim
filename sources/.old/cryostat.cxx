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
double CryoSim::cryostat::EC(double mt, double VaporQuality, double q,
                             double zsrefInput) {

  // homogene model
  // void fraction from homogeneous model
  VF1 = (VaporQuality * rhoL)
        / (VaporQuality * rhoL + (1 - VaporQuality) * rhoG);

  double val = mt * HC * (Tsat - Tsref) + mt * Lv * VaporQuality
               + ((8 * pow(mt, 3)) / (pow(M_PI, 2) * pow(D, 4) * pow(rhoL, 2)))
                     * ((pow(VaporQuality, 3) * pow(rhoL, 2))
                            / (pow(VF1, 2) * pow(rhoG, 2))
                        + ((1 - VaporQuality) * (1 - VaporQuality)
                           * (1 - VaporQuality))
                              / ((1 - VF1) * (1 - VF1))
                        - 1)
               + mt * g * (z - zsref) - q * M_PI * D * (z - zsref);
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
              * (((pow(FinalVaporQuality, 2)) * rhoL) / ((VF1 * rhoG))
                 + ((1 - FinalVaporQuality) * (1 - FinalVaporQuality))
                       / (1 - VF1)
                 - 1)
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

//___Solve cubic equation____
double CryoSim::cryostat::root3(double x) {
  double s = 1.;
  while (x < 1.) {
    x *= 8.;
    s *= 0.5;
  }
  while (x > 8.) {
    x *= 0.125;
    s *= 2.;
  }
  double r = 1.5;
  for (int i = 0; i < 6; i++) {
    r -= 1. / 3. * (r - x / (r * r));
  }

  return r * s;
}

// solve cubic equation x^3 + a*x^2 + b*x + c = 0
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
int CryoSim::cryostat::SolveP3(double* x, double a, double b, double c) {

  double a2 = a * a;
  double q  = (a2 - 3 * b) / 9;
  double r  = (a * (2 * a2 - 9 * b) + 27 * c) / 54;

  // equation x^3 + q*x + r = 0
  double       r2 = r * r;
  double       q3 = q * q * q;
  double       A, B;
  const double eps = 1e-14;
  if (r2 <= (q3 + eps)) { //<<-- FIXED!
    double t = r / sqrt(q3);
    if (t < -1)
      t = -1;
    if (t > 1)
      t = 1;
    t = acos(t);
    a /= 3;
    q    = -2 * sqrt(q);
    x[0] = q * cos(t / 3) - a;
    x[1] = q * cos((t + 2 * M_PI) / 3) - a;
    x[2] = q * cos((t - 2 * M_PI) / 3) - a;
    return (3);
  } else {
    double b = fabs(r) + sqrt(r2 - q3);
    if (b > 0) {
      A = -root3(b);
    } else if (b > 0) {
      A = root3(-b);
    }

    if (r < 0)
      A = -A;
    B = (A == 0 ? 0 : B = q / A);

    a /= 3;
    x[0] = (A + B) - a;
    x[1] = -0.5 * (A + B) - a;
    x[2] = 0.5 * sqrt(3.) * (A - B);
    if (fabs(x[2]) < eps) {
      x[2] = x[1];
      return (2);
    }
    return (1);
  }
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
  // SolveP3(x, a, b, c);
  // std::cout << x[2] << " " << xT[0] << " " << xT[1] << " "<< xT[2] <<
  // std::endl;
  return x[1];
  // return x[2];
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
  double v         = 0;

  // Pressure drop equation
  double eM = 0.5;

  for (int y = 0; y < 1000; y++) {
    if (((abs(eM)) >= 0.0001) || ((abs(previousz - lastz)) >= 0.001)
        || ((x) < 0)) {
      mtot0 = mt;

      // zsref = Vaporization height
      TubeArea = M_PI * pow((D / 2.), 2);
      // TargetArea = M_PI*pow((0.02/2.),2); Target cross section

      // Liquid phase friction coefficient
      double Cflo = 0.079 / pow(((abs(mt)) * D) / (TubeArea * muL), 0.25);

      // Liquid phase friction coeff for the target (different D)
      // double Cflo2 =
      // 0.079/pow(((abs(mt))*0.02)/(TargetArea2*muL),0.25);

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

      // VaporQuality and x have to be close to the expected value otherwise the
      // algorithm diverge
      double VaporQuality = 100; // 0.09
      x                   = 0.1; // 0.01
      double dx           = 0.1;

      z          = zsref;
      int hsteps = 10000;
      dz         = (zmax) / hsteps;

      for (int i = 0; i < hsteps; i++) {
        z = HL;

        double Pz
            = Pe + (2 * Cflo * pow(mt, 2) * z) / (D * rhoL * pow(TubeArea, 2))
              + (rhoL * g * z)
                    * (1 - beta * ((q * M_PI * D) / (2 * mt * HC)) * z);
        Tsat = Temperature(Pz);

        // Energie conservation equation (thesis p111 IV-23)
        double eX = (EC(mt, VaporQuality, q, zsref) * dx)
                    / (EC(mt, VaporQuality + dx, q, zsref)
                       - EC(mt, VaporQuality, q, zsref));

        x = VaporQuality - eX;

        if ((abs(x - VaporQuality)) <= 0.00001) {
          zch = HL;

          // void fraction from homogeneous model
          VF1 = ((abs(x)) * rhoL) / ((abs(x)) * rhoL + (1 - (abs(x))) * rhoG);
          phisquare1 = (1 + (abs(x)) * (rhoL - rhoG) / rhoG)
                       * pow((1 + (abs(x)) * (muL - muG) / muG), -0.25);
          break;
        } else {
          VaporQuality = x;
        }
      }

      // Total mass flux
      double dm         = 0.00001;
      FinalVaporQuality = x;

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
    v = v + 1;
    if (v > 998) {
      // mt    = 0;
      // x     = 0;
      // zsref = 0;
    }
  }
  // std::cout << VaporQuality << std::endl;
}
