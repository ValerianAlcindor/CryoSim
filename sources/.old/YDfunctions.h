/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Yasmine Demane     <ydemane@ikp.tu-darmstadt.de>
//*-- Date: 09/03/2020
//*-- Last Update: 30/07/2020
// --------------------------------------------------------------
// Comments: Calculation of the total mass flux and the vapor
//           quality in a thermosiphon (for Helium)
//
// --------------------------------------------------------------

#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

// Energy Conservation____
double EC(double mt, double x0, double* p, double* var) {

  double PI   = TMath::Pi();
  double HC   = p[0];
  double Lv   = p[1];
  double D    = p[2];
  double rhoL = p[3];
  double rhoG = p[4];
  double g    = p[5];
  double q    = p[6];

  double Tsat  = var[0];
  double Tsref = var[1];
  double z     = var[2];
  double zsref = var[3];

  // homogene model
  double VF1 = (x0 * rhoL) / (x0 * rhoL + (1 - x0) * rhoG);

  double val
      = (mt * HC * (Tsat - Tsref) + mt * Lv * x0
         + ((8 * pow(mt, 3)) / (pow(PI, 2) * pow(D, 4) * pow(rhoL, 2)))
               * ((pow(x0, 3) * pow(rhoL, 2)) / (pow(VF1, 2) * pow(rhoG, 2))
                  + ((1 - x0) * (1 - x0) * (1 - x0)) / ((1 - VF1) * (1 - VF1))
                  - 1)
         + mt * g * (z - zsref) - q * PI * D * (z - zsref));

  return val;
}

// Pressure drop____
double PressureDrop(double mtot0, double* cst, double* variable) {

  double PI           = TMath::Pi();
  double rhoL         = cst[0];
  double rhoG         = cst[1];
  double zmax         = cst[2];
  double D            = cst[3];
  double A            = cst[4];
  double q            = cst[5];
  double g            = cst[6];
  double beta         = cst[7];
  double HC           = cst[8];
  double muL          = cst[9];
  double delta13input = cst[10];

  double zsref      = variable[0];
  double zch        = variable[1];
  double phisquare1 = variable[2];
  double x          = variable[3];
  double VF1        = variable[4];

  double Cflo = 0.079 / pow(mtot0 * D / (A * muL), 0.25);
  double gap  = 0.15;

  double deltaP1phase
      = ((2 * Cflo * pow(mtot0, 2)) * (zsref)) / (D * rhoL * pow(A, 2))
        + ((rhoL * g * (zsref))
           * (1 - beta * ((q * PI * D) / (2 * mtot0 * HC)) * (zsref)));

  double abovetarget
      = ((2 * Cflo * pow(mtot0, 2)) * (gap)) / (D * rhoL * pow(A, 2))
        + ((rhoL * g * (gap))
           * (1 - beta * ((q * PI * D) / (2 * mtot0 * HC)) * (gap)));

  double deltaP2phase = ((2 * Cflo * pow(mtot0, 2)) / (D * rhoL * pow(A, 2)))
                            * ((zch - zsref) * phisquare1)
                        + (pow(mtot0, 2) / ((pow(A, 2)) * rhoL))
                              * (((pow(x, 2)) * rhoL) / ((VF1 * rhoG))
                                 + ((1 - x) * (1 - x)) / (1 - VF1) - 1)
                        + g * (zch - zsref) * (VF1 * rhoG + (1 - VF1) * rhoL);

  double deltaPriser = ((2 * Cflo * pow(mtot0, 2)) * (zmax - zch) * phisquare1)
                           / (D * rhoL * pow(A, 2))
                       + ((VF1 * rhoG + (1 - VF1) * rhoL) * g * (zmax - zch));

  double deltaPS
      = (pow(mtot0, 2)) / (2 * pow(A, 2) * (VF1 * rhoG + (1 - VF1) * rhoL));

  double delta13 = deltaP1phase + deltaP2phase + deltaPriser + deltaPS
                   - delta13input + abovetarget;
  double val = delta13;

  // cout<<"val: "<<val<<endl;
  return val;
}

//___Solve cubic equation____
static double root3(double x) {
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
int SolveP3(double* x, double a, double b, double c) {

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

//____Temperature(Reservoir-Pres)____
double Te(double Pres) {

  double AA = 15.46688;
  double BB = -1.013378E2;
  double C  = 5.432005E-2;
  double D  = -1.105632E-4;

  double a = C / D;
  double b = (AA - log(Pres)) / D;
  double c = BB / D;
  double x[3];

  SolveP3(x, a, b, c);

  // cout<<"x[2]= "<<x[2]<<endl;
  return x[2];
}

//___Temperature(saturation-Pz)___
double Tsatu(double Pz) {

  double AA = 15.46688;
  double BB = -1.013378E2;
  double C  = 5.432005E-2;
  double D  = -1.105632E-4;

  double a = C / D;
  double b = (AA - log(Pz)) / D;
  double c = BB / D;
  double x[3];

  SolveP3(x, a, b, c);

  // cout<<"x[2]= "<<x[2]<<endl;
  return x[2];
}
