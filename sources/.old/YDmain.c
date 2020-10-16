/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Yasmine Demane     <ydemane@ikp.tu-darmstadt.de>
//*-- Date: 09/03/2020
//*-- Last Update: 30/07/2020
// --------------------------------------------------------------
// Comments: Calculation of the total mass flux and the vapor
//           quality in a thermosiphon (for Helium)
//
// --------------------------------------------------------------

// Standard Lib
#include <cerrno>
#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "YDcommon.h"
#include "YDfunctions.h"

void YDmain(double q, double Pres, double HL, double& x, double& zsref,
            double& mt) {

  // What can be modify
  ///////////////////////////////////////////////

  mtot0 = 10;
  mt    = 0.0001;

  D    = 0.006;
  l    = 0.55;
  zmax = 0.55; //(l = zmax)

  L       = 0.3; // Horizontal lenght
  lsupply = 0.7;

  ///////////////////////////////////////////////

  // functionEC
  double var[4];
  double p[7];
  p[0] = HC;
  p[1] = Lv;
  p[2] = D;
  p[3] = rhoL;
  p[4] = rhoG;
  p[5] = g;

  // function PressureDrop
  double variable[5];
  double cst[11];
  cst[0] = rhoL;
  cst[1] = rhoG;
  cst[2] = zmax;
  cst[3] = D;
  cst[6] = g;
  cst[7] = beta;
  cst[8] = HC;
  cst[9] = muL;

  previousz = 10;
  lastz     = 15;

  double v = 0;
  eM       = 0.5;

  for (int y = 0; y < 1000; y++) {

    if (((abs(eM)) >= 0.0001) || ((abs(previousz - lastz)) >= 0.001)
        || ((x) < 0)) {

      mtot0 = mt;

      //__________Zsref = Vaporization
      // height______________________________________________________________________________________

      A = PI
          * pow((D / 2.), 2); // A2 = PI*pow((0.02/2.),2); Target cross section
      Cflo = 0.079
             / pow(((abs(mt)) * D) / (A * muL),
                   0.25); // Cflo2 = 0.079/pow(((abs(mt))*0.02)/(A2*muL),0.25);

      delta13input  = rhoL * g * lsupply;
      frictionforce = (2 * Cflo * pow(mt, 2) * (lsupply + L * 2 + 0.15))
                      / (D * rhoL * pow(A, 2));

      Pe = Pres + delta13input - frictionforce - rhoL * g * 0.15;

      // cout<<"rho*g*h = "<<rhoL*g*lsupply<<"   frictionforce =
      // "<<frictionforce<<endl;

      cst[4]  = A;
      cst[5]  = q;
      p[6]    = q;
      cst[10] = delta13input;

      double z      = 0;
      int    zsteps = 10000;
      double dz     = zmax / zsteps;

      previousz = zsref;

      for (int i = 0; i < zsteps; i++) {

        // Temperature and pressure calculations at z
        Pz = Pe + (2 * Cflo * pow(mt, 2) * z) / (D * rhoL * pow(A, 2))
             + (rhoL * g * z) * (1 - beta * ((q * PI * D) / (2 * mt * HC)) * z);
        Tsatu(Pz);
        Thorizontal = ((2 * PI * D * L * 2) / (mt * HC)); // q=4W/m^2
        Tz = Te(Pres) + ((q * PI * D - mt * g) / (mt * HC)) * (z) + Thorizontal;

        /*    double Cflo2 = 0.079/pow(((abs(0.00011))*D)/(A*muL),0.25);
            double Tz2 = Te(Pres) + ((10*PI*D-0.00011*g)/(0.00011*HC))*(0.35) +
           ((2*PI*D*L*2)/(0.00011*HC)); double Pz2 =
           Pe+(2*Cflo2*pow(0.00011,2)*0.35)/(D*rhoL*pow(A,2))+(rhoL*g*0.35)*(1-beta*((10*PI*D)/(2*0.00011*HC))*0.35);
            Tsatur(Pz2);
            cout<<"Tz2 = "<<Tz2<<" Tsat = "<<Tsatur(Pz2)<<endl; */

        // cout<<"Tz = "<<Tz<<" Tsatu = "<<Tsatu(Pz)<<" z = "<<z<<endl;

        if ((abs(Tsatu(Pz) - Tz)) <= 0.0001) {

          zsref  = z;
          var[3] = zsref;
          Tsref  = Tz;
          var[1] = Tsref;
          // cout<<"zsref = "<<zsref<<" Tsref1 = "<<Tsref<<endl;
          lastz = zsref;
          break;
        }

        if (Tz > (Tsatu(Pz))) {

          zsref  = z;
          var[3] = zsref;
          Tsref  = Tz;
          var[1] = Tsref;
          // cout<<"zsref = "<<zsref<<" Tsref2 = "<<Tsref<<endl;
          break;

        }

        else {

          z = z + dz;
        }
      }

      // cout<<" previousz = "<<previousz<<" lastz =  "<<lastz<<endl;

      //__________Vapor
      // quality_______________________________________________________________________________________________

      mV        = (q * PI * D * HL / Lv);
      double x1 = (mV / mt);

      // cout<<"x1 = "<<x1<<endl;

      // x0 and x have to be close to the expected value otherwise the algorithm
      // diverge
      x0        = 100; // 0.09
      x         = 0.1; // 0.01
      double dx = 0.1;

      z          = zsref;
      int hsteps = 10000;
      dz         = (zmax) / hsteps;

      for (int i = 0; i < hsteps; i++) {

        z = HL;

        var[2] = z;

        Pz = Pe + (2 * Cflo * pow(mt, 2) * z) / (D * rhoL * pow(A, 2))
             + (rhoL * g * z) * (1 - beta * ((q * PI * D) / (2 * mt * HC)) * z);
        var[0] = Tsatu(Pz);

        // Energie equation (thesis p111 IV-23)
        eX = (EC(mt, x0, p, var) * dx)
             / (EC(mt, x0 + dx, p, var) - EC(mt, x0, p, var));

        x = x0 - eX;

        // cout<<" x1 = "<<x<<endl;

        if ((abs(x - x0)) <= 0.00001) {

          // cout<<"zch = "<<zch<<" x = "<<x<<" x0 = "<<x0<<" eX = "<<eX<<"
          // Tsatu = "<<Tsatu(Pz)<<" Tsref = "<<Tsref<<endl;

          zch = HL;

          VF1 = ((abs(x)) * rhoL) / ((abs(x)) * rhoL + (1 - (abs(x))) * rhoG);
          phisquare1 = (1 + (abs(x)) * (rhoL - rhoG) / rhoG)
                       * pow((1 + (abs(x)) * (muL - muG) / muG), -0.25);
          // cout<<"x = "<<x<<" i = "<<i<<endl;

          break;

        }

        else {

          x0 = x;
        }
      }

      //__________Total mass
      // flux__________________________________________________________________________________________________

      double dm = 0.00001;

      variable[0] = zsref;
      variable[1] = zch;
      variable[2] = phisquare1;
      variable[3] = x;
      variable[4] = VF1;

      // Pressure drop equation (thesis p111 IV-22)
      eM = ((PressureDrop(mtot0, cst, variable)) * dm)
           / ((PressureDrop(mtot0 + dm, cst, variable))
              - (PressureDrop(mtot0, cst, variable)));

      mt = mtot0 - eM;

      // cout<<"mt = "<<mt<<"  mtot0 = "<<mtot0<<" eM = "<<eM<<endl;

    }

    else {

      cout << "Total mass flow  = " << mt << " [kg/s], "
           << "  "
           << " Vapor quality  = " << x << ",  "
           << "  Vaporization height = " << zsref << "[m] " << endl;
      break;
    }

    v = v + 1;

    if (v > 998) {

      mt    = 0;
      x     = 0;
      zsref = 0;
    }
  }
}
