/////////////////////////////////////////////////////////////////
//*-- AUTHOR : Yasmine Demane     <ydemane@ikp.tu-darmstadt.de>
//*-- Date: 09/03/2020
//*-- Last Update: 30/07/2020 
// --------------------------------------------------------------
// Comments: Calculation of the total mass flux and the vapor
//     quality in a thermosiphon (for Helium)
//
// --------------------------------------------------------------

//Inputs
double q; //heatflux [W/m^2]=[J/s*m^2] 
double mtot0; //total mass flow rate [kg/s] 
double xmax = 0.4; // Maximum vapor quality              
double Pres; //Reservoir pressure [Pa] 
    
//Thermosiphon properties
double D ;//diameter of the return line [m]
double l; //Lenght of the supply line [m] (careful : length between the heating and the top of the return line)
double A; // tube cross section area 
double A2; //target cross section area
double zmax; // l = zmax
double HL; //heated part of the return line [m]  //(0.1m, 0.2m or 0.3m)
double L; // Horizontal tube lenght [m]
    
//Properties of H2 at 1atm and 20.238K
double beta = 0.0164; // Thermal expansion coefficient [K^-1]
double rhoL = 70.79;  //Liquid density[kg/m^3] 
double rhoG = 1.34; //Gas density [kg/m^3] 
double Lv = 445400; // Latent heat [J/Kg]  
double muL = 13.92E-6;//LViscosity, [Pa*s] = [kg/m*s]
double muG = 10.93E-6;//GViscosity, [Pa*s] = [kg/m*s]
double HC = 10300; ; // specific Heat capacity [J/K*kg] 
    
//Others parameters
double eX; //Energie conservation equation
double eM;//Pressure drop equation
double g = 9.80665; //gravity [m/s^2]   
double PI = TMath::Pi();
double mt; //total mass flow rate [kg/s] 
double mV; //vapor mass flux [kg/s]
double Tz; // Temperature at z [K]
double Tsat; // Saturation temperature [K]
double Thorizontal; //additional radiation heat from the detecteur [K]
double zsref; // Vaporization height [m] 
double Tsref; // Saturation temperature at zsref [K]
double x;//final vapor quality 
double Cflo;// Liquid phase friction coefficient 
double Cflo2;// Liquid phase friction coefficient for the target (different D)
double VF1; //void fraction from homogeneous model
double phisquare1; //Two phase multiplier for homogen model 
double x0; // Vapor quality
double delta13; //Return line pressure drop [Pa]
double delta13input;//
double deltaP1phase; //Pressure variation for the heated liquid 
double deltaP2phase; //Pressure variation for the heated 2phases fluid
double deltaPriser; //Pressure variation in the no heated line
double deltaPS; //Pressure variation at the end of the return line
double Pz; //Pressure at z [Pa]
double Pe; // Pressure at the bottom of the return line [Pa]
double pow ( double base, double exp );
double log (double x);
double zch; //supposed heat return line [m]
double frictionforce; //[Pa]
double lsupply; //Supply line length [m]
double lastz; 
double previousz; 



















