#include "Track_Elsasser2007.h"
#include "Particle.h"

#include <algorithm>
using std::max;
using std::min;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::sqrt;

#include <iomanip>
using std::scientific;
using std::setprecision;

#include <fstream>
using std::ofstream;

#include <iostream>
using std::left;
using std::endl;
using std::cerr;
using std::ios;

#include <cstdlib>
using std::exit;

#include <sstream>
using std::ostringstream;

using std::string;

using namespace Survival;

// constants static variables and precalculated feature indices
const double Track_Elsasser2007::CONV = 160.2177;   //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
const double Track_Elsasser2007::DELTA = 1.7;       // exponent for calculation of maximum radius
const double Track_Elsasser2007::GAMMA = 0.062;     // um/MeV^DELTA : factor for calculation of maximum radius
const double Track_Elsasser2007::R_MIN = 3e-4;      // um
const double Track_Elsasser2007::SIGMA = 4e-3;      // um : radical diffusion length

int Track_Elsasser2007::lengthMasterCurve;
double Track_Elsasser2007::integrationStep;  // um
double Track_Elsasser2007::besselLimit;  // limit between the calculation of Bessel I0 function by means of series development and the calculation with the asintotical exponential approximation
double Track_Elsasser2007::numberOfSigma;
const int Track_Elsasser2007::MAX_LENGTH_MASTER_CURVE;
double Track_Elsasser2007::logRhoMasterCurve[MAX_LENGTH_MASTER_CURVE];
double Track_Elsasser2007::logDoseMasterCurve[MAX_LENGTH_MASTER_CURVE];

double Track_Elsasser2007::tmpLogDoseTail[MAX_LENGTH_MASTER_CURVE];

int Track_Elsasser2007::numberOfElsasser2007Tracks = 0;

Track_Elsasser2007::Track_Elsasser2007(const Particle &particle,
                                       const double rho,  // kg/dm^3 or g/cm^3
                                       const double dose_cutoff,  // Gy
                                       const int numberOfPointsMC,
                                       const double integrationStepFactor,
                                       const double bessel_limit,
                                       const double sigmaNumber,
                                       double tt)
{
    const double AMU2MEV = 931.494027;  //MeV/amu
    particleEnergy = particle.e_c / particle.restEnergy * AMU2MEV;
    particleType = particle.type;
    x_track = particle.x;  // mm
    y_track = particle.y;  // mm
    density = rho;

    r_max = GAMMA * pow(particleEnergy, DELTA);  //um
    if( r_max < R_MIN )  // if r_max is too small, it is set equal to R_MIN
        r_max = R_MIN;

    let = particle.let;
    e_c = particle.e_c;
    lambda = CONV / (2 * M_PI * density * (0.5 + log(r_max/R_MIN))); //Gy * um^3/MeV

    doseCutoff = dose_cutoff;

    weight = particle.weight;

    if( numberOfElsasser2007Tracks == 0 )
        createMasterCurve(numberOfPointsMC,
                          integrationStepFactor,
                          bessel_limit,
                          sigmaNumber);
    else if(numberOfPointsMC != lengthMasterCurve ||
            integrationStepFactor*R_MIN != integrationStep ||
            bessel_limit != besselLimit ||
            sigmaNumber != numberOfSigma)
    {
        cerr << "Track_Elsasser2007: Warning: Indicated track parameters" << endl
             << "are not equal to the corresponding ones of the master curve in use." << endl
             << "Therefore they will be ignored." << endl << endl; 
    }

    bool generateTail = true;

    if( log10(r_max-numberOfSigma*SIGMA) > logRhoMasterCurve[0] )
    {
        lengthMC = int( ceil( (log10(r_max-numberOfSigma*SIGMA) - logRhoMasterCurve[0])/(logRhoMasterCurve[1] - logRhoMasterCurve[0]) ) );

        for( int k = 0; k < lengthMC; k++ )
            if( lambda * let * pow(10, logDoseMasterCurve[k]) < doseCutoff )
            {
                lengthMC = k+1;
                r_eff = pow( 10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k]-logRhoMasterCurve[k-1])*(log10(doseCutoff/(lambda*let))-logDoseMasterCurve[k-1])/(logDoseMasterCurve[k]-logDoseMasterCurve[k-1]) );
                generateTail = false;
                break;
            }
    }
    else
        lengthMC = 0;

    if( generateTail )
    {
        int t = -1;
        int k = lengthMC-1;
        do
        {
            t++;
            k++;
            tmpLogDoseTail[t] = log10( normalizedPunctualDose( pow(10, logRhoMasterCurve[k]) ) );
        }
        while( lambda * let * pow(10, tmpLogDoseTail[t]) > doseCutoff );

        lengthTail = t+1;

        logDoseTail = new double[lengthTail];
        for( int i = 0; i < lengthTail; i++ )
            logDoseTail[i] = tmpLogDoseTail[i];

        if( lengthTail > 1 )
            r_eff = pow( 10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k]-logRhoMasterCurve[k-1])*(log10(doseCutoff/(lambda*let))-logDoseTail[t-1])/(logDoseTail[t]-logDoseTail[t-1]) );
        else if( lengthMC )
            r_eff = pow( 10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k]-logRhoMasterCurve[k-1])*(log10(doseCutoff/(lambda*let))-logDoseMasterCurve[k-1])/(logDoseTail[0]-logDoseMasterCurve[k-1]) );
        else
            r_eff = pow(10, logRhoMasterCurve[0]);
    }
    else
    {
        logDoseTail = 0;
        lengthTail = 0;
    }

    numberOfElsasser2007Tracks++;
    
    time = tt;
}

//------------------------------------------------------------------------------

Track_Elsasser2007::Track_Elsasser2007(const Track_Elsasser2007 &track)  // copy constructor
{
   (*this) = track;

   double* logDoseTail2 = new double[lengthTail];
   for( int t = 0; t < lengthTail; t++ )
      logDoseTail2[t] = logDoseTail[t];
   logDoseTail = logDoseTail2;
}

//------------------------------------------------------------------------------

Track_Elsasser2007::~Track_Elsasser2007()
{
    delete [] logDoseTail;
    
    numberOfElsasser2007Tracks--;
}

//------------------------------------------------------------------------------

Track_Elsasser2007* Track_Elsasser2007::clone() const
{
   return new Track_Elsasser2007(*this);
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getDistance(const double localDose) const  // Gy -> um
{
    if( localDose > doseCutoff )
        return r_eff;
    
    double logLocalDose = log10(localDose) / (lambda * let);
    
    int k = 0;
    while( k < lengthMC )
    {
        if( logDoseMasterCurve[k] < logLocalDose)
        {
            if( k == 0 )
                return -1.;
            else
                return pow(10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k] - logRhoMasterCurve[k-1]) * (logLocalDose - logDoseMasterCurve[k-1]) / (logDoseMasterCurve[k] - logDoseMasterCurve[k-1]) );
        }
        
        k++;
    }
    
    while( k - lengthMC < lengthTail )
    {
        if( logDoseTail[k - lengthMC] < logLocalDose)
        {
            if( k - lengthMC == 0 )
            {
                if( lengthMC == 0 )
                    return -1.;
                else
                    return pow(10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k] - logRhoMasterCurve[k-1]) * (logLocalDose - logDoseMasterCurve[k-1]) / (logDoseTail[0] - logDoseMasterCurve[k-1]) );
            }
            else
                return pow(10, logRhoMasterCurve[k-1] + (logRhoMasterCurve[k] - logRhoMasterCurve[k-1]) * (logLocalDose - logDoseTail[k-1-lengthMC]) / (logDoseTail[k-lengthMC] - logDoseTail[k-1-lengthMC]) );
        }
        k++;
    }
    
    return -1.;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getLocalDose(const double distance) const  // distance from track center um
{
   if( distance > r_eff )
      return 0;

   if( distance < pow(10, logRhoMasterCurve[0]) )
   {
      if( lengthMC )
         return lambda * let * pow(10, logDoseMasterCurve[0]);
      else
         return lambda * let * pow(10, logDoseTail[0]);
   }

   double logDistance = log10(distance);

   int k = int( floor( (logDistance - logRhoMasterCurve[0])/(logRhoMasterCurve[1] - logRhoMasterCurve[0]) ) );

   double logDose, logDose_1, logDose_2;

   if( k < lengthMC )
      logDose_1 = logDoseMasterCurve[k];
   else if( k < lengthMC + lengthTail )
      logDose_1 = logDoseTail[k-lengthMC];
   else
      return 0;

   if( k+1 < lengthMC )
      logDose_2 = logDoseMasterCurve[k+1];
   else if( k+1 < lengthMC + lengthTail )
      logDose_2 = logDoseTail[k+1-lengthMC];
   else
      return 0;

   logDose = logDose_1 + (logDistance - logRhoMasterCurve[k]) * (logDose_2 - logDose_1) / (logRhoMasterCurve[k+1] - logRhoMasterCurve[k]);

   return lambda * let * pow(10, logDose);
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getLocalDoseMeanTime()
{
    const long int numberOfCicles = 1000000;
    double dr = r_eff / numberOfCicles;
    
    double startTime = double(clock()) / CLOCKS_PER_SEC;
    
    for( int i = 0; i <= numberOfCicles; i++ )
        getLocalDose(i * dr);
    
    double stopTime = double(clock()) / CLOCKS_PER_SEC;
    
    return (stopTime - startTime) / numberOfCicles;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getRadialIntegral(const double r_min,
                                             const double r_max) const
{
   if( r_min < 0. || r_max < r_min )
   {
      cerr << "double Track_Elsasser2007::getRadialIntegral(const double, const double) const -- Invalid r_min, r_max specified." << endl;
      exit(1);
   }

   if( r_max == r_min )
      return 0.;

   cerr << "Warning: implementation of double Track_Elsasser2007::getRadialIntegral(const double, const double) const not complete." << endl;

//    double r_end = (r_max < r_penumbra ? r_max : r_penumbra);
//    double dose = 0.;
//    if( r_min < r_core )
//    {
//       if( r_max <= r_core )
//          dose += dose_core * (r_max*r_max - r_min*r_min);
//       else
//       {
//          dose += dose_core * (r_core*r_core - r_min*r_min);
// 
// 
//          dose += 2 * k_p * log(r_end / r_core);
//       }
//    }
//    else if( r_min < r_penumbra )
//    {
//       dose += 2 * k_p * log(r_end / r_min);
//    }
// 
//    return dose / (r_max*r_max - r_min*r_min);
   return -1.;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getRelativePrecision() const
{
    return (getTrackLet() - let)/let;
}

//------------------------------------------------------------------------------

string Track_Elsasser2007::saveTrack() const
{
    ostringstream save_name;
    save_name << "Elsasser2007_" << particleType << "_" << scientific << setprecision(4) << let << "_" << e_c << ".track";
    
    ofstream outDoseFile(save_name.str().c_str(), ios::out);
    outDoseFile << scientific << setprecision(14) << left;
    
    int t;
    for( t = 0; t < lengthMC; t++ )
        outDoseFile << pow(10, logRhoMasterCurve[t]) << "     "
                    << lambda * let * pow(10, logDoseMasterCurve[t]) << endl;
    for( t = 0; t < lengthTail; t++ )
        outDoseFile << pow(10, logRhoMasterCurve[lengthMC + t]) << "     "
                    << lambda * let * pow(10, logDoseTail[t]) << endl;
    
    outDoseFile.close();
    return save_name.str();
}

//------------------------------------------------------------------------------

void Track_Elsasser2007::setPosition(const double x,
                                     const double y)
{
    x_track = x;
    y_track = y;
}

//------------------------------------------------------------------------------

void Track_Elsasser2007::createMasterCurve(const int numberOfPointsMC,
                                           const double integrationStepFactor,
                                           const double bessel_limit,
                                           const double sigmaNumber)
{
    if( numberOfElsasser2007Tracks ) //i.e. if numberOfElsasser2007Tracks !=0
    {
        cerr << endl
             << "Track_Elsasser2007: Error: Master curve already existing. First delete the actual " << numberOfElsasser2007Tracks << " Elsasser2007 tracks and then create a new master curve." << endl;
        exit(1);
    }
    
    if( numberOfPointsMC > MAX_LENGTH_MASTER_CURVE )
    {
        cerr << endl
             << "Track_Elsasser2007: Error: MAX_LENGTH_MASTER_CURVE = " << MAX_LENGTH_MASTER_CURVE
             << ", required " << numberOfPointsMC << "." << endl;
        exit(1);
    }
    
    lengthMasterCurve = numberOfPointsMC;
    integrationStep = integrationStepFactor * R_MIN;
    besselLimit = bessel_limit;
    numberOfSigma = sigmaNumber;
    
    double r_max_save = r_max;
    r_max = 1e10;
    
    double logRhoMin = -4;
    double logRhoMax = 4;
    double stepMC = (logRhoMax - logRhoMin) / (lengthMasterCurve - 1);
    
    for( int i = 0; i < lengthMasterCurve; i++ )
    {
        logRhoMasterCurve[i] = logRhoMin + i * stepMC;
        logDoseMasterCurve[i] = log10( normalizedPunctualDose(pow(10.0, logRhoMasterCurve[i])) );
    }
    
    r_max = r_max_save;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::getTrackLet() const
{
    double trackLet = 0;
    double r1;
    double r2 = 0;
    double i1;
    double i2 = getLocalDose(r2) * r2;
    
    for( double i = 0; i < 12; i += 5e-4 )
    {
        r1 = r2;
        i1 = i2;
        r2 = pow(10.0, -8 + i);
        i2 = getLocalDose(r2) * r2;
        trackLet += (i1 + i2) * (r2 - r1);
    }
    
    return M_PI * density * trackLet / CONV;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::normalizedDoseIntegralArgument(const double r,
                                                          const double r1) const
{
   double ak;
   double arg_bessel;
   double arg_integral;
   double d;

   if( r1 > r_max )
      return 0;

   if (r1 <= R_MIN) {d = 1/ (R_MIN*R_MIN);} // old fashioned dose profile
   else {d = 1/ (r1*r1);}
   
   arg_bessel = (r*r1/(SIGMA*SIGMA)); // argument of Bessel's function
 
   if (arg_bessel < besselLimit)
   {
      ak = 1;  // first element of Bessel's series (k=0)
      arg_integral = 0;

      for( int k = 1; ak != 0; k++)
      {
         arg_integral += ak;
         ak *= arg_bessel*arg_bessel / (4 * k*k);
      }

      arg_integral = 1/(SIGMA*SIGMA) * r1 * d * arg_integral * exp(-(r*r + r1*r1)/(2*SIGMA*SIGMA));
   }
   else
      arg_integral =  1/(SIGMA*SIGMA) * 1/sqrt(2*M_PI*arg_bessel) * r1 * d * exp(-(r*r + r1*r1 - 2*r*r1)/(2*SIGMA*SIGMA));

   return arg_integral;
}

//------------------------------------------------------------------------------

double Track_Elsasser2007::normalizedPunctualDose(const double distance) const
{
   double integralStart = max( 0.0, distance - numberOfSigma * SIGMA );
   double integralStop = min( r_max, distance + numberOfSigma * SIGMA );
   int nr1 = int( ceil( (integralStop-integralStart)/integrationStep ) );
   double dr1;
   if( nr1 > 0 )
      dr1 = (integralStop-integralStart) / nr1;
   else
      return 0;
   double normalizedDose = 0;
   double r1 = integralStart;
   double argument_1 = normalizedDoseIntegralArgument(distance, integralStart);
   double argument_2;

   do
   {
      r1 += dr1;
      argument_2 = normalizedDoseIntegralArgument(distance, r1);
      normalizedDose += (argument_1 + argument_2) * dr1 / 2;  // trapezoidal rule
      argument_1 = argument_2;
   }
   while( r1 < integralStop );

   return normalizedDose;
}
