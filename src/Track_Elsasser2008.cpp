#include "Track_Elsasser2008.h"
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

const double Track_Elsasser2008::CONV = 160.2177;	//(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
const double Track_Elsasser2008::DELTA = 1.7;		// exponent for calculation of maximum radius
const double Track_Elsasser2008::GAMMA = 0.062;		// um/MeV^DELTA :  factor for calculation of maximum radius
const double Track_Elsasser2008::R_C = 4e-2;		// um : maximum track core extension
const double Track_Elsasser2008::SIGMA = 4e-3;		// um : radical diffusion length

const int Track_Elsasser2008::MAX_LENGTH_DOSE_CURVE;
double Track_Elsasser2008::tmpLogRhoCurve[MAX_LENGTH_DOSE_CURVE];
double Track_Elsasser2008::tmpLogDoseCurve[MAX_LENGTH_DOSE_CURVE];

int Track_Elsasser2008::numberOfElsasser2008Tracks = 0;

Track_Elsasser2008::Track_Elsasser2008(const Particle &particle,
                                       const double rho,  // kg/dm^3 or g/cm^3
                                       const double dose_cutoff,  // Gy
                                       const int numberOfCurvePoints,
                                       const double integrationStepFactor,
                                       const double bessel_limit,
                                       const double sigmaNumber,
                                       double t)
{
    const double AMU2MEV = 931.494027;  //MeV/amu
    particleEnergy = particle.e_c / particle.restEnergy * AMU2MEV;
    particleType = particle.type;
    x_track = particle.x;  // mm
    y_track = particle.y;  // mm
    density = rho;

    double gamma = (particle.e_c + particle.restEnergy) / particle.restEnergy;
    double beta = sqrt( 1 - 1 / pow(gamma, 2) );
    r_min = R_C * beta;

    r_max = GAMMA * pow(particleEnergy, DELTA);  //um
    if( r_max < r_min )  // if r_max is too small, it is set equal to r_min
        r_max = r_min;

    let = particle.let;
    e_c = particle.e_c;
    lambda = CONV / (2 * M_PI * density * (0.5 + log(r_max/r_min))); //Gy * um^3/MeV

    doseCutoff = dose_cutoff;

    weight = particle.weight;

    if( numberOfCurvePoints > MAX_LENGTH_DOSE_CURVE )
    {
        cerr << endl
             << "Track_Elsasser2008: Error: MAX_LENGTH_DOSE_CURVE = " << MAX_LENGTH_DOSE_CURVE
             << ", required " << numberOfCurvePoints << "." << endl;
        exit(1);
    }

    lengthDoseCurve = numberOfCurvePoints;
    integrationStep = integrationStepFactor * r_min;
    besselLimit = bessel_limit;
    numberOfSigma = sigmaNumber;

    double logRhoMin = -4;
    double logRhoMax = 4;
    double stepCurve = (logRhoMax - logRhoMin) / (lengthDoseCurve - 1);

    int k = -1;
    do  {
        k++;
        tmpLogRhoCurve[k] = logRhoMin + k * stepCurve;
        tmpLogDoseCurve[k] = log10( normalizedPunctualDose( pow(10, tmpLogRhoCurve[k]) ) );
    } while( lambda * let * pow(10, tmpLogDoseCurve[k]) > doseCutoff );

    lengthDoseCurve = k+1;

    logRhoCurve = new double[lengthDoseCurve];
    logDoseCurve = new double[lengthDoseCurve];
    for( int i = 0; i < lengthDoseCurve; i++)
    {
        logRhoCurve[i] = tmpLogRhoCurve[i];
        logDoseCurve[i] = tmpLogDoseCurve[i];
    }

    if( lengthDoseCurve > 1 )
        r_eff = pow( 10, logRhoCurve[k-1] + (logRhoCurve[k]-logRhoCurve[k-1])*(log10(doseCutoff/(lambda*let))-logDoseCurve[k-1])/(logDoseCurve[k]-logDoseCurve[k-1]) );
    else
        r_eff = pow(10, logRhoCurve[0]);

    numberOfElsasser2008Tracks++;
    
    time = t;
}

//------------------------------------------------------------------------------

Track_Elsasser2008::Track_Elsasser2008(const Track_Elsasser2008 &track)  // copy constructor
{
    (*this) = track;

    double* logRhoCurve2 = new double[lengthDoseCurve];
    for( int t = 0; t < lengthDoseCurve; t++ )
        logRhoCurve2[t] = logRhoCurve[t];
    logRhoCurve = logRhoCurve2;

    double* logDoseCurve2 = new double[lengthDoseCurve];
    for( int t = 0; t < lengthDoseCurve; t++ )
        logDoseCurve2[t] = logDoseCurve[t];
    logDoseCurve = logDoseCurve2;
}

//------------------------------------------------------------------------------

Track_Elsasser2008::~Track_Elsasser2008()
{
    delete [] logRhoCurve;
    delete [] logDoseCurve;
    
    numberOfElsasser2008Tracks--;
}

//------------------------------------------------------------------------------

Track_Elsasser2008* Track_Elsasser2008::clone() const
{
    numberOfElsasser2008Tracks++;

    return new Track_Elsasser2008(*this);
}

//------------------------------------------------------------------------------

double Track_Elsasser2008::getDistance(const double localDose) const  // Gy -> um
{
    if( localDose > doseCutoff )
        return r_eff;
    
    double logLocalDose = log10(localDose) / (lambda * let);
    
    int k = 0;
    while( k < lengthDoseCurve )
    {
        if( logDoseCurve[k] < logLocalDose)
        {
            if( k == 0 )
                return -1.;
            else
                return pow(10, logRhoCurve[k-1] + (logRhoCurve[k] - logRhoCurve[k-1]) * (logLocalDose - logDoseCurve[k-1]) / (logDoseCurve[k] - logDoseCurve[k-1]) );
        }
        k++;
    }
    
    return -1.;
}

//------------------------------------------------------------------------------

double Track_Elsasser2008::getLocalDose(const double distance) const  // distance from track center um
{
    if( distance < pow(10, logRhoCurve[0]) )
        return lambda * let * pow(10, logDoseCurve[0]);

    if( distance > r_eff )
        return 0;

    double logDistance = log10(distance);

    int k = int( floor( (logDistance - logRhoCurve[0])/(logRhoCurve[1] - logRhoCurve[0]) ) );

    double logDose, logDose_1, logDose_2;

    if( k < lengthDoseCurve )
        logDose_1 = logDoseCurve[k];
    else
        return 0;

    if( k+1 < lengthDoseCurve )
        logDose_2 = logDoseCurve[k+1];
    else
        return 0;

    logDose = logDose_1 + (logDistance - logRhoCurve[k]) * (logDose_2 - logDose_1) / (logRhoCurve[k+1] - logRhoCurve[k]);

    return lambda * let * pow(10, logDose);
}

//------------------------------------------------------------------------------

double Track_Elsasser2008::getLocalDoseMeanTime()
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

double Track_Elsasser2008::getRadialIntegral(const double r_min,
                                             const double r_max) const
{
    if( r_min < 0. || r_max < r_min )
    {
        cerr << "double Track_Elsasser2008::getRadialIntegral(const double, const double) const -- Invalid r_min, r_max specified." << endl;
        exit(1);
    }

    if( r_max == r_min )
        return 0.;

    cerr << "Warning: implementation of double Track_Elsasser2008::getRadialIntegral(const double, const double) const not complete." << endl;

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

double Track_Elsasser2008::getRelativePrecision() const
{
    return (getTrackLet() - let)/let;
}

//------------------------------------------------------------------------------

string Track_Elsasser2008::saveTrack() const
{
    ostringstream save_name;
    save_name << "Elsasser2008_" << particleType << "_" << scientific << setprecision(4) << let << "_" << e_c << ".track";
    
    ofstream outDoseFile(save_name.str().c_str(), ios::out);
    outDoseFile << scientific << setprecision(14) << left;
    
    for( int t = 0; t < lengthDoseCurve; t++ )
        outDoseFile << pow(10, logRhoCurve[t]) << "     "
                    << lambda * let * pow(10, logDoseCurve[t]) << endl;
    
    outDoseFile.close();
    return save_name.str();
}

//------------------------------------------------------------------------------

void Track_Elsasser2008::setPosition(const double x,
                                     const double y)
{
    x_track = x;
    y_track = y;
}

//------------------------------------------------------------------------------

double Track_Elsasser2008::getTrackLet() const
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

double Track_Elsasser2008::normalizedDoseIntegralArgument(const double r,
                                                          const double r1) const
{
    double ak;
    double arg_bessel;
    double arg_integral;
    double d;

    if( r1 > r_max )
        return 0;

    if (r1 <= r_min)
        d = 1 / (r_min*r_min); // old fashioned dose profile
    else
        d = 1 / (r1*r1);
   
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

double Track_Elsasser2008::normalizedPunctualDose(const double distance) const
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

    do {
        r1 += dr1;
        argument_2 = normalizedDoseIntegralArgument(distance, r1);
        normalizedDose += (argument_1 + argument_2) * dr1 / 2;  // trapezoidal rule
        argument_1 = argument_2;
    } while( r1 < integralStop );

    return normalizedDose;
}
