#include "Track_KieferChatterjee.h"
#include "Particle.h"

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::sqrt;
using std::exp;

#include <iomanip>
using std::scientific;
using std::setprecision;

#include <fstream>
using std::ofstream;

#include <iostream>
using std::cerr;
using std::left;
using std::endl;
using std::ios;

#include <sstream>
using std::ostringstream;

#include <cstdlib>
using std::exit;

using std::string;

using namespace Survival;

const double Track_KieferChatterjee::DELTA = 1.7;
const double Track_KieferChatterjee::GAMMA = 0.0616;	// um/MeV^DELTA
const double Track_KieferChatterjee::ETA = 0.0116;		// um

Track_KieferChatterjee::Track_KieferChatterjee(const Particle &particle,
                                               const double density,
                                               double t)
{
    particleType = particle.type;

    const double AMU2MEV = 931.494027;  //MeV/amu
    particleEnergy = particle.e_c / particle.restEnergy * AMU2MEV;

    beta = sqrt( 1 - 1/((particle.e_c/particle.restEnergy + 1) * (particle.e_c/particle.restEnergy + 1)) ); // v/c

    r_penumbra = GAMMA * pow(particle.e_c / particle.restEnergy * AMU2MEV, DELTA);  //um

    r_core = ETA * beta;

    z_eff = particle.charge * (1 - exp(-125 * beta / pow(particle.charge, 2.0/3.0)));

    k_p = 1.25 * 0.0001 * (z_eff / beta)*(z_eff / beta); // [Gy um^2] (da controllare se le unità sono giuste)

    if( r_penumbra < r_core )
        r_penumbra = r_core; // da verificare se puo' accadere...

    let = particle.let;
    e_c = particle.e_c;

    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy (ok!)

    dose_core = ( let / density ) * CONV;
    dose_core -= 2 * M_PI * k_p * log(r_penumbra/r_core);
    dose_core /= (M_PI * r_core * r_core);

    x_track = particle.x;  // mm
    y_track = particle.y;  // mm
    weight = particle.weight;
    
    time = t;
}

//------------------------------------------------------------------------------

Track_KieferChatterjee* Track_KieferChatterjee::clone() const
{
    return new Track_KieferChatterjee(*this);
}

//------------------------------------------------------------------------------

double Track_KieferChatterjee::getDistance(const double localDose) const  // Gy -> um
{
    if( localDose < k_p / (r_penumbra*r_penumbra) )
        return r_penumbra;
    if( localDose > k_p / (r_core*r_core) )
        return -1.;
    
    return sqrt(k_p / localDose);
}

//------------------------------------------------------------------------------

double Track_KieferChatterjee::getLocalDose(const double distance) const  // distance from track center
{
    if( r_core < distance )
    {
        if( distance <= r_penumbra )
            return (k_p / (distance*distance));
        else
            return 0;
    }
    else
        return dose_core;
}

//------------------------------------------------------------------------------

void Track_KieferChatterjee::getPosition(double &returnX,
                                         double &returnY) const  // track position referred to the beam axis (mm)
{
    returnX = x_track;
    returnY = y_track;
}

//------------------------------------------------------------------------------

double Track_KieferChatterjee::getRadialIntegral(const double r_min,
                                                 const double r_max) const
{
    if( r_min < 0. || r_max < r_min )
    {
        cerr << "double Track_KieferChatterjee::getRadialIntegral(const double, const double) const -- Invalid r_min, r_max specified." << endl;
        exit(1);
    }
    
    if( r_max == r_min )
        return 0.;
    
    double r_end = (r_max < r_penumbra ? r_max : r_penumbra);
    double dose = 0.;
    if( r_min < r_core )
    {
        if( r_max <= r_core )
            dose += dose_core * (r_max*r_max - r_min*r_min);
        else
        {
            dose += dose_core * (r_core*r_core - r_min*r_min);
            dose += 2 * k_p * log(r_end / r_core);
        }
    }
    else if( r_min < r_penumbra )
        dose += 2 * k_p * log(r_end / r_min);
    
    return dose / (r_max*r_max - r_min*r_min);
}

//------------------------------------------------------------------------------

string Track_KieferChatterjee::saveTrack() const
{
    ostringstream save_name;
    save_name << "KieferChatterjee_" << particleType << "_" << scientific << setprecision(4) << let << "_" << e_c << ".track";
    
    ofstream outDoseFile(save_name.str().c_str(), ios::out);
    outDoseFile << scientific << setprecision(14) << left;
    
    double logDistance_min = log10(r_core/10);
    double logDistance_max = log10(r_penumbra*2);
    double stepLogDistance = (logDistance_max - logDistance_min) / 300;
    
    for( double logDistance = logDistance_min; logDistance < logDistance_max; logDistance += stepLogDistance )
        outDoseFile << pow(10, logDistance) << "     " << Track_KieferChatterjee::getLocalDose(pow(10, logDistance)) << endl;
    
    outDoseFile.close();
    return save_name.str();
}

//------------------------------------------------------------------------------

void Track_KieferChatterjee::setPosition(const double x,
                                         const double y)  // track position referred to the beam axis (mm)
{
    x_track = x;  // mm
    y_track = y;  // mm  
}
