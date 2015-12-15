#include "Track_Scholz2000.h"
#include "Particle.h"

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
using std::cerr;
using std::left;
using std::endl;
using std::ios;

#include <cstdlib>
using std::exit;

#include <sstream>
using std::ostringstream;

using std::string;

using namespace Survival;

const double Track_Scholz2000::DELTA = 1.7;
const double Track_Scholz2000::GAMMA = 0.062;	// um/MeV^DELTA
const double Track_Scholz2000::R_MIN = 0.01;	// um

Track_Scholz2000::Track_Scholz2000(const Particle &particle,
                                   const double density,
                                   double t)
{
    particleType = particle.type;

    const double AMU2MEV = 931.494027;  //MeV/amu
    particleEnergy = particle.e_c / particle.restEnergy * AMU2MEV;
    r_max = GAMMA * pow(particle.e_c / particle.restEnergy * AMU2MEV, DELTA);  //um

    if( r_max < R_MIN )
        r_max = R_MIN;

    let = particle.let;
    e_c = particle.e_c;

    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    lambda = CONV / (2 * M_PI * density * (0.5 + log(r_max/R_MIN))); //Gy * um^3/MeV

    x_track = particle.x;  // mm
    y_track = particle.y;  // mm
    weight = particle.weight;
    
    time = t;
}

//------------------------------------------------------------------------------

Track_Scholz2000* Track_Scholz2000::clone() const
{
    return new Track_Scholz2000(*this);
}

//------------------------------------------------------------------------------

double Track_Scholz2000::getDistance(const double localDose) const  // Gy -> um
{
    if( localDose < lambda * let / (r_max*r_max) )
        return r_max;
    if( localDose > lambda * let / (R_MIN*R_MIN) )
        return -1.;
    
    return sqrt(lambda * let / localDose);
}

//------------------------------------------------------------------------------

double Track_Scholz2000::getLocalDose(const double distance) const  // distance from track center
{
    if( R_MIN < distance )
    {
        if( distance <= r_max )
            return (lambda * let / (distance*distance));
        else
            return 0;
    }
    else
        return (lambda * let / (R_MIN*R_MIN));
}

//------------------------------------------------------------------------------

void Track_Scholz2000::getPosition(double &returnX,
                                   double &returnY) const
{
    returnX = x_track;
    returnY = y_track;
}

//------------------------------------------------------------------------------

double Track_Scholz2000::getRadialIntegral(const double r_begin,
                                           const double r_end) const
{
    if( r_begin < 0. || r_end < r_begin )
    {
        cerr << "double Track_Scholz2000::getRadialIntegral(const double, const double) const -- Invalid r_min, r_max specified." << endl;
        exit(1);
    }
    
    if( r_end == r_begin )
        return 0.;
    
    double r_stop = (r_end < r_max ? r_end : r_max);
    double dose = 0.;
    if( r_begin < R_MIN )
    {
        if( r_stop <= R_MIN )
            dose += lambda * let * (r_stop*r_stop - r_begin*r_begin);
        else
        {
            dose += lambda * let * (R_MIN*R_MIN - r_begin*r_begin);
            dose += 2 * lambda * let * log(r_stop / R_MIN);
        }
    }
    else if( r_begin < r_max )
        dose += 2 * lambda * let * log(r_stop / r_begin);
    
    return dose / (r_end*r_end - r_begin*r_begin);
}

//------------------------------------------------------------------------------

string Track_Scholz2000::saveTrack() const
{
    ostringstream save_name;
    save_name << "Scholz2000_" << particleType << "_" << scientific << setprecision(4) << let << "_" << e_c << ".track";
    
    ofstream outDoseFile(save_name.str().c_str(), ios::out);
    outDoseFile << scientific << setprecision(14) << left;
    
    double logDistance_min = -4.0;
    double logDistance_max = 4.0;
    double stepLogDistance = (logDistance_max - logDistance_min) / 300;
    
    for( double logDistance = logDistance_min; logDistance < logDistance_max; logDistance += stepLogDistance )
        outDoseFile << pow(10, logDistance) << "     "
                    << Track_Scholz2000::getLocalDose(pow(10, logDistance)) << endl;
    
    outDoseFile.close();
    return save_name.str();
}

//------------------------------------------------------------------------------

void Track_Scholz2000::setPosition(const double x,
                                   const double y)
{
    x_track = x;  // mm
    y_track = y;  // mm
}
