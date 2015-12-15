#include "Tracks.h"
#include "Track_Scholz2000.h"
#include "Track_Elsasser2007.h"
#include "Track_Elsasser2008.h"
#include "Track_KieferChatterjee.h"
#include "Particles.h"
#include "Particle.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <iomanip>
using std::setprecision;
using std::showpoint;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::fixed;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;

using namespace Survival;

Tracks::Tracks(const Particles &particles,
               const string trackType,
			   const double massDensity ) // kg/dm^3 or g/cm^3
{
    string text("Converting Particles in Tracks... ");

    cout << showpoint << fixed << setprecision( floor(log10(particles.size()/100)) );
    cout << text << "0%\r" << text << flush;

    density = massDensity ; // kg/dm^3 or g/cm^3

    if (particles.size()!=0)
        trackVector.reserve(particles.size());
    else {
        cerr << endl << "Error: 0 particle to be converted." << endl;
        exit(1);
    }
    

    for( int i = 0; i < particles.size(); i++ )
    {
        if(particles[i].charge > 0 &&
           particles[i].restEnergy > 1 &&
           particles[i].e_c > 0 &&
           particles[i].let > 0)//excluded( photons, e-, mu-, pi-, n, neutrins, e+ );   accepted( mu+, pi+, ions+ )
        {
            if( trackType == "Scholz2000" )
                trackVector.push_back( new Track_Scholz2000(particles[i], density) );
            else if( trackType == "Elsasser2007" )
                trackVector.push_back( new Track_Elsasser2007(particles[i], density) );
            else if( trackType == "Elsasser2008" )
                trackVector.push_back( new Track_Elsasser2008(particles[i], density) );
            else if( trackType == "KieferChatterjee" )
                trackVector.push_back( new Track_KieferChatterjee(particles[i], density) );
            else
            {
                cerr << "The indicated track type doesn't exist." << endl;
                exit(1);
            }

            cout << double(i+1)/particles.size()*100
                 << "%\r" << text << flush;
        }
    }
    
    spectrum_file = particles.getSpectrumFile();

    cout << trackVector.size() << " tracks of type " << trackType << " created" << endl
         << endl << flush;
}

//------------------------------------------------------------------------------

Tracks::Tracks(const int numberOfTracks,
               const double massDensity)
{
    density = massDensity;
    trackVector.reserve(numberOfTracks);
    spectrum_file = "no_File";
}

//------------------------------------------------------------------------------

Tracks::~Tracks()
{
    eraseAll();
}

//------------------------------------------------------------------------------

void Tracks::operator<<(const Track &track)
{
    trackVector.push_back( track.clone() );
}

//------------------------------------------------------------------------------

void Tracks::operator<<(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        trackVector.push_back( tracks[i].clone() );
}

//------------------------------------------------------------------------------

const Track& Tracks::operator[](const int index) const
{
    return *trackVector[index];
}

//------------------------------------------------------------------------------

void Tracks::eraseAll()
{
    for( size_t i = 0; i < trackVector.size(); i++ )
        delete trackVector[i];

    trackVector.erase(trackVector.begin(), trackVector.end());
}

//------------------------------------------------------------------------------

double Tracks::getDoseAveragedLet() const
{
    double num = 0.0;
    double den = 0.0;
    double LET = 0.0;
    for (size_t i = 0; i < trackVector.size(); i++)
    {
        LET = trackVector[i]->getLet();
        num += trackVector[i]->getWeight()*LET*LET;
        den += trackVector[i]->getWeight()*LET;
    }
    return num/den;
}

//------------------------------------------------------------------------------

double Tracks::getMeanEnergy() const
{
    double totEnergy = 0.0;
    for (size_t i=0; i<trackVector.size(); ++i)
        totEnergy += trackVector[i]->getWeight() * trackVector[i]->getKineticEnergy();
    return totEnergy/getTotalWeight();
}

//------------------------------------------------------------------------------

double Tracks::getMeanLet() const
{
    double totalLet = 0.0;

    for( size_t i = 0; i < trackVector.size(); i++ )
        totalLet += trackVector[i]->getWeight() * trackVector[i]->getLet();

    return totalLet / getTotalWeight();
}

//------------------------------------------------------------------------------

double Tracks::getSigmaDoseAveragedLet() const
{
    double num = 0.0;
    double den = 0.0;
    double LET = 0.0;
    double LETd = getDoseAveragedLet();
    for (size_t i = 0; i < trackVector.size(); i++)
    {
        LET = trackVector[i]->getLet();
        num += (LET-LETd)*(LET-LETd)*trackVector[i]->getWeight()*LET;
        den += trackVector[i]->getWeight()*LET;
    }
    return sqrt(num/den);
}

//------------------------------------------------------------------------------

double Tracks::getSigmaMeanEnergy() const
{
    double num = 0.0;
    double energy = 0.0;
    double mEnergy = getMeanEnergy();
    for (size_t i=0; i<trackVector.size(); ++i)
    {
        energy = trackVector[i]->getKineticEnergy();
        num += (energy-mEnergy)*(energy-mEnergy)*trackVector[i]->getWeight();
    }
    return sqrt(num/getTotalWeight());
}

//------------------------------------------------------------------------------

double Tracks::getSigmaMeanLet() const
{
    double num = 0.0;
    double LET = 0.0;
    double mLET = getMeanLet();
    for (size_t i = 0; i < trackVector.size(); i++)
    {
        LET = trackVector[i]->getLet();
        num += (LET-mLET)*(LET-mLET)*trackVector[i]->getWeight();
    }
    return sqrt(num/getTotalWeight());
}

//------------------------------------------------------------------------------

double Tracks::getTotalWeight() const
{
    double totalWeight = 0.0;

    for( size_t i = 0; i < trackVector.size(); i++ )
        totalWeight += trackVector[i]->getWeight();

    return totalWeight;
}

//------------------------------------------------------------------------------

bool Tracks::isMonoenergetic() const
{
    if (trackVector.size()==1)
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------

int Tracks::size() const
{
    return trackVector.size();
}
