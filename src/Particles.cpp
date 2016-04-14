#include "Particles.h"
#include "Particle.h"
#include "usefulFunctions.h"

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;

#include <iostream>
using std::cout;
using std::endl;
using std::flush;
using std::cerr;

#include <fstream>
using std::ifstream;

#include <cctype>
using std::isspace;

#include <algorithm>
using std::remove_if;

using std::string;

using namespace Survival;

Particles::Particles(const int numberOfParticles)
{
    particleVector.reserve(numberOfParticles);
    
    spectrum_file = "no_File";
}

//------------------------------------------------------------------------------

Particles::Particles(const string file_name)
{
    loadSpectrum(file_name);
}

//------------------------------------------------------------------------------

void Particles::operator<<(const Particle &particle)
{
    particleVector.push_back(particle);
}

//------------------------------------------------------------------------------

void Particles::operator<<(const Particles &particles)
{
    for( int i = 0; i < particles.size(); i++ )
        particleVector.push_back( particles[i] );
}

//------------------------------------------------------------------------------

Particle& Particles::operator[](const int index)
{
    return particleVector[index];
}

//------------------------------------------------------------------------------

const Particle& Particles::operator[](const int index) const
{
    return particleVector[index];
}

//------------------------------------------------------------------------------

double Particles::getDoseAveragedLet() const
{
    double numerator = 0.0;
    double denominator = 0.0;
    double let;
    
    for( size_t i = 0; i < particleVector.size(); i++ )
    {
        let = particleVector[i].let;
        numerator += particleVector[i].weight * let*let;
        denominator += particleVector[i].weight * let;
    }
    
    if( denominator != 0. )
        return numerator / denominator;
    else
        return 0.;
}

//------------------------------------------------------------------------------
/*
Particles* Particles::getIons()
{
    cout << "Extracting ions... " << flush;
    
    Particles *ions = new Particles(1);
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        if( particleVector[i].charge >= 0 && particleVector[i].A >= 1 )
            *ions << particleVector[i];
    
    cout << ions->size() << " ions extracted" << endl
         << endl << flush;
    
    return ions;
}
*/
//------------------------------------------------------------------------------

Particles Particles::getIons()
{
    cout << "Extracting ions... " << flush;
    
    Particles ions;
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        if( particleVector[i].charge >= 0 && particleVector[i].A >= 1 )
            ions << particleVector[i];
    
    cout << ions.size() << " ions extracted" << endl
         << endl << flush;
    
    return ions;
}

//------------------------------------------------------------------------------

Particles* Particles::getIons(const int charge)
{
    cout << "Extracting particles with charge = " << charge << "... " << flush;
    
    Particles* ions = new Particles(1);
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        if( particleVector[i].charge == charge )
            *ions << particleVector[i];
    
    cout << ions->size() << " particles with charge = " << charge << " extracted" << endl
         << endl << flush;
    
    return ions;
}

//------------------------------------------------------------------------------

Particles* Particles::getIons(const int charge, const int A)
{
    cout << "Extracting particles with charge = " << charge
         << " and mass number = " << A << "... " << flush;
    
    Particles *ions = new Particles(1);
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        if( particleVector[i].charge == charge && particleVector[i].A == A )
            *ions << particleVector[i];
    
    cout << ions->size() << " particles with charge = " << charge
         << " and mass number = " << A << " extracted" << endl
         << endl << flush;
    
    return ions;
}

//------------------------------------------------------------------------------
/*
void Particles::getMeanDepositedEnergyPerEvent(const int eventNumber,
                                               double &meanEd,
                                               double &uncertaintyMeanEd) const
{
    meanEd = 0;
    double meanEdSquared = 0;
    double edCurrentEvent = 0;
    int currentEvent = -1;
    for( size_t i = 0; i < particleVector.size(); i++ )
    {
        if( particleVector[i].eventID != currentEvent )
        {
            meanEd += edCurrentEvent / eventNumber;
            meanEdSquared += pow(edCurrentEvent, 2) / eventNumber;
            
            edCurrentEvent = 0;
            currentEvent = particleVector[i].eventID;
        }
        
        edCurrentEvent += particleVector[i].e_d;
    }
    meanEd += edCurrentEvent / eventNumber;
    meanEdSquared += pow(edCurrentEvent, 2) / eventNumber;
    
    uncertaintyMeanEd = sqrt( (meanEdSquared - pow(meanEd, 2)) / eventNumber );
}
*/
//------------------------------------------------------------------------------

double Particles::getMeanLet() const
{
    double totalLet = 0.0;
    double totalWeight = 0.0;
    
    for( size_t i = 0; i < particleVector.size(); i++ )
    {
        totalLet += particleVector[i].weight * particleVector[i].let;
        totalWeight += particleVector[i].weight;
    }
    
    if( totalWeight != 0. )
        return totalLet / totalWeight;
    else
        return 0.;
}

//------------------------------------------------------------------------------
/*
void Particles::getMeanLetPerEvent(const int eventNumber,
                                   double &meanLet,
                                   double &uncertaintyMeanLet) const
{
    meanLet = 0;
    double meanLetSquared = 0;
    double letCurrentEvent = 0;
    int currentEvent = -1;
    for( size_t i = 0; i < particleVector.size(); i++ )
    {
        if( particleVector[i].eventID != currentEvent )
        {
            meanLet += letCurrentEvent / eventNumber;
            meanLetSquared += pow(letCurrentEvent, 2) / eventNumber;
            
            letCurrentEvent = 0;
            currentEvent = particleVector[i].eventID;
        }
        
        letCurrentEvent += particleVector[i].let;
    }
    meanLet += letCurrentEvent / eventNumber;
    meanLetSquared += pow(letCurrentEvent, 2) / eventNumber;
    
    uncertaintyMeanLet = sqrt( (meanLetSquared - pow(meanLet, 2)) / eventNumber );
}
*/
//------------------------------------------------------------------------------

double Particles::getTotalLet() const
{
    double totalLet = 0.0;
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        totalLet += particleVector[i].let;
    
    return totalLet;
}

//------------------------------------------------------------------------------

double Particles::getTotalWeight() const
{
    double totalWeight = 0.0;
    for( size_t i = 0; i < particleVector.size(); i++ )
        totalWeight += particleVector[i].weight;

    return totalWeight;
}

//------------------------------------------------------------------------------

Particles* Particles::getWithCoordinatesBetween(const double x_min,  // mm
                                                const double x_max,  // mm
                                                const double y_min,  // mm
                                                const double y_max)  // mm
{
    cout << "Extracting particles with " << x_min << " mm < x < " << x_max << " mm"
         << " and " << y_min << " mm < y < " << y_max << " mm... " << flush;
    
    Particles *particlesVoxel = new Particles(1);
    
    for( size_t i = 0; i < particleVector.size(); i++ )
        if( x_min <= particleVector[i].x && particleVector[i].x < x_max &&
            y_min <= particleVector[i].y && particleVector[i].y < y_max )
            *particlesVoxel << particleVector[i];
    
    cout << particlesVoxel->size() << " particles extracted" << endl
         << endl << flush;
    
    return particlesVoxel;
}

//------------------------------------------------------------------------------

Particles* Particles::getWithDistanceBetween(const double distance_min,  // mm
                                             const double distance_max)  // mm
{
    Particles *particlesAnnulus = new Particles(1);

    double distance;
    for( size_t i = 0; i < particleVector.size(); i++ )
    {
        distance = sqrt( pow(particleVector[i].x, 2) +
                         pow(particleVector[i].y, 2) );

        if( distance_min <= distance && distance < distance_max )
            *particlesAnnulus << particleVector[i];
    }

    return particlesAnnulus;
}

//------------------------------------------------------------------------------

void Particles::loadSpectrum(const string file_name)
{
    cout << "Loading spectrum file..." << endl;
    
    const double AMU2MEV = 931.494027;
    
    spectrum_file = file_name;
    
    ifstream file(spectrum_file.c_str(), ifstream::in);
    if (!file) {
        cerr << "File " << spectrum_file << " could not be opened." << endl;
        exit(1);
    }
    
    // check if header is correct
    string header;
    getline(file,header);
    header.erase(remove_if(header.begin(), header.end(), [](char x){return isspace(x);}), header.end());
    if (header!="typechargeAenergyLETweight")
    {
        cerr << "Wrong header in input file " << spectrum_file << endl;
        exit(1);
    }
    
    string ion_type = "";
    int charge = 0;
    int A = 0;
    double energy = 0.0;
    double LET = 0.0;
    double w = 0.0;
    int nLine=1;
    
    while (file >> ion_type >> charge >> A >> energy >> LET >> w)
    {
        if (A<0) {
            cerr << "Wrong input detected in spectrum file " << spectrum_file << " at line " << nLine << ":" << endl
                 << "  Negative mass number not acceptable." << endl;
            exit(1);
        }
        if (w<=0) {
            cerr << "Wrong input detected in spectrum file " << spectrum_file << " at line " << nLine << ":" << endl
                 << "  Negative or NULL weight not acceptable." << endl;
            exit(1);
        }
        
        Particle particle;
        particle.type = ion_type;
        particle.charge = charge;
        particle.A = A;
        particle.restEnergy = A * AMU2MEV;  // MeV
        particle.let = LET;
        particle.e_c = energy;
        particle.x = 0.0; // mm
        particle.y = 0.0; // mm
        particle.z = 0.0; // mm
        particle.weight = w;
        particleVector.push_back(particle);
        nLine++;
    }
    file.close();
}

//------------------------------------------------------------------------------

void Particles::reconstructIonLETandEnergy()
{
    for (size_t i=0; i<particleVector.size(); ++i) {
        if (particleVector[i].e_c<0)
        {
            if (particleVector[i].let>=0) {
                particleVector[i].let /= 1000.;
                particleVector[i].e_c = betheBloch_inv_Srim(particleVector[i].type, particleVector[i].let)*particleVector[i].A;
            }
            else {
                cerr << "Invalid particle construction: negative both let and energy." << endl;
                exit(1);
            }
        }
        else if (particleVector[i].let<0)
        {
            if (particleVector[i].e_c>=0)
                particleVector[i].let = betheBloch_Srim(particleVector[i].type, particleVector[i].e_c/particleVector[i].A);
        }
    }
}

//------------------------------------------------------------------------------

int Particles::size() const
{
    return particleVector.size();
}
