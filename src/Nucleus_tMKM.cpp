#include "Nucleus_tMKM.h"
#include "Nucleus_Pixel.h"
#include "Nucleus_Integral_t.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Track.h"

#include <fstream>
using std::ofstream;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::abs;
using std::sqrt;

#include <iomanip>
using std::setprecision;

#include <iostream>
using std::cout;
using std::cerr;
using std::clog;
using std::endl;
using std::flush;
using std::scientific;
using std::left;
using std::ios;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;

using namespace Survival;

    
Nucleus_tMKM::Nucleus_tMKM(const CellLine &cellLineRef,
                           const double xPosition,
                           const double yPosition)
: cellLine( cellLineRef ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition ),
  r_nucleus( 0.0 )
{
    //clog << "Creating Nucleus of type tMKM... " << endl;

    double tmp;
    domainRadius = cellLine.getDomainRadius();  // um
    double nucleusRadius = cellLine.getNucleusRadius();  // um
    numberOfDomains = nucleusRadius*nucleusRadius / (domainRadius*domainRadius);
    cellLine.getParameters_LQ_noDt_T(alpha_d, beta_d, tmp);
    
    alpha_d /= numberOfDomains;
    beta_d /= numberOfDomains;
    
    createDomains();

    cleanNucleus();

    //clog << "done" << endl;
}

//------------------------------------------------------------------------------

Nucleus_tMKM::Nucleus_tMKM(const CellLine &cellLineRef,
                           const double domainRadius_,
                           const int numberOfDomains_,
                           const double xPosition,
                           const double yPosition)
: cellLine( cellLineRef ),
  domainRadius( domainRadius_ ),
  numberOfDomains( numberOfDomains_ ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition ),
  r_nucleus( 0.0 )
{
    //clog << "Creating Nucleus of type tMKM... " << flush;

    double tmp;
    
    cellLine.getParameters_LQ_noDt_T(alpha_d, beta_d, tmp);
    
    alpha_d /= numberOfDomains;
    beta_d /= numberOfDomains;
    
    createDomains();

    cleanNucleus();

    //clog << "done" << endl
    //     << endl << flush;
}

//------------------------------------------------------------------------------

Nucleus_tMKM::~Nucleus_tMKM()
{
    for( int i = 0; i < numberOfDomains; i++ )
        delete domains[i];

    delete [] domains;
    domains = NULL;

    delete domainCell;
    domainCell = NULL;
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::addBackgroundDose(const double dose,
                                     const double t)
{
    for (int indexOfDomain = 0; indexOfDomain < numberOfDomains ; indexOfDomain++)
        domains[indexOfDomain]->addBackgroundDose(dose, t);
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::cleanNucleus()
{
    for( int i = 0; i < numberOfDomains; i++ )
        domains[i]->cleanNucleus();
    
    inNucleusCount = 0;
    intersectionCount = 0;
}

//------------------------------------------------------------------------------

Nucleus_tMKM* Nucleus_tMKM::clone(const CellLine& cellLine) // NOTE: the "clone" is a NEW CLEAN nucleus with the same characteristics,
                                                            //       not really a clone.
{
    Nucleus_tMKM* nucleo = new Nucleus_tMKM(cellLine); //VA CANCELLATO DA QUALCHE PARTE!!!
    return nucleo;
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::createDomains()
{
    domains = new Nucleus_Integral_t*[numberOfDomains];
    
    domainCell = new CellLine(cellLine.getCellType(), domainRadius);
    domainCell->addParametrization_LQ_noDt_T(1., 1.);
    domainCell->setParametrization("LQ_noDt_T");
    
    vector<double> x(numberOfDomains, 0.);  // um
    vector<double> y(numberOfDomains, 0.);  // um
    domains[0] = new Nucleus_Integral_t(*domainCell, x[0]*1e-3, y[0]*1e-3);  // um -> mm
    double xTranslation = -2. * domainRadius;
    double yTranslation = 0.;
    bool isNewPosition;
    
    for( int i = 1; i < numberOfDomains; i++ )
    {
        xTranslation = -xTranslation;
        yTranslation = -yTranslation;
        do {
            Nucleus_tMKM::rotate(xTranslation, yTranslation);
            x[i] = x[i-1] + xTranslation;
            y[i] = y[i-1] + yTranslation;
            
            isNewPosition = true;
            for( int j = 0; j < i; j++ )
                if( abs(x[i] - x[j]) < 1e-12 && abs(y[i] - y[j]) < 1e-12 ) // ATTENZIONE ALLA TOLLERANZA NUMERICA!!!
                    isNewPosition = false;
        } while( !isNewPosition );
        
        domains[i] = new Nucleus_Integral_t(*domainCell, x[i]*1e-3, y[i]*1e-3);  // um -> mm
        
        //valuta l'r_nucleus effettivo, dal dominio più distante
        double r_nucleus_temp=sqrt(x[i]*x[i]+y[i]*y[i])+domainRadius;
        if(r_nucleus_temp>r_nucleus){r_nucleus=r_nucleus_temp;}
    }
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::distributeDose(const Track &track)
{
    double x_track;
    double y_track;
    
    track.getPosition(x_track, y_track);
    x_track = (x_track - x_nucleus) * 1e3;  // mm -> um
    y_track = (y_track - y_nucleus) * 1e3;  // mm -> um
    
    if( sqrt(pow(x_track,2) + pow(y_track,2)) <= r_nucleus )
        inNucleusCount++;
    else if( sqrt(pow(x_track,2) + pow(y_track,2)) <= r_nucleus + track.getRadius() )
        intersectionCount++;
    else
        return;
    
    for( int i = 0; i < numberOfDomains; i++ )
        domains[i]->distributeDose(track);
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::distributeDose(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        Nucleus_tMKM::distributeDose(tracks[i]);
}

//------------------------------------------------------------------------------
/*
double Nucleus_tMKM::evaluateG(double alpha_ion,
                               double beta_ion)
{
    double surv=0.0, dose=0.0, tmp2=0.0, tmp3=0.0;
    getDoseAndSurvival(dose, tmp2, surv, tmp3);
    return (-log(surv)-alpha_ion * dose)/(beta_ion * dose*dose);
}
*/
//------------------------------------------------------------------------------

string Nucleus_tMKM::getCellType() const
{
    return cellLine.getCellType();
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::getDoseAndSurvival(double &dose,
                                      double &doseUncertainty,
                                      double &survival,
                                      double &survivalUncertainty) const
{
    double sumDose = 0.0, sumLethal = 0.0, lethal=0.0;
    vector<double> doses, times;
    
    for( int dom = 0; dom < numberOfDomains; dom++ )
    {
        doses = domains[dom]->getDoses();
        times = domains[dom]->getTimes();
        if (times.size()!=0) {
            domains[dom]->getDose(dose);
            lethal = alpha_d*dose + beta_d*dose*dose;
        
            for (unsigned long int i=0; i<times.size()-1; ++i)
                for (unsigned long int j=i+1; j<times.size(); ++j)
                    lethal -= 2*beta_d * doses[i]*doses[j] * ( 1-exp(-cellLine.getAC()*(times[j]-times[i])) );
        }
        else {
            dose=0.0;
            lethal=0.0;
        }
        sumDose += dose;
        sumLethal += lethal;
    }
    dose = sumDose / numberOfDomains; // average dose.
    survival = exp(-sumLethal);
    doseUncertainty = -1;
    survivalUncertainty = -1;
}

//------------------------------------------------------------------------------

double Nucleus_tMKM::getDoseForDomain(int indexOfDomain) const
{
    double dose=0.0;
    
    if (indexOfDomain < numberOfDomains) {
        domains[indexOfDomain]->getDose(dose);
        return dose;
    } else
        return -1;
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::getPosition(double &returnX,
                               double &returnY) const
{
    returnX = x_nucleus;
    returnY = y_nucleus;
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::rotate(double &xTranslation,
                          double &yTranslation)
{
    double x = xTranslation;
    double y = yTranslation;
    xTranslation = x * 0.5 + y * sqrt(3.)/2.;
    yTranslation = - x * sqrt(3.)/2. + y * 0.5;
}

//------------------------------------------------------------------------------

void Nucleus_tMKM::saveLocalDose(const string fileName) const
{
    ofstream outLocalDoseFile(fileName.c_str(), ios::out);
    if(!outLocalDoseFile)
    {
        cerr << "File " << fileName << " could not be opened." << endl;
        exit(1);
    }
    
    outLocalDoseFile << setprecision(6) << scientific << left;
    
    /*
    double xDomainPosition = 0.0, yDomainPosition = 0.0; // mm
    double dose = 0.0, lethals = 0.0, dummy1, dummy2;
    for (int indexOfDomain = 0; indexOfDomain < numberOfDomains; indexOfDomain++) {
        domains[indexOfDomain]->getPosition(xDomainPosition, yDomainPosition);
        getDoseAndLethalForDomain(indexOfDomain, dose, dummy1, lethals, dummy2);
        outLocalDoseFile << indexOfDomain << " " << domains[indexOfDomain]->getRadius()
        << " " << xDomainPosition * 1e3 << " " << yDomainPosition * 1e3 //um
        << " " << dose << " " << lethals << endl;
    }
     */
    
    outLocalDoseFile.close();
    
    // code to visualize the output file with gnuplot:
    // set size ratio -1; # in questo modo le unità hanno la stessa lunghezza sull'asse x e y (um)
    // set xlabel "x [um]";
    // set ylabel "y [um]";
    // set cblabel "Dose [Gy^-1]"; # color bar label
    // plot "<output_file>" u 3:4:5 with points pt 7 ps 3 lt palette notitle;
}
