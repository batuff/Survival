#include "Nucleus_MKM.h"
#include "Nucleus_Pixel.h"
#include "Nucleus_Integral.h"
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

Nucleus_MKM::Nucleus_MKM(const CellLine &cellLine_,
                         const double xPosition,
                         const double yPosition)
: cellLine( cellLine_ ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition ),
  r_nucleus( 0.0 )
{
    //clog << "Creating Nucleus of type MKM... " << endl;

    domainRadius = cellLine.getDomainRadius();  // um
    double nucleusRadius = cellLine.getNucleusRadius();  // um
    numberOfDomains = nucleusRadius*nucleusRadius / (domainRadius*domainRadius);
    cellLine.getParameters_LQ_noDt(alpha_d, beta_d);
    
    alpha_d /= numberOfDomains;
    beta_d /= numberOfDomains;

    createDomains();

    cleanNucleus();

    //clog << "done" << endl;
}

//------------------------------------------------------------------------------

Nucleus_MKM::Nucleus_MKM(const CellLine &cellLine_,
                         const double domainRadius_,
                         const int numberOfDomains_,
                         const double xPosition,
                         const double yPosition)
: cellLine( cellLine_ ),
  domainRadius( domainRadius_ ),
  numberOfDomains( numberOfDomains_ ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition ),
  r_nucleus( 0.0 )
{
    //clog << "Creating Nucleus of type MKM... " << flush;

    cellLine.getParameters_LQ_noDt(alpha_d, beta_d);
    
    alpha_d /= numberOfDomains;
    beta_d /= numberOfDomains;

    createDomains();

    cleanNucleus();

    //clog << "done" << endl
    //     << endl << flush;
}

//------------------------------------------------------------------------------

Nucleus_MKM::~Nucleus_MKM()
{
    for( int i = 0; i < numberOfDomains; i++ )
        delete domains[i];
    
    delete [] domains;
    domains = NULL;
    
    delete domainCell;
    domainCell = NULL;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::createDomains()
{
    domains = new Nucleus_Integral*[numberOfDomains];
    
    domainCell = new CellLine(cellLine.getCellType(), domainRadius);
    domainCell->addParametrization_LQ_noDt(1., 1.);
    domainCell->setParametrization("LQ_noDt");
    
    vector<double> x(numberOfDomains, 0.);  // um
    vector<double> y(numberOfDomains, 0.);  // um
    domains[0] = new Nucleus_Integral(*domainCell, x[0]*1e-3, y[0]*1e-3);  // um -> mm
    double xTranslation = -2. * domainRadius;
    double yTranslation = 0.;
    bool isNewPosition;
    
    r_nucleus=0.0;
    
    for( int i = 1; i < numberOfDomains; i++ )
    {
        xTranslation = -xTranslation;
        yTranslation = -yTranslation;
        do
        {
            Nucleus_MKM::rotate(xTranslation, yTranslation);
            x[i] = x[i-1] + xTranslation;
            y[i] = y[i-1] + yTranslation;
            
            isNewPosition = true;
            for( int j = 0; j < i; j++ )
                if( abs(x[i] - x[j]) < 1e-12 && abs(y[i] - y[j]) < 1e-12 ) // ATTENZIONE ALLA TOLLERANZA NUMERICA!!!
                    isNewPosition = false;
        } while( !isNewPosition );
        
        domains[i] = new Nucleus_Integral(*domainCell, x[i]*1e-3, y[i]*1e-3);  // um -> mm
        
        //valuta l'r_nucleus effettivo, dal dominio più distante
        double r_nucleus_temp=sqrt(x[i]*x[i]+y[i]*y[i])+domainRadius;
        if(r_nucleus_temp>r_nucleus){r_nucleus=r_nucleus_temp;}
    }
}

//------------------------------------------------------------------------------

void Nucleus_MKM::rotate(double &xTranslation,
                         double &yTranslation)
{
    double x = xTranslation;
    double y = yTranslation;
    xTranslation = x * 0.5 + y * sqrt(3.)/2.;
    yTranslation = - x * sqrt(3.)/2. + y * 0.5;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::addBackgroundDose(const double dose)  // Gy
{
    for (int indexOfDomain = 0; indexOfDomain < numberOfDomains ; indexOfDomain++)
        domains[indexOfDomain]->addBackgroundDose(dose); // adds uniform background in domain of type integral.
}

//------------------------------------------------------------------------------

void Nucleus_MKM::addNucleusDoses(Nucleus_MKM &nucleus)
{
    if (nucleus.getRadius() != r_nucleus) {
        cout << "Warning: summing local dose of MKM nuclei with different radii." << endl;
        return;
    }
    
    if ((nucleus.getNumberOfDomains() != numberOfDomains) & (nucleus.getDomainRadius() != domainRadius)) {
        cout << "Warning: summing local dose of MKM nuclei with different number and size of the domains." << endl;
        return;
    }
    
    double dose = 0.0, lethals = 0.0, dummy;
    for (int indexOfDomain = 0; indexOfDomain < numberOfDomains ; indexOfDomain++) {
        nucleus.getDoseAndLethalForDomain(indexOfDomain, dose, dummy, lethals, dummy);
        domains[indexOfDomain]->addBackgroundDose(dose); // adds integral dose (uniform background) in domain of type integral.
    }
    
}

//------------------------------------------------------------------------------

void Nucleus_MKM::cleanNucleus()
{
    for( int i = 0; i < numberOfDomains; i++ )
        domains[i]->cleanNucleus();
    
    inNucleusCount = 0;
    intersectionCount = 0;
}

//------------------------------------------------------------------------------
//      nuova
Nucleus_MKM* Nucleus_MKM::clone(const CellLine& cellLine)
{
    Nucleus_MKM* nucleo=new Nucleus_MKM(cellLine); // VA CANCELLATO DA QUALCHE PARTE!!!
    return nucleo;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::distributeDose(const Track &track)
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

void Nucleus_MKM::distributeDose(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        Nucleus_MKM::distributeDose(tracks[i]);
}

//------------------------------------------------------------------------------

string Nucleus_MKM::getCellType() const
{
    return cellLine.getCellType();
}

//------------------------------------------------------------------------------

void Nucleus_MKM::getDoseAndLethalForDomain(int domainIndex,
                                            double &dose,
                                            double &doseUncertainty,
                                            double &lethal,
                                            double &lethalUncertainty) const
{
    double tmp1, tmp2, tmp3;
    
    if (domainIndex < numberOfDomains) {
        domains[domainIndex]->getDoseAndSurvival(dose, tmp1, tmp2, tmp3);
        lethal = alpha_d * dose + beta_d * dose*dose;
    } else {
        dose = -1;
        lethal = -1;
    }
    doseUncertainty = -1;
    lethalUncertainty = -1;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::getDosesAndLethals(vector<double> &doses,
                                     vector<double> &dosesUncertainty,
                                     vector<double> &lethals,
                                     vector<double> &lethalsUncertainty) const
{
    for (int i = 0; i < numberOfDomains; i++)
        Nucleus_MKM::getDoseAndLethalForDomain(i, doses[i], dosesUncertainty[i], lethals[i], lethalsUncertainty[i]);
}

//------------------------------------------------------------------------------

void Nucleus_MKM::getDoseAndSurvival(double &dose,
                                     double &doseUncertainty,
                                     double &survival,
                                     double &survivalUncertainty) const
{
    double tmp1, tmp2, tmp3;
    //vector<double> doseDomains(numberOfDomains, 0.);
    //vector<double> lethalDomains(numberOfDomains, 0.);
    double sumDose = 0;
    double sumLethal = 0;
    
    for( int i = 0; i < numberOfDomains; i++ )
    {
        domains[i]->getDoseAndSurvival(dose, tmp1, tmp2, tmp3);
        //doseDomains[i] = dose;
        sumDose += dose;
        //lethalDomains[i] = alpha_d * dose + beta_d * dose*dose;
        //sumLethal += lethalDomains[i];
        sumLethal += alpha_d * dose + beta_d * dose*dose;//aggiunto
    }
    
    dose = sumDose / numberOfDomains; // average dose.
    survival = exp(-sumLethal);
    doseUncertainty = -1;
    survivalUncertainty = -1;
}

//------------------------------------------------------------------------------

double Nucleus_MKM::getDoseForDomain(int indexOfDomain) const
{
    double tmp1, tmp2, tmp3, dose;
    
    if (indexOfDomain < numberOfDomains) {
        domains[indexOfDomain]->getDoseAndSurvival(dose, tmp1, tmp2, tmp3);
        return dose;
    } else
        return -1;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::getPosition(double &returnX,
                              double &returnY) const
{
    returnX = x_nucleus;
    returnY = y_nucleus;
}

//------------------------------------------------------------------------------

void Nucleus_MKM::saveLocalDose(const string fileName) const
{
    ofstream outLocalDoseFile(fileName.c_str(), ios::out);
    if(!outLocalDoseFile)
    {
        cerr << "File " << fileName << " could not be opened." << endl;
        exit(1);
    }
    
    outLocalDoseFile << setprecision(6) << scientific << left;
    
    double xDomainPosition = 0.0, yDomainPosition = 0.0; // mm
    double dose = 0.0, lethals = 0.0, dummy1, dummy2;
    for (int indexOfDomain = 0; indexOfDomain < numberOfDomains; indexOfDomain++) {
        domains[indexOfDomain]->getPosition(xDomainPosition, yDomainPosition);
        Nucleus_MKM::getDoseAndLethalForDomain(indexOfDomain, dose, dummy1, lethals, dummy2);
        outLocalDoseFile << indexOfDomain << " " << domains[indexOfDomain]->getRadius()
                         << " " << xDomainPosition * 1e3 << " " << yDomainPosition * 1e3 //um
                         << " " << dose << " " << lethals << endl;
    }
    
    outLocalDoseFile.close();
    
    // code to visualize the output file with gnuplot:
    // set size ratio -1; # in questo modo le unità hanno la stessa lunghezza sull'asse x e y (um)
    // set xlabel "x [um]";
    // set ylabel "y [um]";
    // set cblabel "Dose [Gy^-1]"; # color bar label
    // plot "<output_file>" u 3:4:5 with points pt 7 ps 3 lt palette notitle;
}
