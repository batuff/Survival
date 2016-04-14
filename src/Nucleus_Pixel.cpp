#include "Nucleus_Pixel.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Track.h"

#include <fstream>
using std::ofstream;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::abs;

#include <iomanip>
using std::setprecision;
using std::setw;

#include <iostream>
using std::cout;
using std::cerr;
using std::clog;
using std::endl;
using std::scientific;
using std::left;
using std::ios;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;

using namespace Survival;

Nucleus_Pixel::Nucleus_Pixel(const CellLine &cellLineRef,
                             const double xPosition,
                             const double yPosition,
                             const double pixelSide1,
                             const int scale1,
                             const double radius1,
                             const int scale2,
                             const double radius2,
                             const int scale3,
                             const double radius3)
: cellLine( cellLineRef ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition ),
  pixelSide_1( pixelSide1 ),
  scale_1( scale1 ),
  radius_1( radius1 ),
  pixelSide_2( pixelSide_1 * scale_1 ),
  scale_2( scale2 ),
  radius_2( radius2 ),
  pixelSide_3( pixelSide_2 * scale_2 ),
  scale_3( scale3 ),
  radius_3( radius3 ),
  pixelSide_4( pixelSide_3 * scale_3 )
{
    //clog << "Creating Nucleus of type Pixel... " << endl;

    r_nucleus = cellLine.getNucleusRadius();

    createPixels();

    cleanNucleus();

    //clog << "done" << endl;
}

//------------------------------------------------------------------------------

Nucleus_Pixel::~Nucleus_Pixel()
{
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
                delete [] pixelVector[k].v[l].v[m].v;
            delete [] pixelVector[k].v[l].v;
        }
        delete [] pixelVector[k].v;
    }
    delete [] pixelVector;
    pixelVector = 0;
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::addBackgroundDose(const double dose)  // Gy
{
    for( int k = 0; k < numberOfBiggestPixels; k++ )
        pixelVector[k].dose += dose;
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::addNucleusDoses(Nucleus_Pixel &nucleus)
{
    if (nucleus.getRadius() != r_nucleus) {
        cout << "Warning: summing local dose of nuclei with different radii." << endl;
        return;
    }
    
    if ((nucleus.getNumberOfSmallestPixels() != numberOfSmallestPixels) || (nucleus.getNumberOfBiggestPixels() != numberOfBiggestPixels)) {
        cout << "Warning: summing local dose of nuclei with different number of pixels." << endl;
        return;
    }
    
    // occorre riempire per prima la grigliatura più piccola per entrambi i nuclei:
    double dummy1, dummy2, dummy3, dummy4;
    nucleus.getDoseAndSurvival(dummy1, dummy2, dummy3, dummy4);
    getDoseAndSurvival(dummy1, dummy2, dummy3, dummy4);
    
    vector<double> doses, dosesUncertainty, lethals, lethalsUncertainty;
    nucleus.getDosesAndLethals(doses, dosesUncertainty, lethals, lethalsUncertainty);
    
    int index = 0;
    
    for( int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for( int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            for( int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
            {
                for( int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                {
                    pixelVector[k].v[l].v[m].v[n].dose += doses[index];
                    index++;
                }
                pixelVector[k].v[l].v[m].dose = 0;
            }
            pixelVector[k].v[l].dose = 0;
        }
        pixelVector[k].dose = 0;
    }
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::cleanNucleus()
{
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        pixelVector[k].dose = 0;
        
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            pixelVector[k].v[l].dose = 0;
            
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
            {
                pixelVector[k].v[l].v[m].dose = 0;
                
                for(int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                    pixelVector[k].v[l].v[m].v[n].dose = 0;
            }
        }
    }
    inNucleusCount = 0;
    intersectionCount = 0;
}

//------------------------------------------------------------------------------

Nucleus_Pixel* Nucleus_Pixel::clone(const CellLine& cellLine)
{
    Nucleus_Pixel* nucleo=new Nucleus_Pixel(cellLine); //VA CANCELLATO DA QUALCHE PARTE!!!
    return nucleo;
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::distributeDose(const Track &track)
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
    
    Pixel *v_3;
    Pixel *v_2;
    Pixel *v_1;
    
    double distance;
    
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        if( !intersection(pixelVector[k].x, pixelVector[k].y, pixelSide_4, x_track, y_track, radius_3, distance) )
        {
            if( distance <= track.getRadius() )
                pixelVector[k].dose += track.getLocalDose(distance);
        }
        else
        {
            v_3 = pixelVector[k].v;
            for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++)
            {
                if( !intersection(v_3[l].x, v_3[l].y, pixelSide_3, x_track, y_track, radius_2, distance) )
                {
                    if( distance <= track.getRadius() )
                        v_3[l].dose += track.getLocalDose(distance);
                }
                else
                {
                    v_2 = v_3[l].v;
                    for(int m = 0; m < v_3[l].numberOfSubPixels; m++)
                    {
                        if( !intersection(v_2[m].x, v_2[m].y, pixelSide_2, x_track, y_track, radius_1, distance) )
                        {
                            if( distance <= track.getRadius() )
                                v_2[m].dose += track.getLocalDose(distance);
                        }
                        else
                        {
                            v_1 = v_2[m].v;
                            for(int n = 0; n < v_2[m].numberOfSubPixels; n++)
                            {
                                distance = sqrt( pow(v_1[n].x - x_track, 2) + pow(v_1[n].y - y_track, 2) );
                                if( distance <= track.getRadius() )
                                    v_1[n].dose += track.getLocalDose(distance);
                            }
                        }
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::distributeDose(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        Nucleus_Pixel::distributeDose(tracks[i]);
}

//------------------------------------------------------------------------------

string Nucleus_Pixel::getCellType() const
{
    return cellLine.getCellType();
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::getDosesAndLethals(vector<double> &doses,
                                       vector<double> &dosesUncertainty,
                                       vector<double> &lethals,
                                       vector<double> &lethalsUncertainty) const
{
    double nucleusFraction = 1 / static_cast<double>(numberOfSmallestPixels);
    doses.resize(numberOfSmallestPixels, 0.0);
    dosesUncertainty.resize(numberOfSmallestPixels, 0.0);
    lethals.resize(numberOfSmallestPixels, 0.0);
    lethalsUncertainty.resize(numberOfSmallestPixels, 0.0);
    
    int index = 0;
    
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
            {
                for(int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                {
                    pixelVector[k].v[l].v[m].v[n].dose += pixelVector[k].v[l].v[m].dose +
                                                          pixelVector[k].v[l].dose +
                                                          pixelVector[k].dose;
                    
                    doses[index] = pixelVector[k].v[l].v[m].v[n].dose;
                    lethals[index] = - cellLine.getLogSurvival_X( doses[index] ) * nucleusFraction;
                    dosesUncertainty[index] = -1;
                    lethalsUncertainty[index] = -1;
                    index++;
                }
                pixelVector[k].v[l].v[m].dose = 0;
            }
            pixelVector[k].v[l].dose = 0;
        }
        pixelVector[k].dose = 0;
    }
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::getDoseAndSurvival(double &dose,
                                       double &doseUncertainty,
                                       double &survival,
                                       double &survivalUncertainty) const
{
    dose = 0;
    double sum_lethal = 0;
    double nucleusFraction = 1 / static_cast<double>(numberOfSmallestPixels);
    
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
            {
                for(int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                {
                    pixelVector[k].v[l].v[m].v[n].dose += pixelVector[k].v[l].v[m].dose +
                                                          pixelVector[k].v[l].dose +
                                                          pixelVector[k].dose;
                    
                    dose += pixelVector[k].v[l].v[m].v[n].dose * nucleusFraction;
                    sum_lethal += - cellLine.getLogSurvival_X( pixelVector[k].v[l].v[m].v[n].dose ) * nucleusFraction;
                }
                pixelVector[k].v[l].v[m].dose = 0;
            }
            pixelVector[k].v[l].dose = 0;
        }
        pixelVector[k].dose = 0;
    }
    
    survival = exp(-sum_lethal);
    
    doseUncertainty = -1;
    survivalUncertainty = -1;
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::getPosition(double &returnX,
                                double &returnY) const
{
    returnX = x_nucleus;
    returnY = y_nucleus;
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::saveLocalDose(const string fileName) const
{
    ofstream outLocalDoseFile(fileName.c_str(), ios::out);
    if(!outLocalDoseFile)
    {
        cerr << "File " << fileName << " could not be opened." << endl;
        exit(1);
    }
    outLocalDoseFile << setprecision(6) << scientific << left;
    
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
        {
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
            {
                for(int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                {
                    pixelVector[k].v[l].v[m].v[n].dose += pixelVector[k].v[l].v[m].dose +
                                                          pixelVector[k].v[l].dose +
                                                          pixelVector[k].dose;
                    outLocalDoseFile << setw(15) << pixelVector[k].v[l].v[m].v[n].x
                                     << setw(15) << pixelVector[k].v[l].v[m].v[n].y
                                     << setw(15) << pixelVector[k].v[l].v[m].v[n].dose << endl;
                }
                pixelVector[k].v[l].v[m].dose = 0;
            }
            pixelVector[k].v[l].dose = 0;
        }
        pixelVector[k].dose = 0;
    }
    outLocalDoseFile.close();
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::writeDoses(vector<double> &doses)
{
    cleanNucleus();
    
    int index = 0;
    
    for( int k = 0; k < numberOfBiggestPixels; k++ ) {
        for( int l = 0; l < pixelVector[k].numberOfSubPixels; l++ ) {
            for( int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ ) {
                for( int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ ) {
                    pixelVector[k].v[l].v[m].v[n].dose = doses[index];
                    index++;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------

void Nucleus_Pixel::createPixels()
{
    int i, j, k, l, m, n;
    int numberOfX = static_cast<int>( 2 * ceil(r_nucleus / pixelSide_4) + 1 );

    double x_tmp;
    double y_tmp;
    double *x_4 = new double[ numberOfX * numberOfX ];
    double *y_4 = new double[ numberOfX * numberOfX ];
    double *x_3 = new double[ scale_3 * scale_3 ];
    double *y_3 = new double[ scale_3 * scale_3 ];
    double *x_2 = new double[ scale_2 * scale_2 ];
    double *y_2 = new double[ scale_2 * scale_2 ];
    double *x_1 = new double[ scale_1 * scale_1 ];
    double *y_1 = new double[ scale_1 * scale_1 ];

    k = 0;
    for( i = 0; i < numberOfX; i++ )
    {
        x_tmp = (-(numberOfX - 1)/2 + i) * pixelSide_4;
        for( j = 0; j < numberOfX; j++ )
        {
            y_tmp = (-(numberOfX - 1)/2 + j) * pixelSide_4;
            if( intersection(x_tmp, y_tmp, pixelSide_4) )
            {
                x_4[k] = x_tmp;
                y_4[k] = y_tmp;
                k++;
            }
        }
    }

    numberOfBiggestPixels = k;
    pixelVector = new Pixel[ numberOfBiggestPixels ];
    numberOfSmallestPixels = 0;

    for( k = 0; k < numberOfBiggestPixels; k++)
    {
        pixelVector[k].x = x_4[k];
        pixelVector[k].y = y_4[k];
        pixelVector[k].dose = 0;

        l = 0;
        for( i = 1; i <= scale_3; i++ )
        {
            x_tmp = -(pixelSide_3 + pixelSide_4) / 2 + i * pixelSide_3 + pixelVector[k].x;
            for( j = 1; j <= scale_3; j++ )
            {
                y_tmp = -(pixelSide_3 + pixelSide_4) / 2 + j * pixelSide_3 + pixelVector[k].y;
                if( intersection(x_tmp, y_tmp, pixelSide_3) )
                {
                    x_3[l] = x_tmp;
                    y_3[l] = y_tmp;
                    l++;
                }
            }
        }

        pixelVector[k].numberOfSubPixels = l;
        pixelVector[k].v = new Pixel[ pixelVector[k].numberOfSubPixels ];

        for( l = 0; l < pixelVector[k].numberOfSubPixels; l++)
        {
            pixelVector[k].v[l].x = x_3[l];
            pixelVector[k].v[l].y = y_3[l];
            pixelVector[k].v[l].dose = 0;

            m = 0;
            for( i = 1; i <= scale_2; i++ )
            {
                x_tmp = -(pixelSide_2 + pixelSide_3) / 2 + i * pixelSide_2 + pixelVector[k].v[l].x;
                for( j = 1; j <= scale_2; j++ )
                {
                    y_tmp = -(pixelSide_2 + pixelSide_3) / 2 + j * pixelSide_2 + pixelVector[k].v[l].y;
                    if( sqrt(pow(x_tmp,2) + pow(y_tmp,2)) <= r_nucleus )
                    {
                        x_2[m] = x_tmp;
                        y_2[m] = y_tmp;
                        m++;
                    }
                }
            }

            pixelVector[k].v[l].numberOfSubPixels = m;
            pixelVector[k].v[l].v = new Pixel[ pixelVector[k].v[l].numberOfSubPixels ];

            for( m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++)
            {
                pixelVector[k].v[l].v[m].x = x_2[m];
                pixelVector[k].v[l].v[m].y = y_2[m];
                pixelVector[k].v[l].v[m].dose = 0;
            
                n = 0;
                for( i = 1; i <= scale_1; i++ )
                {
                    x_tmp = -(pixelSide_1 + pixelSide_2) / 2 + i * pixelSide_1 + pixelVector[k].v[l].v[m].x;
                    for( j = 1; j <= scale_1; j++ )
                    {
                        y_tmp = -(pixelSide_1 + pixelSide_2) / 2 + j * pixelSide_1 + pixelVector[k].v[l].v[m].y;
                        if( sqrt(pow(x_tmp,2) + pow(y_tmp,2)) <= r_nucleus )
                        {
                            x_1[n] = x_tmp;
                            y_1[n] = y_tmp;
                            n++;
                        }
                    }
                }

                pixelVector[k].v[l].v[m].numberOfSubPixels = n;
                numberOfSmallestPixels += n;
                pixelVector[k].v[l].v[m].v = new Pixel[ pixelVector[k].v[l].v[m].numberOfSubPixels ];

                for( n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++)
                {
                    pixelVector[k].v[l].v[m].v[n].x = x_1[n];
                    pixelVector[k].v[l].v[m].v[n].y = y_1[n];
                    pixelVector[k].v[l].v[m].v[n].dose = 0;
                    pixelVector[k].v[l].v[m].v[n].numberOfSubPixels = 0;
                    pixelVector[k].v[l].v[m].v[n].v = 0;
                }
            }
        }
    }

    delete [] x_4;
    x_4 = 0;
    delete [] y_4;
    y_4 = 0;
    delete [] x_3;
    x_3 = 0;
    delete [] y_3;
    y_3 = 0;
    delete [] x_2;
    x_2 = 0;
    delete [] y_2;
    y_2 = 0;
    delete [] x_1;
    x_1 = 0;
    delete [] y_1;
    y_1 = 0;

    cout << "Pixel structure created (number of smallest pixels = " << numberOfSmallestPixels << ")... " << endl;
}

//------------------------------------------------------------------------------

inline bool Nucleus_Pixel::intersection(const double x_pixel,
                                        const double y_pixel,
                                        const double pixel_side) const
{
    if( sqrt(pow(x_pixel, 2) + pow(y_pixel, 2)) <= r_nucleus)
        return true;

    if( sqrt(pow(fabs(x_pixel) - pixel_side/2, 2) + pow(fabs(y_pixel) - pixel_side/2, 2)) < r_nucleus )
        return true;

    if( (fabs(x_pixel) < pixel_side/2 && fabs(y_pixel) - pixel_side/2 < r_nucleus) ||
        (fabs(y_pixel) < pixel_side/2 && fabs(x_pixel) - pixel_side/2 < r_nucleus) )
        return true;

    return false;
}

//------------------------------------------------------------------------------

inline bool Nucleus_Pixel::intersection(const double x_pixel,
                                        const double y_pixel,
                                        const double pixel_side,
                                        const double x_track,
                                        const double y_track,
                                        const double radius,
                                        double &distance) const
{
    distance = sqrt(pow(x_pixel - x_track, 2) + pow(y_pixel - y_track, 2));

    if( distance <= radius)
        return true;
    if( sqrt(pow(fabs(x_pixel - x_track) - pixel_side/2, 2) + pow(fabs(y_pixel - y_track) - pixel_side/2, 2)) < radius )
        return true;
    if( (fabs(x_pixel - x_track) < pixel_side/2 && fabs(y_pixel - y_track) - pixel_side/2 < radius) ||
        (fabs(y_pixel - y_track) < pixel_side/2 && fabs(x_pixel - x_track) - pixel_side/2 < radius) )
        return true;
    return false;
}
