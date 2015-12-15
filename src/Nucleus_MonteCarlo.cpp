#include "Nucleus_MonteCarlo.h"
#include "CellLine.h"
#include "Track.h"

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;

using namespace Survival;

Nucleus_MonteCarlo::Nucleus_MonteCarlo(const CellLine &cellLineRef,
                                       const double precision,
                                       const double xPosition,
                                       const double yPosition,
                                       const double pixelSide1,
                                       const int scale1,
                                       const double radius1,
                                       const int scale2,
                                       const double radius2,
                                       const int scale3,
                                       const double radius3)
: Nucleus_Pixel( cellLineRef, xPosition, yPosition,
                 pixelSide1, scale1, radius1,
                 scale2, radius2,
                 scale3, radius3 ),
  distributedTracks( 1, 1.0 )
{
    numberOfIterations = 10000;
    relativeStdDeviation = 1e+6;
    double tmp;
    if( 0.0 < precision && precision < 1.0 )
        relativeStdDeviation = precision;
    else if( 1.0 <= precision && modf(precision, &tmp) == 0.0)
        numberOfIterations = static_cast<int>(precision);
    else
    {
        cerr << "Precision has not been set correctly." << endl;
        exit(1);
    }

    cout << "Creating Nucleus of type MonteCarlo... " << flush;
    cout << "done" << endl
         << endl << flush;
}

//------------------------------------------------------------------------------

void Nucleus_MonteCarlo::cleanNucleus()
{
    Nucleus_Pixel::cleanNucleus();

    distributedTracks.eraseAll();
}

//------------------------------------------------------------------------------

void Nucleus_MonteCarlo::distributeDose(const Track &track)
{
    Nucleus_Pixel::distributeDose(track);

    distributedTracks << track;
}

//------------------------------------------------------------------------------

void Nucleus_MonteCarlo::distributeDose(const Tracks &tracks)
{
    Nucleus_Pixel::distributeDose(tracks);

    distributedTracks << tracks;
}

//------------------------------------------------------------------------------

void Nucleus_MonteCarlo::getDoseAndSurvival(double &dose,
                                            double &doseUncertainty,
                                            double &survival,
                                            double &survivalUncertainty)
{
    addBackgroundDose(1);

    Nucleus_Pixel::getDoseAndSurvival(dose, doseUncertainty, survival, survivalUncertainty);
    cerr << "d = " << dose << "   s = " << survival << endl;
    int i = 0;

    double nucleusFraction = 1 / static_cast<double>(numberOfSmallestPixels);

    double *x = new double[ numberOfSmallestPixels ];
    double *y = new double[ numberOfSmallestPixels ];
    double *probability = new double[ numberOfSmallestPixels ];
    double *cumulative = new double[ numberOfSmallestPixels ];
    double tmp = 0;
    for(int k = 0; k < numberOfBiggestPixels; k++ )
    {
        for(int l = 0; l < pixelVector[k].numberOfSubPixels; l++ )
            for(int m = 0; m < pixelVector[k].v[l].numberOfSubPixels; m++ )
                for(int n = 0; n < pixelVector[k].v[l].v[m].numberOfSubPixels; n++ )
                {
                    probability[i] = pixelVector[k].v[l].v[m].v[n].dose * nucleusFraction / dose;
                    tmp += probability[i];
                    cumulative[i] = tmp;
                    x[i] = pixelVector[k].v[l].v[m].v[n].x;
                    y[i] = pixelVector[k].v[l].v[m].v[n].y;
                    i++;
                }
    }
    gsl_rng *randomGenerator;
    randomGenerator = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(randomGenerator, time(0));

    double random;
    int low, middle, high;
    double x_extracted, y_extracted, x_track, y_track, distance;
    dose = 0;
    double localDose;
    double doseSquared = 0;
    double numberOfLethal = 0, numberOfLethalSquared = 0;
    i = 1;
    do {
        random = gsl_rng_uniform(randomGenerator);
        low = 0;
        high = numberOfSmallestPixels - 1;
        do {
            middle = (low + high) / 2;
            if( cumulative[middle] < random )
                low = middle;
            else
                high = middle;
        } while( high - low > 1 );

        x_extracted = x[high] + pixelSide_1 * ( gsl_rng_uniform(randomGenerator) - 0.5);
        y_extracted = y[high] + pixelSide_1 * ( gsl_rng_uniform(randomGenerator) - 0.5);

        localDose = 0;
        for(int j = 0; j < distributedTracks.size(); j++ )
        {
            distributedTracks[j].getPosition(x_track, y_track);
            x_track = (x_track - x_nucleus) * 1e3;  // mm -> um
            y_track = (y_track - y_nucleus) * 1e3;  // mm -> um
            distance = sqrt(pow(x_extracted-x_track, 2) + pow(y_extracted-y_track, 2));
            localDose += distributedTracks[j].getLocalDose(distance);
        }
        dose += localDose / probability[high] * nucleusFraction;
        doseSquared += pow(localDose / probability[high] * nucleusFraction, 2);

        numberOfLethal += - cellLine.getLogSurvival_X( localDose ) / probability[high] * nucleusFraction;
        numberOfLethalSquared += pow( cellLine.getLogSurvival_X( localDose ) / probability[high] * nucleusFraction, 2 );

        survival = exp(-numberOfLethal/i);
        survivalUncertainty = survival * sqrt( (numberOfLethalSquared / i - pow(numberOfLethal/i, 2)) / (i==1 ? i : (i - 1)) );

        i++;
    } while( i <= numberOfIterations || survivalUncertainty/survival > relativeStdDeviation );

    dose /= --i;
    doseUncertainty = sqrt( (doseSquared / i - pow(dose, 2)) / (i - 1) );

    gsl_rng_free(randomGenerator);

    delete [] x;
    delete [] y;
    delete [] probability;
    delete [] cumulative;
}
