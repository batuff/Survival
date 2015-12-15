#include "Nucleus_Integral.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Track.h"

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::ceil;
using std::acos;
using std::min;

#include <iostream>
using std::cout;
using std::endl;

using std::string;
using std::vector;

using namespace Survival;

Nucleus_Integral::Nucleus_Integral(const CellLine &cellLineRef,
                             const double xPosition,
                             const double yPosition)
: cellLine( cellLineRef ),
  x_nucleus( xPosition ),
  y_nucleus( yPosition )
{
    r_nucleus = cellLine.getNucleusRadius();

    cleanNucleus();
}

//------------------------------------------------------------------------------

void Nucleus_Integral::addBackgroundDose(const double dose)  // Gy
{
    totalNucleusDose += dose;
}

//------------------------------------------------------------------------------

void Nucleus_Integral::cleanNucleus()
{
    totalNucleusDose = 0;
    inNucleusCount = 0;
    intersectionCount = 0;
}

//------------------------------------------------------------------------------

Nucleus_Integral* Nucleus_Integral::clone(const CellLine& cellLine)
{
    Nucleus_Integral* nucleo=new Nucleus_Integral(cellLine); //VA CANCELLATO DA QUALCHE PARTE!!!
    return nucleo;
}

/*          vecchia
 Nucleus_Integral* Nucleus_Integral::clone() const
 {
 return new Nucleus_Integral(*this);
 }
 */

//------------------------------------------------------------------------------

void Nucleus_Integral::distributeDose(const Track &track)
{
    double x_track;
    double y_track;
    double b, area1, area2, area3, theta1, theta2, dose = 0.0;
    double rMin, rMax, r_intersection;
    
    track.getPosition(x_track, y_track);
    x_track = (x_track - x_nucleus) * 1e3;  // mm -> um
    y_track = (y_track - y_nucleus) * 1e3;  // mm -> um
    b = sqrt( x_track*x_track + y_track*y_track );
    
    rMax = min(track.getRadius(), b + r_nucleus);
    
    area1 = area2 = area3 = 0.0;
    
    if( b <= r_nucleus ) {
        inNucleusCount++;
        rMin = 0.;
        
        if (b + track.getRadius() < r_nucleus)
            r_intersection = track.getRadius();
        else
            r_intersection = r_nucleus - b;
        
        area1 = M_PI * r_intersection * r_intersection;
        dose  += track.getRadialIntegral(0.0, r_intersection) * area1;
        
        if (rMax > r_intersection)
            dose += IntegrateWeightedRadialTrack(track, r_intersection, rMax, b, area2, 0.01) * area2;
        
        if (rMax == track.getRadius()) {
            if (track.getRadius() > r_nucleus - b) {
                theta1 = acos( (b*b + rMax*rMax - r_nucleus*r_nucleus)/(2*b*rMax) );
                theta2 = acos( (b*b - rMax*rMax + r_nucleus*r_nucleus)/(2*b*r_nucleus) );
                area3 = M_PI * (r_nucleus*r_nucleus) - (theta1*rMax*rMax + theta2*r_nucleus*r_nucleus - rMax*b*sin(theta1));
            } else
                area3 = M_PI * (r_nucleus*r_nucleus - r_intersection*r_intersection);
        }
        
        dose  /= area1 + area2 + area3;
        
        totalNucleusDose += dose;
        
    } else if( b <= r_nucleus + track.getRadius() ) {
        intersectionCount++;
        rMin = b - r_nucleus;
        dose = IntegrateWeightedRadialTrack(track, rMin, rMax, b, area2, 0.01) * area2;
        
        if (rMax == track.getRadius()) {
            theta1 = acos( (b*b + rMax*rMax - r_nucleus*r_nucleus)/(2*b*rMax) );
            theta2 = acos( (b*b - rMax*rMax + r_nucleus*r_nucleus)/(2*b*r_nucleus) );
            area3 = M_PI * (r_nucleus*r_nucleus) - (theta1*rMax*rMax + theta2*r_nucleus*r_nucleus - rMax*b*sin(theta1));
        }
        dose /= area2 + area3;
        totalNucleusDose += dose;
    } else
        return;
}

//------------------------------------------------------------------------------

void Nucleus_Integral::distributeDose(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        Nucleus_Integral::distributeDose(tracks[i]);
}

//------------------------------------------------------------------------------

string Nucleus_Integral::getCellType() const
{
    return cellLine.getCellType();
}

//------------------------------------------------------------------------------

void Nucleus_Integral::getDoseAndSurvival(double &dose,
                                          double &doseUncertainty,
                                          double &survival,
                                          double &survivalUncertainty) const
{
    dose = totalNucleusDose;
    double sum_lethal = cellLine.getLogSurvival_X( totalNucleusDose );
    survival = exp(-sum_lethal);
    doseUncertainty = -1;
    survivalUncertainty = -1;
}

//------------------------------------------------------------------------------

void Nucleus_Integral::getPosition(double &returnX,
                                   double &returnY) const
{
    returnX = x_nucleus;
    returnY = y_nucleus;
}

//------------------------------------------------------------------------------

double Nucleus_Integral::ArcIntersectionWeight(double r, double b)
{
    double arg;
    
    if (b < r_nucleus) {
        if (r <= r_nucleus - b){
            return 2 * M_PI;
        } else if (r < b + r_nucleus) {
            arg = b/(2*r) + r/(2*b) - r_nucleus*r_nucleus / (2*b*r);
            return 2 * acos( arg );
        }
    } else {
        if (r <= b - r_nucleus) {
            return 0.;
        } else if (r < b + r_nucleus) {
            arg = b/(2*r) + r/(2*b) - r_nucleus*r_nucleus / (2*b*r);
            if (arg > 1.)
                arg = 1.;
            return 2 * acos( arg );
        }
    }
    return 0.;
}

//------------------------------------------------------------------------------

/*
 double Nucleus_Integral::IntegrateRadialTrack(const Track &track, double rMin, double rMax, double step)
 {
 double r1, r2;
	double f1, f2;
	double integral = 0.;
	int nSteps = (int)ceil( (rMax-rMin)/step );
 
	f2 = track.getLocalDose(rMin) * rMin;
	for (int i = 0; i < nSteps; i++) {
 r1 = rMin + step*i;
 r2 = rMin + step*(i + 1);
 f1 = f2;
 f2 = track.getLocalDose(r2) * r2;
 integral += step * (f1/2. + f2/2.);
	}
 
 return 2*M_PI*integral;
 }
 */

//------------------------------------------------------------------------------

double Nucleus_Integral::IntegrateWeightedRadialTrack(const Track &track, double rMin, double rMax, double b, double &area, double step)
{
    double r1, r2, log_r2, log_rMin, log_rMax, log_step;
    //double log_r1;
    double f1, f2, f, arc_weight1, arc_weight2;
    double integral = 0.;
    if (rMin > 0)
        log_rMin = log10(rMin);
    else
        log_rMin = -5;
    
    log_rMax = log10(rMax);
    log_step = step;
    
    int nSteps = (int)ceil( (log_rMax - log_rMin)/log_step );
    
    if (nSteps < 3) {
        log_step = (log_rMax - log_rMin)/3;
        nSteps = 3;
    }
    
    area = 0;
    
    /*
     f2 = track.getLocalDose(rMin) * rMin * Nucleus_Integral::ArcIntersectionWeight(rMin, b);
     //if (isnan(f2)) cout << "isnan!1" << endl;
     for (int i = 0; i < nSteps; i++) {
     //r1 = rMin + step*i;
     r2 = rMin + step*(i + 1);
     f1 = f2;
     f2 = track.getLocalDose(r2) * r2 * Nucleus_Integral::ArcIntersectionWeight(r2, b);
     //if (isnan(f2)) cout << "isnan!2" << endl;
     f = step * (f1/2. + f2/2.);
     integral += f;
     }
     */
    
    arc_weight2 = Nucleus_Integral::ArcIntersectionWeight(rMin, b);
    f2 = track.getLocalDose(rMin) * rMin * arc_weight2;
    r2 = rMin;
    for (int i = 0; i < nSteps -1; i++) {   // morale: lo divide in step, calcola il valore nei vari step, e prende la media in mezzo...
        //... in pratica una funzione a gradini!
        log_r2 = log_rMin + log_step*(i + 1);
        f1 = f2;
        r1 = r2;
        arc_weight1 = arc_weight2;
        
        r2 = pow(10, log_r2);
        arc_weight2 = Nucleus_Integral::ArcIntersectionWeight(r2, b);
        f2 = track.getLocalDose(r2) * r2 * arc_weight2;
        f =  (r2 - r1)* (f1/2. + f2/2.);
        integral += f;
        area += (r2 - r1)* (arc_weight1 * r1 /2. + arc_weight2 * r2 /2.);
    }
    f1 = f2;
    r1 = r2;
    arc_weight1 = arc_weight2;
    
    r2 = rMax;
    arc_weight2 = Nucleus_Integral::ArcIntersectionWeight(r2, b);
    f2 = track.getLocalDose(r2) * r2 * arc_weight2;
    f =  (r2 - r1)* (f1/2. + f2/2.);
    integral += f;
    area += (r2 - r1) * (arc_weight1 * r1/2. + arc_weight2 * r2/2.);
    
    return integral/area;
}
