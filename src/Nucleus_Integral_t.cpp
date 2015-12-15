#include "Nucleus_Integral_t.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Track.h"

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::ceil;
using std::acos;
using std::min;

using std::string;
using std::vector;

using namespace Survival;

Nucleus_Integral_t::Nucleus_Integral_t(const CellLine &cellLineRef,
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
/*  Questo è scritto bene e compila
Nucleus_Integral_t::Nucleus_Integral_t(const Nucleus_Integral_t& nn,
                                       const CellLine &cellLineRef)
: totalNucleusDose( nn.totalNucleusDose ),
  cellLine( cellLineRef ),
  r_nucleus( nn.r_nucleus ),
  x_nucleus( nn.x_nucleus ),
  y_nucleus( nn.y_nucleus ),
  inNucleusCount( nn.inNucleusCount ),
  intersectionCount( nn.intersectionCount )
{
    times = nn.times;
    doses = nn.doses;
}
*/

//------------------------------------------------------------------------------

void Nucleus_Integral_t::addBackgroundDose(const double dose,   // Gy
                                           const double t)
{
    totalNucleusDose += dose;
    doses.push_back(dose);
    times.push_back(t);
}

//------------------------------------------------------------------------------

double Nucleus_Integral_t::ArcIntersectionWeight(double r,
                                                 double b)
{
    double arg;
    
    if (b < r_nucleus) {
        if (r <= r_nucleus - b) {
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

void Nucleus_Integral_t::cleanNucleus()
{
    totalNucleusDose = 0.0;
    inNucleusCount = 0;
    intersectionCount = 0;
    times.clear();
    doses.clear();
}

//------------------------------------------------------------------------------
// NOTA: Non usare ... buggy ...
Nucleus_Integral_t* Nucleus_Integral_t::clone(const CellLine& cellLine)
{
    Nucleus_Integral_t* nucleo=new Nucleus_Integral_t(cellLine); //VA CANCELLATO DA QUALCHE PARTE!!!
    return nucleo;
}

//------------------------------------------------------------------------------

void Nucleus_Integral_t::distributeDose(const Track &track)
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
            } else {
                area3 = M_PI * (r_nucleus*r_nucleus - r_intersection*r_intersection);
            }
        }
        
        dose  /= area1 + area2 + area3;
        
        totalNucleusDose += dose;
        doses.push_back(dose);
        times.push_back(track.getTime());
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
        doses.push_back(dose);
        times.push_back(track.getTime());
    } else
        return;
}

//------------------------------------------------------------------------------

void Nucleus_Integral_t::distributeDose(const Tracks &tracks)
{
    for( int i = 0; i < tracks.size(); i++ )
        Nucleus_Integral_t::distributeDose(tracks[i]);
}

//------------------------------------------------------------------------------

string Nucleus_Integral_t::getCellType() const
{
    return cellLine.getCellType();
}

//------------------------------------------------------------------------------

void Nucleus_Integral_t::getDoseAndSurvival(double &dose,
                                            double &doseUncertainty,
                                            double &survival,
                                            double &survivalUncertainty) const
{
    dose = totalNucleusDose;
    double sum_lethal = cellLine.getLogSurvival_X( doses, times );
    survival = exp(-sum_lethal);
    doseUncertainty = -1;
    survivalUncertainty = -1;
}

//------------------------------------------------------------------------------

void Nucleus_Integral_t::getDoseAndLethals(double &dose,
                                           double &doseUncertainty,
                                           double &lethals,
                                           double &lethalsUncertainty) const
{
    dose = totalNucleusDose;
    lethals = cellLine.getLogSurvival_X( doses, times );
    doseUncertainty = -1;
    lethalsUncertainty = -1;
}

//------------------------------------------------------------------------------

void Nucleus_Integral_t::getPosition(double &returnX,
                                     double &returnY) const
{
    returnX = x_nucleus;
    returnY = y_nucleus;
}

//------------------------------------------------------------------------------

double Nucleus_Integral_t::IntegrateWeightedRadialTrack(const Track &track,
                                                        double rMin,
                                                        double rMax,
                                                        double b,
                                                        double &area,
                                                        double step)
{
    double r1, r2, log_r2, log_rMin, log_rMax, log_step;
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

    arc_weight2 = Nucleus_Integral_t::ArcIntersectionWeight(rMin, b);
	f2 = track.getLocalDose(rMin) * rMin * arc_weight2;
	r2 = rMin;
	for (int i = 0; i < nSteps -1; i++) {   // morale: lo divide in step, calcola il valore nei vari step, e prende la media in mezzo...
                                            //... in pratica una funzione a gradini!
		log_r2 = log_rMin + log_step*(i + 1);
		f1 = f2;
		r1 = r2;
		arc_weight1 = arc_weight2;

		r2 = pow(10, log_r2);
		arc_weight2 = Nucleus_Integral_t::ArcIntersectionWeight(r2, b);
		f2 = track.getLocalDose(r2) * r2 * arc_weight2;
		f =  (r2 - r1)* (f1/2. + f2/2.);
		integral += f;
		area += (r2 - r1)* (arc_weight1 * r1 /2. + arc_weight2 * r2 /2.);
	}
	f1 = f2;
	r1 = r2;
	arc_weight1 = arc_weight2;

	r2 = rMax;
	arc_weight2 = Nucleus_Integral_t::ArcIntersectionWeight(r2, b);
	f2 = track.getLocalDose(r2) * r2 * arc_weight2;
	f =  (r2 - r1)* (f1/2. + f2/2.);
	integral += f;
	area += (r2 - r1) * (arc_weight1 * r1/2. + arc_weight2 * r2/2.);
	
	return integral/area;
}
