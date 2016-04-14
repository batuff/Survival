#include "Calculus.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Track.h"
#include "Nucleus_Integral.h"
#include "Nucleus_Integral_t.h"
#include "Nucleus_tMKM.h"
#include "usefulFunctions.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_gamma.h>

#include <sstream>
using std::ostringstream;

#include <fstream>
using std::ofstream;

#include <cstdlib>
using std::exit;

#include <ctime>
using std::time;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::abs;
#ifdef OSX
using std::isnan;
#endif

#include <iomanip>
using std::setprecision;
using std::setw;

#include <iostream>
using std::cout;
using std::clog;
using std::cerr;
using std::ios;
using std::endl;
using std::fixed;
using std::scientific;
using std::left;

#include <algorithm>
using std::sort;

using std::string;
using std::vector;

#include <sys/types.h>
#include <unistd.h>
#include <omp.h>

using namespace Survival;

Calculus::Calculus(const Tracks &tracksRef,
                   const CellLine &cellLineRef,
                   Nucleus &nucleusRef,
                   string save_prefix,
                   string model_,
                   int p_type,
                   long randomSeed)
: tracks( tracksRef ),
  cellLine( cellLineRef ),
  nucleus( nucleusRef ),
  savePrefix( save_prefix ),
  model( model_ )
{
    randomGenerator = gsl_rng_alloc(gsl_rng_taus);
    long seed = randomSeed;
    cout << "Calculus object created -- ";
    if (randomSeed == 0) {
        seed = time(NULL) * getpid();
        cout << "Setting random seed: " << seed << endl;
    } else
        cout << "Using external random seed: " << seed << endl;
    
    if (p_type!=0)
        nThreads=p_type;
    else
        nThreads=omp_get_num_procs();
    
    gsl_rng_set(randomGenerator, seed);
}

//------------------------------------------------------------------------------

Calculus::~Calculus()
{
   gsl_rng_free(randomGenerator);
}

//------------------------------------------------------------------------------

void Calculus::evaluateG(const string trackMode,
                         double totalDose,
                         int nFrac,
                         double timeSpacing,
                         double fracDeliveryTime,
                         double precision,
                         double alpha,
                         double beta)
{
    clog << endl << "------------------------------------------------------------------------------------" << endl << endl
         << "Evaluating G (Lea-Catcheside) factor for fractionating treatment:" << endl << endl
         << "    Total dose: " << totalDose << " Gy" << endl
         << "    Number of fraction: " << nFrac << endl
         << "    Time spacing between fraction: " << timeSpacing << " hours" << endl
         << "    Fraction delivery time: " << fracDeliveryTime << " hours" << endl;
    
    bool num=false, dev=false, abort=false;
    
    if (trackMode != "random")
    {
        cerr << "EvaluateG: invalid trackMode. Only \"histogram\" supported." << endl;
        exit(1);
    }
    
    ostringstream GFile;
    GFile << savePrefix << setprecision(2) << fixed << timeSpacing << ".txt";
    ofstream outGFile(GFile.str().c_str(), ios::app);
    if(!outGFile)
    {
        cerr << "File " << GFile.str() << " could not be opened." << endl;
        exit(1);
    }
    outGFile << setprecision(6) << scientific << left;
    
    int numberOfIterations = 10;
    double relativeStdDeviation = 1e+6;
    double dump;
    if( 0.0 < precision && precision < 1.0 )
        relativeStdDeviation = precision;
    else if( 1.0 <= precision && modf(precision, &dump) == 0.0){
        numberOfIterations = static_cast<int>(precision);
        relativeStdDeviation=0.0;
        dev=true;
    }
    else
    {
        cerr << "Precision has not been set correctly." << endl;
        exit(1);
    }

    omp_set_num_threads(nThreads);
    
    vector<Nucleus*> nucleus_vec(nThreads);
    for (int it=0; it<nThreads; ++it)
        nucleus_vec[it]=nucleus.clone(cellLine);
    
    double doseFrac = totalDose/(double)nFrac, time=0.0;
    double dose=0.0, meanDose=0.0, survival=0.0, meanSurvival=0.0, tmp1=0.0, tmp2=0.0;
    int it=0, thread=0;
    double event=0.0;
    double sumSSquared=0.0, varS=0.0, meanSUncertainty=0.0;
    
    ostringstream SFile;
    SFile   << savePrefix << nFrac << "_" << setprecision(2) << fixed << timeSpacing << "_"
            << setprecision(2) << fixed << fracDeliveryTime << ".dat";
    ofstream outSFile(SFile.str().c_str(), ios::app);
    if(!outSFile)
    {
        cerr << "File " << SFile.str() << " could not be opened." << endl;
        exit(1);
    }
    outSFile << setprecision(6) << scientific << left;

    #pragma omp parallel private(dose, survival, it, time, thread)
    {
        #pragma omp critical
        {
            clog << "Thread " << omp_get_thread_num() << " entering parallel loop" << endl;
        }
        while ((!num)||(!dev)) {
            thread = omp_get_thread_num();
            for (it=0; it<nFrac; ++it) {
                time = timeSpacing * it;
                histogram_dose_survival_t(*nucleus_vec[thread], doseFrac, time, fracDeliveryTime);
            }
            nucleus_vec[thread]->getDoseAndSurvival(dose,tmp1,survival,tmp2);
            
            #pragma omp critical
            {
                if(!isnan(dose) && !isnan(survival) && !abort){
                    event+=1.0;
                    
                    outSFile << setprecision(8) << fixed << setw(15) << dose << setw(15) << survival << endl;
                    
                    meanDose = (meanDose * (event-1.0) + dose)/event;
                    meanSurvival = (meanSurvival * (event-1.0) + survival)/event;
                    
                    sumSSquared += survival*survival;
                    varS = sumSSquared / event - meanSurvival*meanSurvival;  //biased estimation
                    meanSUncertainty = sqrt(varS / event);
                    
                    if (event >= numberOfIterations)
                        num=true;
                    if (num && meanSUncertainty/meanSurvival <= relativeStdDeviation && meanDose>1e-4)
                        dev=true;
                    
                    cout << (int)event << ")"
                         << "   Mean dose = " << setprecision(8) << meanDose
                         << "   Mean survival = " << setprecision(8) << meanSurvival
                         << "   Precision = " << setprecision(2) << meanSUncertainty/meanSurvival*100 << "%" << endl;
                    
                    if(num && dev) abort=true;
                }//isnan
            }//critical

            nucleus_vec[thread]->cleanNucleus();
        }//while
        #pragma omp critical
        {
            clog << "Thread " << omp_get_thread_num() << " exiting parallel loop" << endl;
        }
    }//parallel
    
    outSFile.close();
    
    double G_sim = (-log(meanSurvival)-alpha*totalDose)/(beta*totalDose*totalDose);
    double G_macro = 1.0;
    for (int i=0; i<nFrac-1; ++i)
        for (int j=i+1; j<nFrac; ++j)
            G_macro -= 2/(totalDose*totalDose) * ( 1-exp(-cellLine.getAC()*timeSpacing*(j-i)) ) * doseFrac*doseFrac;
    
    outGFile << setw(8) << fixed << nFrac << setw(12) << setprecision(3) << fracDeliveryTime
             << setw(12) << setprecision(4) << fixed << G_sim << setprecision(4) << fixed << G_macro << endl;
    outGFile.close();
    
    cout << endl
         << setprecision(6) << fixed << "Total mean dose = " << meanDose
         << "   Total mean survival = " << meanSurvival
         << "   G (simulated) = " << G_sim
         << "   G (macroscopical) = " << G_macro
         << endl << endl;
    
    for (int iter=0; iter<nThreads; ++iter)
        delete nucleus_vec[iter];
}

//------------------------------------------------------------------------------

vector<double> Calculus::generateSequence(int nEv,
                                          double begin,
                                          double duration)
{
    if(duration<0){
        cerr << "Timeline error in function \"generateSequence\"." << endl;
        exit(1);
    }
    vector<double> tt(nEv);
    for (int i=0; i<nEv; ++i)
        tt[i] = begin + duration * gsl_rng_uniform(randomGenerator);
    sort(tt.begin(), tt.end());
    return tt;
}

//------------------------------------------------------------------------------

void Calculus::histogram_dose_survival_p(const double doseImposed, // implementata per il calcolo parallelo
                                         double &dose,
                                         double &doseUncertainty,
                                         double &survival,
                                         double &survivalUncertainty,
                                         Nucleus &nucleus_cp,       // gli si passa anche il nucleo su cui lavorare
                                         bool clean)
{
    Track *extractedTrack;
    int numberOfExtractions;
    double rho_track, phi_track;
    double x_track, y_track;
    
    double meanLet = tracks.getMeanLet();
    double totalWeight = tracks.getTotalWeight();
    
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    double totalFluence = tracks.getDensity() * doseImposed / (CONV * meanLet);  //um^(-2)
    double r_nucleus = nucleus_cp.getRadius();
    double area;  //um^2
    double meanNumber;

    for(int k = 0; k < tracks.size(); k++)
    {
        area = M_PI * pow(tracks[k].getRadius() + r_nucleus, 2);
        meanNumber = totalFluence * tracks[k].getWeight()/totalWeight * area;
        numberOfExtractions = gsl_ran_poisson(randomGenerator, meanNumber);
        
        for(int j = 0; j < numberOfExtractions; j++)
        {
            rho_track = sqrt( area / M_PI ) * sqrt( gsl_rng_uniform(randomGenerator) );
            phi_track = 2 * M_PI * gsl_rng_uniform(randomGenerator);
            x_track = rho_track * cos(phi_track);
            y_track = rho_track * sin(phi_track);
            
            extractedTrack = tracks[k].clone();
            extractedTrack->setPosition(x_track * 1e-3, y_track * 1e-3); // um -> mm
            extractedTrack->setTime(0.0);
     
            nucleus_cp.distributeDose(*extractedTrack);
            delete extractedTrack;
        }
    }
    
    nucleus_cp.getDoseAndSurvival(dose, doseUncertainty, survival, survivalUncertainty);
    
    if (clean) nucleus_cp.cleanNucleus();
}

//------------------------------------------------------------------------------

void Calculus::histogram_dose_survival_t(Nucleus &nucleus_cp,
                                         const double doseImposed,
                                         const double time,
                                         const double fracDeliveryTime)
{
    Track *extractedTrack;
    vector<int> numberOfExtractions(tracks.size());
    int totTracks=0;
    double rho_track, phi_track;
    double x_track, y_track;
    
    double meanLet = tracks.getMeanLet();
    double totalWeight = tracks.getTotalWeight();
    
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    double totalFluence = tracks.getDensity() * doseImposed / (CONV * meanLet);  //um^(-2)
    double r_nucleus = nucleus_cp.getRadius();
    double area=0.0;  //um^2
    double meanNumber;
    
    for(int k = 0; k < tracks.size(); k++) {
        area = M_PI * pow(tracks[k].getRadius() + r_nucleus, 2);
        meanNumber = totalFluence * tracks[k].getWeight()/totalWeight * area;
        numberOfExtractions[k] = gsl_ran_poisson(randomGenerator, meanNumber);
        totTracks+=numberOfExtractions[k];
    }
    
    vector<double> TIMES = generateSequence(totTracks,time,fracDeliveryTime);
    int it=0;
    for(int k = 0; k < tracks.size(); k++) {
        for(int j = 0; j < numberOfExtractions[k]; j++) {
            rho_track = sqrt( area / M_PI ) * sqrt( gsl_rng_uniform(randomGenerator) );
            phi_track = 2 * M_PI * gsl_rng_uniform(randomGenerator);
            x_track = rho_track * cos(phi_track);
            y_track = rho_track * sin(phi_track);
            
            extractedTrack = tracks[k].clone();
            extractedTrack->setPosition(x_track * 1e-3, y_track * 1e-3); // um -> mm
            extractedTrack->setTime(TIMES[it]);
            
            nucleus_cp.distributeDose(*extractedTrack);
            
            delete extractedTrack;
            ++it;
        }
    }
}

//------------------------------------------------------------------------------

void Calculus::histogram_dose_survival_with_domains(const double doseImposed,
                                                    double &dose,
                                                    double &doseUncertainty,
                                                    double &survival,
                                                    double &survivalUncertainty,
                                                    vector<double> &doses,
                                                    vector<double> &lethals,
                                                    vector<double> &dosesUncertainty,
                                                    vector<double> &lethalsUncertainty,
                                                    bool clean)
{
    Track *extractedTrack;
    int numberOfExtractions;
    double rho_track, phi_track;
    double x_track, y_track;
    
    double meanLet = tracks.getMeanLet();
    double totalWeight = tracks.getTotalWeight();
    
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    double totalFluence = tracks.getDensity() * doseImposed / (CONV * meanLet);  //um^(-2)
    double r_nucleus = nucleus.getRadius();
    double area;  //um^2
    double meanNumber;
    
    for(int k = 0; k < tracks.size(); k++)
    {
        area = M_PI * pow(tracks[k].getRadius() + r_nucleus, 2);
        meanNumber = totalFluence * tracks[k].getWeight()/totalWeight * area;
        numberOfExtractions = gsl_ran_poisson(randomGenerator, meanNumber);
        
        for(int j = 0; j < numberOfExtractions; j++)
        {
            rho_track = sqrt( area / M_PI ) * sqrt( gsl_rng_uniform(randomGenerator) );
            phi_track = 2 * M_PI * gsl_rng_uniform(randomGenerator);
            x_track = rho_track * cos(phi_track);
            y_track = rho_track * sin(phi_track);
            
            extractedTrack = tracks[k].clone();
            extractedTrack->setPosition(x_track * 1e-3, y_track * 1e-3); // um -> mm
            nucleus.distributeDose(*extractedTrack);
            delete extractedTrack;
        }
    }
    
    nucleus.getDoseAndSurvival(dose, doseUncertainty, survival, survivalUncertainty);
    
    nucleus.getDosesAndLethals(doses, dosesUncertainty, lethals, lethalsUncertainty);
    
    if (clean) nucleus.cleanNucleus();
}

//------------------------------------------------------------------------------

void Calculus::random_dose_survival_p(const double doseImposed,
                                    double &dose,
                                    double &doseUncertainty,
                                    double &survival,
                                    double &survivalUncertainty,
                                    Nucleus &nucleus_cp,
                                    bool clean)
{
    double let_tot = 0.0;
    double max_r_max = 0.0;
    for(int k = 0; k < tracks.size(); k++)
    {
        let_tot += tracks[k].getLet();
        if( max_r_max < tracks[k].getRadius() )
            max_r_max = tracks[k].getRadius();
    }
    double meanLet = let_tot / tracks.size();
    double interactionRadius = nucleus_cp.getRadius() + max_r_max;
    
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    double totalFluence = tracks.getDensity() * doseImposed / (CONV * meanLet);  //um^(-2)
    double meanNumber = totalFluence * M_PI * pow(interactionRadius, 2);
    int numberOfExtractions = gsl_ran_poisson(randomGenerator, meanNumber);
    
    int extractedTrackIndex;
    Track *extractedTrack;
    double rho_track, phi_track;
    double x_track, y_track;
    
    for(int k = 0; k < numberOfExtractions; k++)
    {
        extractedTrackIndex = gsl_rng_get(randomGenerator) % tracks.size();
        rho_track = interactionRadius * sqrt( gsl_rng_uniform(randomGenerator) );
        if( rho_track < nucleus_cp.getRadius() + tracks[extractedTrackIndex].getRadius() )
        {
            phi_track = 2 * M_PI * gsl_rng_uniform(randomGenerator);
            x_track = rho_track * cos(phi_track);
            y_track = rho_track * sin(phi_track);
            
            extractedTrack = tracks[extractedTrackIndex].clone();
            extractedTrack->setPosition(x_track * 1e-3, y_track * 1e-3); // um -> mm
            nucleus.distributeDose(*extractedTrack);
            delete extractedTrack;
        }
    }
    
    nucleus_cp.getDoseAndSurvival(dose, doseUncertainty, survival, survivalUncertainty);
    
    cout << "Tracks in nucleus: " << setw(4) << nucleus_cp.getInNucleusCount()
    << "   tracks intersecting nucleus: " << setw(6) << nucleus_cp.getIntersectionCount()
    << "     dose: " << setprecision(6) << dose
    << " (" << setprecision(2) << doseUncertainty/dose * 100.0 << "%)"
    << "     survival: " << setprecision(6) << survival
    << " (" << setprecision(2) << survivalUncertainty/survival * 100.0 << "%)" << endl;
    
    if (clean) nucleus_cp.cleanNucleus();
}

//------------------------------------------------------------------------------

void Calculus::rapidINFN_alphaIon_betaIon(double &alphaIon,
                                          double &betaIon)
{
    // Questo dev'essere quello che nella tesi di Germano viene identificato con "INFN approximate implementation"
    if(true)
    {
        //dummy initialization to avoid warnings when compiling.
        alphaIon=0.0;
        betaIon=0.0;
        cerr << "rapidINFN_alphaIon_betaIon() method hasn't been implemented yet." << endl;
        exit(1);
    }
    /*
    double alpha_X, beta_X, D_t;
    cellLine.getParameters(alpha_X, beta_X, D_t);
    double s = alpha_X + 2 * beta_X * D_t;  //Gy^(-1)
 
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
 
    double dummy1, dummy2;
    double doseSingleImpact;
    double S_1;
    double alpha_D, beta_D;
    double alpha_tot = 0, sqrt_beta_tot = 0;
    double meanLet = 0;
 
    for( int k = 0; k < tracks.size(); k++ )
    {
        Track *track = tracks[k].clone();
        double r_t = track->getDistance(D_t);
        if( r_t < 0. ) r_t = 0.;
        double r_max = track->getRadius();
 
        double d1_lessThanDt = 0.;
        double d1_squared_lessThanDt = 0.;
        if(r_t < r_max)
        {
            //Compute d1_lessThanDt
            double r1;
            double r2 = r_t;
            double log10_r2 = (r2 != 0.) ? log10(r2) : -8.;
            double d2 = track->getLocalDose(r2);
            double i1;
            double i2 = d2 * r2;
            do
            {
                r1 = r2;
                i1 = i2;
                log10_r2 += 5e-4;
                r2 = pow(10.0, log10_r2);
                d2 = track->getLocalDose(r2);
                i2 = d2 * r2;
                d1_lessThanDt += (i1 + i2) * (r2 - r1) / 2;
            } while( r2 < r_max );
            d1_lessThanDt /= (r_max*r_max - r_t*r_t)/2;
 
            //Compute d1_squared_lessThanDt
            r2 = r_t;
            log10_r2 = (r2 != 0.) ? log10(r2) : -8.;
            d2 = track->getLocalDose(r2);
            i2 = d2*d2 * r2;
            do
            {
                r1 = r2;
                i1 = i2;
                log10_r2 += 5e-4;
                r2 = pow(10.0, log10_r2);
                d2 = track->getLocalDose(r2);
                i2 = d2*d2 * r2;
                d1_squared_lessThanDt += (i1 + i2) * (r2 - r1) / 2;
            } while( r2 < r_max );
            d1_squared_lessThanDt /= (r_max*r_max - r_t*r_t)/2;
        }
    }
     */
}

//------------------------------------------------------------------------------

void Calculus::rapidMKM_Hawkins_alphaIon_betaIon(double &alphaIon,
                                                 double &betaIon)
{
    double dummy1, dummy2, dummy3;
    double alpha_X, beta_X;
    cellLine.getParameters(alpha_X, beta_X, dummy1);

    double alpha_D, beta_D;
    double alpha_tot = 0.0, sqrt_beta_tot = 0.0;
    double meanLet = 0.0;
    double r_domain = cellLine.getDomainRadius();
    double r_nucleus = cellLine.getNucleusRadius();
    
    CellLine fakeD("domain", r_domain);
    fakeD.addParametrization_LQ_noDt(alpha_X, beta_X);
    fakeD.setParametrization("LQ_noDt");
    Nucleus_Integral domain(fakeD);
    
    CellLine fakeN("nucleus", r_nucleus);
    fakeN.addParametrization_LQ_noDt(alpha_X, beta_X);
    fakeN.setParametrization("LQ_noDt");
    Nucleus_Integral nucleus_int(fakeN);
    
    for(int k = 0; k < tracks.size(); k++)
    {
        Track *track = tracks[k].clone();
        double r_max = track->getRadius();
        double let = tracks[k].getLet();
        double weight = tracks[k].getWeight();
        
        //Compute gamma
        double r1;
        double r2 = 1e-4;
        double log10_r2 = log10(r2);
        track->setPosition(0.0, r2*1e-3);
        domain.distributeDose(*track);
        double d1, d2;
        domain.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
        domain.cleanNucleus();
        
        double z_1 = 0.;
        double z_1D = 0.;
        do
        {
            r1 = r2;
            d1 = d2;
            log10_r2 += 1e-3; // orig 1e-2
            r2 = pow(10.0, log10_r2);
            
            track->setPosition(0.0, r2*1e-3);
            domain.distributeDose(*track);
            domain.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
            domain.cleanNucleus();
            z_1 += (d1*r1 + d2*r2) * (r2 - r1);
            z_1D += (d1*d1*r1 + d2*d2*r2) * (r2 - r1);
        } while( r2 < r_domain + r_max );
        double gamma = z_1D / z_1;
        
        //Compute gamma_nucleus
        r2 = 1e-4;
        log10_r2 = log10(r2);
        track->setPosition(0.0, r2*1e-3);
        nucleus_int.distributeDose(*track);
        nucleus_int.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
        nucleus_int.cleanNucleus();
        
        z_1 = 0.;
        z_1D = 0.;
        do
        {
            r1 = r2;
            d1 = d2;
            log10_r2 += 1e-3; // orig 1e-2
            r2 = pow(10.0, log10_r2);
            
            track->setPosition(0.0, r2*1e-3);
            nucleus_int.distributeDose(*track);
            nucleus_int.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
            nucleus_int.cleanNucleus();
            z_1 += (d1*r1 + d2*r2) * (r2 - r1);
            z_1D += (d1*d1*r1 + d2*d2*r2) * (r2 - r1);
        } while( r2 < r_nucleus + r_max );
        double gamma_nucleus = z_1D / z_1;
        
        double alpha_P = alpha_X + beta_X * gamma;
        double beta_P = beta_X;
        alpha_D = (1 - exp( - alpha_P * gamma_nucleus )) / gamma_nucleus;
        beta_D = beta_P;
        
        alpha_tot += weight * let * alpha_D;
        sqrt_beta_tot += weight * let * sqrt(beta_D);
        meanLet += weight * let;
        
    }
    
    alphaIon = alpha_tot / meanLet;
    betaIon = pow(sqrt_beta_tot/meanLet, 2);
    
    cout << setprecision(6)
         << "alpha ion = " << alphaIon << "   beta ion = " << betaIon << endl
         << endl;
}

//------------------------------------------------------------------------------

void Calculus::rapidMKM_Kase_alphaIon_betaIon(double &alphaIon,
                                              double &betaIon)
{
    double dummy1, dummy2, dummy3;
    double alpha_X, beta_X;
    cellLine.getParameters(alpha_X, beta_X, dummy1);
    
    double alpha_D, beta_D;
    double alpha_tot = 0.0, sqrt_beta_tot = 0.0;
    double meanLet = 0.0;
    double r_domain = cellLine.getDomainRadius();
    double r_nucleus = cellLine.getNucleusRadius();
    
    CellLine fakeD("domain", r_domain);
    fakeD.addParametrization_LQ_noDt(alpha_X, beta_X);
    fakeD.setParametrization("LQ_noDt");
    Nucleus_Integral domain(fakeD);
    
    for(int k = 0; k < tracks.size(); k++)
    {
        Track *track = tracks[k].clone();
        double r_max = track->getRadius();
        double let = tracks[k].getLet();
        double weight = tracks[k].getWeight();
        
        //Compute gamma
        double r1;
        double r2 = 1e-4;
        double log10_r2 = log10(r2);
        track->setPosition(0.0, r2*1e-3);
        domain.distributeDose(*track);
        double d1, d2;
        domain.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
        domain.cleanNucleus();
        
        double z_1 = 0.;
        double z_1D = 0.;
        do
        {
            r1 = r2;
            d1 = d2;
            log10_r2 += 1e-3; // orig 1e-2
            r2 = pow(10.0, log10_r2);
            
            track->setPosition(0.0, r2*1e-3);
            domain.distributeDose(*track);
            domain.getDoseAndSurvival(d2, dummy1, dummy2, dummy3);
            domain.cleanNucleus();
            z_1 += (d1*r1 + d2*r2) * (r2 - r1);
            z_1D += (d1*d1*r1 + d2*d2*r2) * (r2 - r1);
        } while( r2 < r_domain + r_max );
        double gamma = z_1D / z_1;
        
        //Compute gamma_nucleus
        const double CONV = 160.2177;
        double gamma_nucleus = CONV*track->getLet()/(tracks.getDensity()*M_PI*r_nucleus*r_nucleus);
        
        double alpha_P = alpha_X + beta_X * gamma;
        double beta_P = beta_X;
        alpha_D = (1 - exp( - alpha_P * gamma_nucleus )) / gamma_nucleus;
        beta_D = beta_P;
        
        alpha_tot += weight * let * alpha_D;
        sqrt_beta_tot += weight * let * sqrt(beta_D);
        meanLet += weight * let;
        
    }
    
    alphaIon = alpha_tot / meanLet;
    betaIon = pow(sqrt_beta_tot/meanLet, 2);
    
    cout << setprecision(6)
         << "alpha ion = " << alphaIon << "   beta ion = " << betaIon << endl
         << endl;
}

//------------------------------------------------------------------------------

void Calculus::rapidRusso_alphaIon_betaIon(double &alphaIon,
                                           double &betaIon)
{
    // Nella tesi di Germano è quello che viene indicato con "A refinement of the Rapid GSI approach"
    double alpha_X, beta_X, D_t;
    cellLine.getParameters(alpha_X, beta_X, D_t);
    double s = alpha_X + 2 * beta_X * D_t;  //Gy^(-1)
    
    const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
    
    Track *track;
    double dummy1, dummy2;
    double doseSingleImpact;
    double survivalSingleImpact;
    double alpha_D, beta_D;
    double alpha_tot = 0.0, sqrt_beta_tot = 0.0;
    double meanLet = 0.0;
    double alpha_z = 0.0, beta_z = 0.0;
    
    for(int k = 0; k < tracks.size(); k++)
    {
        track = tracks[k].clone();
        track->setPosition(0.0, 0.0);
        
        nucleus.cleanNucleus();
        nucleus.distributeDose(*track);
        nucleus.getDoseAndSurvival(doseSingleImpact, dummy1,
                                   survivalSingleImpact, dummy2);
        nucleus.cleanNucleus();
        
        delete track;
        
        double area_nucleus = M_PI * pow(nucleus.getRadius(), 2);
        double let = tracks[k].getLet();
        double density = tracks.getDensity();
        double weight = tracks[k].getWeight();
        
        alpha_z = alpha_X * (1 - (density * area_nucleus * doseSingleImpact) / (CONV * let)) -
                    log(survivalSingleImpact) * (density * area_nucleus) / (CONV * let);
        alpha_D = alpha_X * (1 - (density * area_nucleus * doseSingleImpact) / (CONV * let)) +
                    (1 - survivalSingleImpact) * (density * area_nucleus) / (CONV * let);
        
        if( D_t == 0 )
            beta_z = 0;
        else
        {
            beta_z = (s - alpha_z) / (2 * D_t);
            if( beta_z < 0 )
                beta_z = 0;
        }
        beta_D = pow(alpha_D/alpha_z, 2) * beta_z;
        
        alpha_tot += weight * let * alpha_D;
        sqrt_beta_tot += weight * let * sqrt(beta_D);
        meanLet += weight * let;
    }
    
    alphaIon = alpha_tot / meanLet;
    betaIon = pow(sqrt_beta_tot/meanLet, 2);
    
    cout << setprecision(6)
         << "alpha ion = " << alphaIon << "   beta ion = " << betaIon << endl
         << endl;
}

//------------------------------------------------------------------------------

void Calculus::rapidScholz_alphaIon_betaIon(double &alphaIon,
                                            double &betaIon)
{
    // Nella tesi di Germano è quello che viene indicato con "GSI approximate implementation"
    double alpha_X, beta_X, D_t;
    cellLine.getParameters(alpha_X, beta_X, D_t);
    double s = alpha_X + 2 * beta_X * D_t;  //Gy^(-1)
    
    Track *track;
    double doseSingleImpact, survivalSingleImpact;
    double tmp1, tmp2, tmp3;
    double alpha_D, beta_D;
    double alpha_tot = 0.0, sqrt_beta_tot = 0.0;
    double meanLet = 0.0;
    double alpha_z = 0.0, beta_z = 0.0;
    
    for(int k = 0; k < tracks.size(); k++)
    {
        track = tracks[k].clone();
        track->setPosition(0.0, 0.0);
        nucleus.distributeDose(*track);
        
        nucleus.getDoseAndSurvival(tmp1, tmp2, survivalSingleImpact, tmp3);
        nucleus.cleanNucleus();
        
        const double CONV = 160.2177;  //(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
        
        doseSingleImpact = CONV * track->getLet() / (tracks.getDensity() * M_PI * pow(nucleus.getRadius(), 2));
        
        alpha_z = -log(survivalSingleImpact) / doseSingleImpact;
        
        if( D_t == 0 )
            beta_z = 0;
        else
        {
            beta_z = (s - alpha_z) / (2 * D_t);
            if( beta_z < 0 )
                beta_z = 0;
        }
        
        alpha_D = (1 - survivalSingleImpact) / doseSingleImpact;
        beta_D = pow(alpha_D/alpha_z, 2) * beta_z;
        
        alpha_tot += track->getWeight() * track->getLet() * alpha_D;
        sqrt_beta_tot += track->getWeight() * track->getLet() * sqrt(beta_D);
        meanLet += track->getWeight() * track->getLet();
        
        delete track;
    }
    
    alphaIon = alpha_tot / meanLet;
    betaIon = pow(sqrt_beta_tot/meanLet, 2);
    
    cout << setprecision(6)
         << "alpha ion = " << alphaIon << "   beta ion = " << betaIon << endl
         << endl;
}

//------------------------------------------------------------------------------
/*
void Calculus::prova_write_means(ofstream* file,
                                 const string calculusType,
                                 const int nFraction, // va passato
                                 const double timeSpacing, // va passato
                                 const double fracDeliveryTime,
                                 const vector<double> dosesImposed, // Va passata
                                 const vector<double> meanDose,
                                 const vector<double> meanDoseUncertainty,
                                 const vector<double> meanSurvival,
                                 const vector<double> meanSurvivalUncertainty)
{
    bool isMono = tracks.isMonoenergetic();
    MeansFile << model << ", " << calculusType << ", " << cellLine.getCellType() << ", ";
    double p0=0.0;
    double p1=0.0;
    double p2=0.0;
    double p3=0.0;
    double p4=0.0;
    double tmp1;
    long int tmp2;
    if (model=="LEMI") {
        cellLine.getParameters(p0,p1,p3);
        p3=cellLine.getNucleusRadius();
    }
    else if (model=="LEMII") {
        cellLine.getParameters_LQ2(p0,p1,p3,tmp1,tmp1,tmp1,tmp2);
        p3=cellLine.getNucleusRadius();
    }
    else if (model=="LEMIII") {
        cellLine.getParameters_LQ3(p0,p1,p3,tmp1,tmp1,tmp1,tmp2);
        p3=cellLine.getNucleusRadius();
    }
    else if (model=="MKM") {
        cellLine.getParameters_LQ_noDt(p0,p1);
        p2=cellLine.getNucleusRadius();
        p3=cellLine.getDomainRadius();
    }
    else if (model=="tMKM") {
        cellLine.getParameters_LQ_noDt(p0,p1,p4);
        p2=cellLine.getNucleusRadius();
        p3=cellLine.getDomainRadius();
    }
    MeansFile << p0 << ", " << p1 << ", " << p2 << ", " << p3 << ", ";
    if (model=="tMKM")
        MeansFile << p4 << ", ";
    if (isMono)
        MeansFile << "Monoenergetic, " << tracks[0].getParticleType();
    else
        MeansFile << "Spectrum, " << tracks.getSpectrumFile();
    MeansFile << ", " << tracks.getMeanEnergy() << ", " << tracks.getSigmaMeanEnergy() << ", " << tracks.getMeanLet() << ", "
              << tracks.getSigmaMeanLet() << ", " << tracks.getDoseAveragedLet() << ", " << tracks.getSigmaDoseAveragedLet() << ", "
              << nFraction << ", " << timeSpacing << ", " << fracDeliveryTime << ", " << dosesImposed[i] << ", ";
    MeansFile << scientific << meanDose[i] << ", " << meanSurvival[i] << ", "
              << meanDoseUncertainty[i] << ", " << meanSurvivalUncertainty[i] << endl;
    MeansFile.unsetf(std::ios::scientific);
    
}
*/
//------------------------------------------------------------------------------

void Calculus::slow_alphaIon_betaIon(const string trackMode,
                                     const vector<double> parameters,
                                     const vector<double> dosesImposed,
                                     const double precision,
                                     double& alphaIon,
                                     double& alphaIonUncertainty,
                                     double& betaIon,
                                     double& betaIonUncertainty,
                                     const int nFraction,
                                     const double timeSpacing,
                                     const double fracDeliveryTime,
                                     const bool saveAlphaBeta,
                                     const bool saveMeans,
                                     const bool saveCell,
                                     const string title_means)
{
    vector<double> meanDose(dosesImposed.size()), meanDoseUncertainty(dosesImposed.size());
    vector<double> meanSurvival(dosesImposed.size()), meanSurvivalUncertainty(dosesImposed.size());

    bool isMono = tracks.isMonoenergetic();
    
    ofstream MeansFile;
    if (saveMeans){
        MeansFile.open(title_means.c_str(), ios::app);
        if(!MeansFile)
        {
            cerr << "File " << title_means << " could not be opened." << endl;
            exit(1);
        }
    }
    
    for (size_t i=0; i<dosesImposed.size(); i++) {
        
        slow_meanDose_meanSurvival(trackMode, dosesImposed[i], precision,
                                   meanDose[i], meanDoseUncertainty[i], meanSurvival[i], meanSurvivalUncertainty[i],
                                   nFraction, timeSpacing, fracDeliveryTime, saveCell);
        
        // Save Mean Data
        if (saveMeans) {
            double meanEn = tracks.getMeanEnergy();
            double meanEnSig = tracks.getSigmaMeanEnergy();
            double meanL = tracks.getMeanLet();
            double meanLSig = tracks.getSigmaMeanLet();
            double meanDoseAvL = tracks.getDoseAveragedLet();
            double meanDoseAvLSig = tracks.getSigmaDoseAveragedLet();
            
            MeansFile << model << ",MonteCarlo," << cellLine.getCellType() << ",";
            for (size_t a=0; a<parameters.size(); a++)
                MeansFile << parameters[a] << ",";
            if (isMono)
                MeansFile << "Monoenergetic," << (tracks.getSpectrumFile()=="no_File" ? tracks[0].getParticleType() : tracks.getSpectrumFile());
            else
                MeansFile << "Spectrum," << tracks.getSpectrumFile();
            MeansFile << "," << meanEn << "," << (meanEnSig<1e-12 ? 0 : meanEnSig) << ","
                      << meanL << "," << (meanLSig<1e-12 ? 0 : meanLSig) << ","
                      << meanDoseAvL << "," << (meanDoseAvLSig<1e-12 ? 0 : meanDoseAvLSig) << ","
                      << nFraction << "," << timeSpacing << "," << fracDeliveryTime << "," << dosesImposed[i] << ",";
            MeansFile << scientific << meanDose[i] << "," << meanSurvival[i] << ","
                      << meanDoseUncertainty[i] << "," << meanSurvivalUncertainty[i] << endl;
            MeansFile.unsetf(std::ios::scientific);
        }
    }
    if (saveMeans)
        MeansFile.close();
    
    // fit data
    if (saveAlphaBeta) {
        double chiSquared;
        double incompleteGammaQ;
    
        fit_LQ(dosesImposed,
               meanSurvival, meanSurvivalUncertainty,
               alphaIon, alphaIonUncertainty,
               betaIon, betaIonUncertainty,
               chiSquared, incompleteGammaQ);
    
        cout << endl
             << "alpha ion = " << setprecision(6) << alphaIon
             << " (" << setprecision(2) << alphaIonUncertainty/abs(alphaIon) * 100.0 << "%)"
             << "      beta ion = " << setprecision(6) << betaIon
             << " (" << setprecision(2) << betaIonUncertainty/abs(betaIon) * 100.0 << "%)"
             << "      chi squared = " << chiSquared
             << "      goodness of fit = " << incompleteGammaQ << endl
             << endl;
    }
}

//------------------------------------------------------------------------------

void Calculus::slow_alphaIon_betaIon_with_Domains(const string trackMode,
                                                  const double minDose,
                                                  const double maxDose,
                                                  const int numberOfDoses,
                                                  const double precision,
                                                  double &alphaIon,
                                                  double &alphaIonUncertainty,
                                                  double &betaIon,
                                                  double &betaIonUncertainty)
{
    if( minDose >= maxDose )
    {
        cerr << "Mininum dose and maximum dose are not set correctly in function slow_alphaIon_betaIon" << endl;
        exit(1);
    }
    
    vector<double> doseImposed(numberOfDoses);
    vector<double> meanDose(numberOfDoses), meanDoseUncertainty(numberOfDoses);
    vector<double> meanSurvival(numberOfDoses), meanSurvivalUncertainty(numberOfDoses);
    
    for(int i = 0; i < numberOfDoses; i++ )
    {
        // logarithmically spaced doses
        //    doseImposed = minDose * pow(maxDose/minDose, static_cast<double>(i)/(numberOfDoses-1));
        doseImposed[i] = minDose + i/(numberOfDoses-1.0) * (maxDose - minDose);
        
        slow_meanDose_meanSurvival_with_Domains(trackMode,
                                                doseImposed[i], precision,
                                                meanDose[i], meanDoseUncertainty[i],
                                                meanSurvival[i], meanSurvivalUncertainty[i]);
    }
    
    double chiSquared;
    double incompleteGammaQ;
    
    fit_LQ(doseImposed,
           meanSurvival, meanSurvivalUncertainty,
           alphaIon, alphaIonUncertainty,
           betaIon, betaIonUncertainty,
           chiSquared, incompleteGammaQ);
    
    cout << endl
         << "alpha ion = " << setprecision(6) << alphaIon
         << " (" << setprecision(2) << alphaIonUncertainty/abs(alphaIon) * 100.0 << "%)"
         << "      beta ion = " << setprecision(6) << betaIon
         << " (" << setprecision(2) << betaIonUncertainty/abs(betaIon) * 100.0 << "%)"
         << "      chi squared = " << chiSquared
         << "      goodness of fit = " << incompleteGammaQ << endl
         << endl;
}

//------------------------------------------------------------------------------

void Calculus::slow_meanDose_meanSurvival(const string trackMode,
                                          const double doseImposed,
                                          const double precision,
                                          double &meanDose,
                                          double &meanDoseUncertainty,
                                          double &meanSurvival,
                                          double &meanSurvivalUncertainty,
                                          const int nFraction,
                                          const double timeSpacing,
                                          const double fracDeliveryTime,
                                          const bool saveCell)
{
    bool num=false, dev=false, abort=false;
    
    int numberOfIterations = 10;
    double relativeStdDeviation = 1e+6;
    double dump;
    if( 0.0 < precision && precision < 1.0 )
        relativeStdDeviation = precision;
    else if( 1.0 <= precision && modf(precision, &dump) == 0.0){
        numberOfIterations = static_cast<int>(precision);
        relativeStdDeviation=0.0;
        dev=true;
    }
    else
    {
        cerr << "Precision has not been set correctly." << endl;
        exit(1);
    }
    
    double dose=0.0, doseUncertainty=0.0, survival=0.0, survivalUncertainty=0.0;
    double sumSurvivalSquared, varSurvival;
    double sumDoseSquared, varDose;
    
    cout << setprecision(6) << fixed << left;
    
    ostringstream survivalFile_name;
    ofstream survivalFile;
    if(saveCell) {
        survivalFile_name << savePrefix << setprecision(2) << fixed << doseImposed << "(Gy).dat";
        survivalFile.open(survivalFile_name.str().c_str(), ios::app);
        if(!survivalFile)
        {
            cerr << "File " << survivalFile_name.str() << " could not be opened." << endl;
            exit(1);
        }
        survivalFile << setprecision(6) << scientific << left;
    }
    
    cout << endl
         << "Beginning to compute mean survival of " << nucleus.getCellType() << " cells" << endl
         << "following the irradiation with a pencil beam" << endl
         << "that deposite an amount of dose equal to " << setprecision(2) << doseImposed << " Gy" << endl
         << "in " << nFraction << " fractions." << endl
         << "Time spacing between fractions: " << timeSpacing << " hours." << endl
         << "Single fraction delivery time: " << fracDeliveryTime << " hours." << endl
         << endl;
    
    meanDose = 0.0;
    sumDoseSquared = 0.0;
    meanSurvival = 0.0;
    sumSurvivalSquared = 0.0;
    
    double dosePerFraction = doseImposed/(double)nFraction;
    double start_delivery_t=0.0;
    int it=0, thread=0;
    double tmp1 = 0.0, tmp2 = 0.0;
    
    double shIter=0.0;   // shared iterator
    
    omp_set_num_threads(nThreads);
    vector<Nucleus*> nucleus_vec(nThreads);
    for (int it=0; it<nThreads; ++it)
        nucleus_vec[it]=nucleus.clone(cellLine);
    
    time_t start = time(NULL);
    
#pragma omp parallel private(dose, doseUncertainty, survival, survivalUncertainty, it, start_delivery_t, thread)
    {
        while (!abort) {
            thread = omp_get_thread_num();
            if (model!="tMKM")
            {
                if (trackMode == "random")
                    random_dose_survival_p(doseImposed,
                                           dose, doseUncertainty,
                                           survival, survivalUncertainty,
                                           *nucleus_vec[thread]);
                else if (trackMode == "histogram")
                    histogram_dose_survival_p(doseImposed,
                                              dose, doseUncertainty,
                                              survival, survivalUncertainty,
                                              *nucleus_vec[thread]);
                else {
                    cerr << "The selected track mode does not exist" << endl;
                    exit(1);
                }
            }
            else {
                for (it=0; it<nFraction; ++it) {
                    start_delivery_t = timeSpacing * it;
                    histogram_dose_survival_t(*nucleus_vec[thread], dosePerFraction, start_delivery_t, fracDeliveryTime);
                }
                nucleus_vec[thread]->getDoseAndSurvival(dose,tmp1,survival,tmp2);
            }
#pragma omp critical
            {
                if(!isnan(dose) && !isnan(survival) && !abort){
                    shIter+=1.0;
                    meanDose = (meanDose * (shIter-1.0) + dose) / shIter;
                    sumDoseSquared += dose*dose;
                    varDose = sumDoseSquared / shIter - meanDose*meanDose;  //biased estimation
                    meanDoseUncertainty = sqrt(varDose / shIter);
                    
                    meanSurvival = (meanSurvival * (shIter-1.0) + survival) / shIter;
                    sumSurvivalSquared += survival*survival;
                    varSurvival = sumSurvivalSquared / shIter - meanSurvival*meanSurvival;  //biased estimation
                    meanSurvivalUncertainty = sqrt(varSurvival / shIter);
                    
                    if(difftime(time(NULL), start) > 5) {
                        cout << "iter. n. " << (int)shIter << ")"
                             << "   Mean dose = " << setprecision(8) << meanDose
                             << "   Mean survival = " << setprecision(8) << meanSurvival
                             << "   Precision = " << setprecision(2) << meanSurvivalUncertainty/meanSurvival*100 << "%" << endl;
                        start = time(NULL);
                    }
                    
                    if (saveCell) survivalFile << setw(15) << dose << setw(15) << survival << endl;
                    if (shIter >= numberOfIterations)
                        num=true;
                    if (num && meanSurvivalUncertainty/meanSurvival <= relativeStdDeviation && abs(meanDose-doseImposed)/doseImposed<precision)
                        dev=true;
                    if(num && dev) abort=true;
                }//isnan
            }//critical
            if (model=="tMKM")
                nucleus_vec[thread]->cleanNucleus();
        }//while
    }//parallel
    
    if (saveCell) survivalFile.close();
    
    cout << endl
         << "Mean dose: " << setprecision(6) << meanDose
         << " (" << setprecision(2) << meanDoseUncertainty/meanDose * 100.0 << "%)"
         << "          mean survival: " << setprecision(6) << meanSurvival
         << " (" << setprecision(2) << meanSurvivalUncertainty/meanSurvival * 100.0 << "%)" << endl
         << endl;
    if(saveCell) cout << "Output saved in " << survivalFile_name.str() << endl;
    cout << endl
         << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
         << endl;
    
    for (int it=0; it<nThreads; ++it)
        delete nucleus_vec[it];
}

//------------------------------------------------------------------------------

void Calculus::slow_meanDose_meanSurvival_with_Domains(const string trackMode,
                                                       const double doseImposed,
                                                       const double precision,
                                                       double &meanDose,
                                                       double &meanDoseUncertainty,
                                                       double &meanSurvival,
                                                       double &meanSurvivalUncertainty)
{
    vector<double> doses(nucleus.getNumberOfDomains(), 0.0);
    vector<double> dosesUncertainty(nucleus.getNumberOfDomains(), 0.0);
    vector<double> lethals(nucleus.getNumberOfDomains(), 0.0);
    vector<double> lethalsUncertainty(nucleus.getNumberOfDomains(), 0.0);
    
    int numberOfIterations = 10;
    double relativeStdDeviation = 1e+6;
    double dump;
    if( 0.0 < precision && precision < 1.0 )
        relativeStdDeviation = precision;
    else if( 1.0 <= precision && modf(precision, &dump) == 0.0)
        numberOfIterations = static_cast<int>(precision);
    else
    {
        cerr << "Precision has not been set correctly." << endl;
        exit(1);
    }
    
    double dose, doseUncertainty, survival, survivalUncertainty;
    double sumSurvivalSquared, varSurvival;
    double sumDoseSquared, varDose;
    
    cout << setprecision(6) << fixed << left;
    
    ostringstream survivalFile;
    survivalFile << setprecision(2) << fixed;
    survivalFile << savePrefix << doseImposed << "Gy.dat";
    ofstream outSurvivalFile(survivalFile.str().c_str(), ios::app);
    if(!outSurvivalFile)
    {
        cerr << "File " << survivalFile.str() << " could not be opened." << endl;
        exit(1);
    }
    outSurvivalFile << setprecision(6) << scientific << left;
    
    cout << endl
         << "Beginning to compute mean survival (multi domain mode) of " << nucleus.getCellType() << " cells" << endl
         << "following the irradiation with a pencil beam" << endl
         << "that deposite an amount of dose equal to " << setprecision(2) << doseImposed << " Gy" << endl
         << endl;
    
    meanDose = 0.0;
    sumDoseSquared = 0.0;
    meanSurvival = 0.0;
    sumSurvivalSquared = 0.0;
    int j = 1;
    
    do
    {
        if( trackMode == "histogram" )
            histogram_dose_survival_with_domains(doseImposed, dose, doseUncertainty, survival, survivalUncertainty, doses, lethals, dosesUncertainty, lethalsUncertainty);
        else
        {
            cerr << "The selected track mode does not exist" << endl;
            exit(1);
        }
        
        meanDose = (meanDose * (j-1) + dose) / j;
        sumDoseSquared += dose*dose;
        varDose = sumDoseSquared / j - meanDose*meanDose;  //biased estimation
        meanDoseUncertainty = sqrt(varDose / j);
        
        meanSurvival = (meanSurvival * (j-1) + survival) / j;
        sumSurvivalSquared += survival*survival;
        varSurvival = sumSurvivalSquared / j - meanSurvival*meanSurvival;  //biased estimation
        meanSurvivalUncertainty = sqrt(varSurvival / j);
        
        double totalLethal = 0.0;
        for( size_t k = 0; k < (unsigned long)nucleus.getNumberOfDomains(); k++ )
            totalLethal += lethals[k];
        
        outSurvivalFile << setw(15) << dose << setw(15) << survival << setw(15) << totalLethal << endl;
        
        j++;
    } while( j <= numberOfIterations || meanSurvivalUncertainty/meanSurvival > relativeStdDeviation );
    
    outSurvivalFile.close();
    
    cout << endl
         << "Mean dose (domain mode): " << setprecision(6) << meanDose
         << " (" << setprecision(2) << meanDoseUncertainty/meanDose * 100.0 << "%)"
         << "          mean survival: " << setprecision(6) << meanSurvival
         << " (" << setprecision(2) << meanSurvivalUncertainty/meanSurvival * 100.0 << "%)" << endl
         << endl
         << "Output saved in " << survivalFile.str() << endl
         << endl
         << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl
         << endl;
}

//------------------------------------------------------------------------------

void Calculus::verbatim_dose_survival(double &dose,
                                      double &doseUncertainty,
                                      double &survival,
                                      double &survivalUncertainty,
                                      bool clean)
{
    cout << setprecision(6) << fixed << left;

    nucleus.distributeDose(tracks);

    nucleus.getDoseAndSurvival(dose, doseUncertainty, survival, survivalUncertainty);

    cout << "Tracks in nucleus: " << setw(4) << nucleus.getInNucleusCount()
         << "   tracks intersecting nucleus: " << setw(6) << nucleus.getIntersectionCount()
         << "     dose: " << setprecision(6) << dose
         << " (" << setprecision(2) << doseUncertainty/dose * 100.0 << "%)"
         << "     survival: " << setprecision(6) << survival
         << " (" << setprecision(2) << survivalUncertainty/survival * 100.0 << "%)" << endl;

    if (clean) nucleus.cleanNucleus();
}
