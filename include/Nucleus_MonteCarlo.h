#ifndef NUCLEUS_MONTECARLO_H
#define NUCLEUS_MONTECARLO_H

#include "Nucleus_Pixel.h"
#include "Tracks.h"

#include <gsl/gsl_rng.h>

namespace Survival {

    //! Inherited from the Nucleus_Pixel class, it performs the integration of the dose deposited by the track via the Monte Carlo \a importance \a sampling method.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2008
     
        The class was thought to be a test class, useful to perform a simulation making use of the Nucleus_Pixel class and verifying the correctness and the precision of the nucleus structure and the way it integrates the track radial profile.
     
        The approximation of the integrand is given by the pixel-wise constant dose profile generated with the Nucleus_Pixel class; to avoid having zero values in some points of the nucleus, an arbitrary local dose background of 1 Gy is added. This approximation of dose profile is then normalized by its integral and is interpreted as a probability distribution; the related cumulative distribution is used to extract a position in the nucleus where to compute the exact local dose value and the consequent local number of lethal events. The average of the obtained lethal events values divided by their probability of extraction gives an estimate of the integral. The precision of the estimate is of course dependent on the number of extractions; the algorithm cycles until the user defined precision is reached (in the probabilistic convergence sense).
     */
    class Nucleus_MonteCarlo : public Nucleus_Pixel
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Calls explicitly Nucleus_Pixel constructor. It differs only in the \c precision parameter, used to set the ending condition of the Monte Carlo simulation. There are two way to set such a condition:
                - A fixed number of iterations, hence the \c precision has to be an integer value greater (or at least equal) to 1.
                - A constraint on the precision to reach in the simulation in the evaluation of the cell survival (precisely the relative error on the survival), hence the \c precision has to be a \c double in (0, 1).
         
            \warning The execution of the program will be terminated if the precision is not set correctly.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param precision Fix the ending condition of the Monte Carlo simulation.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c y coordinate of the center) referred to the beam axis, expressed in mm.
            \param pixelSide1 The side of the smallest (or the third) sub-grid of pixels (#pixelSide_1), expressed in um.
            \param scale1 The scale factor between the second and the third subgrid of pixels (#scale_1).
            \param radius1 The radius of the smallest circumference that defines the sampling of the track, expressed in um.
            \param scale2 The scale factor between the first and the second subgrid of pixels (#scale_2).
            \param radius2 The radius of the second circumference that defines the sampling of the track, expressed in um.
            \param scale3 The scale factor between the biggest grid of pixel and the first sub-grid (#scale_3).
            \param radius3 The radius of the biggest circumference that defines the sampling of the track, expressed in um.
         */
        Nucleus_MonteCarlo(const CellLine &cellLineRef,
                           const double precision = 3e-3,
                           const double xPosition = 0.0,
                           const double yPosition = 0.0,
                           const double pixelSide1 = 0.005,
                           const int scale1 = 2,
                           const double radius1 = 0.1,
                           const int scale2 = 10,
                           const double radius2 = 1.0,
                           const int scale3 = 10,
                           const double radius3 = 10.0);
        
        //! Destructor.
        ~Nucleus_MonteCarlo() {};
        
        //! Clean the object by calling the Tracks::eraseAll() and Nucleus_Pixel::cleanNucleus() methods.
        void cleanNucleus();
        
        //! It simply calls Nucleus_Pixel::distributeDose(const Track &track) to evaluate the dose deposited in the nucleus and appends the Track object to #distributedTracks.
        void distributeDose(const Track &track);
        
        //! Overload of distributeDose(const Track &track) to manage a Tracks object.
        /*!
            It simply calls the method Nucleus_Pixel::distributeDose(const Tracks &tracks) and appends the Tracks object to #distributedTracks.
         */
        void distributeDose(const Tracks &tracks);
        
        //! Perform a Monte Carlo simulation to integrate over the nucleus the dose deposited by a track distribution in a stochastic way.
        /*!
            The approximation of the integrand is given by the pixel-wise constant dose profile generated with the Nucleus_Pixel class; to avoid having zero values in some points of the nucleus, an arbitrary local dose background of 1 Gy is added. This approximation of dose profile is then normalized by its integral and is interpreted as a probability distribution; the related cumulative distribution is used to extract a position in the nucleus where to compute the exact local dose value and the consequent local number of lethal events. The average of the obtained lethal events values divided by their probability of extraction gives an estimate of the integral. The precision of the estimate is of course dependent on the number of extractions; the algorithm cycles until the user defined precision is reached (see Nucleus_MonteCarlo()).
         
            \param dose The total dose deposited in the nucleus by the radiation, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival observed related to the dose absorbed.
            \param survivalUncertainty The uncertainty associated to the cellular survival.
         */
        void getDoseAndSurvival(double &dose,
                                double &doseUncertainty,
                                double &survival,
                                double &survivalUncertainty);
        
    private:
        
        //! A Tracks object storing all Track objects interacting with the nucleus.
        Tracks distributedTracks;
        
        //! One of the two ending conditions of the Monte Carlo simulation. Fix the maximum number of iterations executable.
        long int numberOfIterations;
        
        //! One of the two ending conditions of the Monte Carlo simulation. Fix a constraint on the precision that is the maximum relative error on the cell survival evaluated.
        double relativeStdDeviation;
    };

}


#endif /* NUCLEUS_MONTECARLO_H */
