#ifndef TRACK_ELSASSER2008_H
#define TRACK_ELSASSER2008_H

#include "Track.h"

namespace Survival {
    
    class Particle;

    //! Inherited from the Track class, it implements the LEM III track model.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2009
     
        To further increase the level of agreement with carbon ion data with respect to the LEM II formulation (see Track_Elsasser2007) a LEM III version (\ref LEMIII "1") was implemented where the track core radius #r_min has been made proportional to the particle velocity, according to the following relation:
        \f[
            \rho_{min}=\rho_c\cdot\beta
        \f]
        where \f$\rho_c=40\,nm\f$ (#R_C) represents the maximum core radius, corresponding to ions traveling at the speed of light, and \f$beta\f$ is the relativistic ion velocity.
     
        In this case the \a master \a curve approach (see Track_Elsasser2007) is not applicable, because the shape of track core deposition varies with the kinetic energy of the particle (due to the previous relation). Therefore the track profile is evaluated at the time of the object construction (object by object) and stored it as a member variable. This could bring issues with memory consumption if a large number of Track_Elsasser2008 objects needs to be instantiate at the same time.
     
        \anchor LEMIII 1. T. Elsässer, M. Krämer and M. Scholz, "Accuracy of the local effect model for the prediction of biologic effects of carbon ion beams \em in \em vitro and \em in \em vivo", \a International \a Journal \a of \a Radiation \a Oncology-\a Biology-\a Physics \b 71, 866-872 (2008).
     */
    class Track_Elsasser2008 : public Track
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Converts a particle (object of class Particle), passed by reference, in a track according to the LEM III amorphous track model.
            Some of the data members are instantiated on the basis of the informations stored in the Particle object. (For a more detailed description of the instantiation of each member respectively look at its specific documentation).
         
            It constructs also the track radial profile via a process similar to the one defined in the Track_Elssasser2007 implementation but without the general purpose of the \em master \em curve.
         
            \param particle The particle generating the track in the medium, passed by reference.
            \param density The density of the medium expressed in \f$\frac{g}{cm^3}\f$. The default value is the density of water.
            \param doseCutoff The minimum possible dose deposited evaluable, expressed in Gy (see #doseCutoff).
            \param lengthDoseCurve The length of the curve representing the track radial profile, expressed in um.
            \param integrationStepFactor Dimensionless integration factor (multiplied by #r_min gives the #integrationStep)
            \param besselLimit The limit between the calculation of Bessel function by means of series development and the calculation with the asymptotic exponential approximation (see #besselLimit).
            \param numberOfSigma Half width of non-zero window in units of sigma (see #numberOfSigma).
            \param t The time corresponding to the generation of the track in the target. The default value is 0. (See #time).
         */
        Track_Elsasser2008(const Particle& particle,
                           const double density,
                           const double doseCutoff = 1e-8,
                           const int lengthDoseCurve = 300.0,
                           const double integrationStepFactor = 1e-2,
                           const double besselLimit = 400.0,
                           const double numberOfSigma = 20.0,
                           double t = 0.0);
        
        //! Copy constructor. Instantiates a new Track_Elsasser2008 object copying an existent one, including the #logRhoCurve and #logDoseCurve values.
        Track_Elsasser2008(const Track_Elsasser2008 &track);
        
        //! Destructor.
        /*!
            Delete the object, deallocating also the memory occupied by #logRhoCurve and #logDoseCurve, and decrements the counter of Track_Elsasser2008 objects (#numberOfElsasser2008Tracks).
         */
        virtual ~Track_Elsasser2008();
        
        //! Returns a pointer to a new Track_Elsasser2008 object created as a copy of an existent one by means of the copy constructor.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         */
        virtual Track_Elsasser2008* clone() const;
        
        //! Returns the distance from the center of the track, knowing the local dose deposited.
        /*!
            It's the inverse function of getLocalDose(). The function runs in a loop over the precalculated values of the track profile stored in #logRhoCurve and #logDoseCurve data members until it finds a value smaller than the required dose deposited, then it interpolates the nearest neighbors to get the correct value.
         
            \note If the dose deposited is smaller than the #doseCutoff it returns #r_eff.
         
            \param localDose The local dose deposited, expressed in Gy.
         
            \return The distance from the center of the track, expressed in um.
         
            \sa getLocalDose()
         */
        virtual double getDistance(const double localDose) const;
        
        //! Returns the kinetic energy of the particle generating the track expressed in MeV.
        /*!
            \return The kinetic energy of the particle generating the track expressed in MeV.
         
            \sa #e_c
         */
        virtual double getKineticEnergy() const {return e_c;};
        
        //! Returns the LET in water of the particle generating the track expressed in MeV/um.
        /*!
            \return The LET in water of the particle generating the track expressed in MeV/um (according to the Bethe-Bloch formula).
         
            \sa #let
         */
        virtual double getLet() const {return let;};
        
        //! Returns the local dose (in Gy) deposited by the track at a certain distance from the track center (expressed in um).
        /*!
            The function evaluates some possible cases:
                - If the distance is smaller than the minimum radius stored in the #logRhoCurve it returns the dose at the minimum radius evaluated
                - If the distance is greater than the effective radius of the track (#r_eff) it returns 0
                - Else it returns the precalculated local dose at the required distance obtained by an interpolation of the nearest neighbors.
         
            \param distance The distance from the track center expressed in um.
         
            \return The local dose (in Gy) deposited by the track at a certain distance from the track center (expressed in um).
         */
        virtual double getLocalDose(const double distance) const;
        
        //! Function created to calculate the mean time required to the evaluation the getLocalDose() method.
        /*!
            It cyclically calls getLocalDose() 1000000 times, timing the total elapsed time and dividing it by 1000000.
         
            \return The mean time needed to a complete evaluation of the getLocalDose() method, expressed in s.
         */
        double getLocalDoseMeanTime();
        
        //! Returns the specific energy of the particle generating the track, expressed in MeV/u.
        /*!
            \return The specific energy of the particle generating the track, expressed in MeV/u.
         
            \sa #particleEnergy
         */
        virtual double getParticleEnergy() const {return particleEnergy;};
        
        //! Returns the type of particle generating the track (e.g. Chemical symbol for ions: H, He, Li, ...).
        /*!
            \return The type of particle generating the track (e.g. Chemical symbol for ions: H, He, Li, ...).
         
            \sa Particle::type
        */
        virtual std::string getParticleType() const {return particleType;};
        
        //! Returns the track position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the track referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the track, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c x coordinate of the track, expressed in mm, passed by reference to be overwritten.
         
            \sa setPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const
        {
            returnX = x_track;
            returnY = y_track;
        };
        
        //! Evaluates the radial integral of the track profile in \f$\left[r_{min}, r_{max}\right]\f$.
        /*!
            \warning Not yet implemented.
         */
        virtual double getRadialIntegral(const double r_min,
                                         const double r_max) const;  // extremes of integration (um)
        
        //! Returns the effective radius of the track, expressed in um.
        /*!
            \return The effective radius of the track (#r_eff) expressed in um.
         */
        virtual double getRadius() const {return r_eff;};
        
        //! Evaluates the relative precision of the calculated LET with respect of the original particle LET.
        /*!
            The relative precision is evaluated by the difference between the calculated and the "original" LETs divided by "original" LET.
         
            \return The relative precision of the calculated LET with respect of the original particle LET.
         
            \sa getTrackLet()
         */
        double getRelativePrecision() const;
        
        //! Returns the time associated to a particular event expressed in hours.
        /*!
            \return The time associated to a particular event expressed in hours.
         
            \sa #time and setTime()
         */
        virtual double getTime() const {return time;};
        
        //! Returns the weight in the beam of the particle generating the track. Useful in the case of "mixed fields".
        /*!
            \return The weight in the beam of the particle generating the track. Useful in the case of "mixed fields".
         
            \sa Particle::weight
         */
        virtual double getWeight() const {return weight;};
        
        //! Saves the local dose deposited by the track at different distances from its center to reconstruct the track profile.
        /*!
            It execute a \c for loop saving the values stored in the #logRhoCurve and #logDoseCurve data members; that is the local dose deposited and corresponding radii along the whole radial profile.
         
            \return The name of the file created.
         */
        virtual std::string saveTrack() const;
        
        //! Sets the track position (\c x and \c y coordinates) referred to the beam axis.
        /*!
            \param x The \c x coordinate of the track to be set, referred to the beam axis and expressed in mm.
            \param y The \c x coordinate of the track to be set, referred to the beam axis and expressed in mm.
         
            \sa #x_track, #y_track and getPosition()
         */
        virtual void setPosition(const double x,
                                 const double y);
        
        //! Sets the time associated to a particular event.
        /*!
            \param t The time to be set expressed in hours.
         
            \sa #time
         */
        virtual void setTime(double t) {time=t;};
        
        
    private:
        
        //! Evaluates and returns the LET of the track.
        /*!
            Since it was observed a minimal discrepancy from the imposed particle LET, this function calculates the real observed LET of the particle integrating the radial profile.
         
            \return The calculated LET of the track, expressed in \f$\frac{MeV}{um}\f$
         */
        double getTrackLet() const;
        
        //! Evaluates the argument of the normalized integral in the construction of the track profile.
        /*!
            The calculation is divided in two cases:
                - If the argument of the Bessel's function \f$\rho=\frac{r\,r'}{\sigma^2}\f$ is smaller than the fixed #besselLimit then the evaluation is based on a series development of the Bessel function
                - If the argument of the Bessel's function \f$\rho=\frac{r\,r'}{\sigma^2}\f$ is grater than the fixed #besselLimit then it's used an asymptotic exponential approximation.
         
            \param r The radial coordinate of the track profile, expressed in um.
            \param r1 The radial coordinate of the gaussian function, expressed in um.
         
            \return The argument of the integral defining the convolution between the standard radial profile and the gaussian function.
         
            \sa normalizedPunctualDose()
         */
        double normalizedDoseIntegralArgument(const double r,
                                              const double r1) const; // normalized argument of dose integral
        
        //! Evaluates the local dose along the radial profile of the track.
        /*!
            The integration process is based on the trapezoidal rule by Newton-Cotes and, step by step, the argument of the integral is evaluated by means of the normalizedDoseIntegralArgument() method.
         
            \param distance The distance from the track center, expressed in um.
         
            \return The local dose at a fixed distance from the track center.
         */
        double normalizedPunctualDose(const double distance) const;
        
        //! The specific energy of the particle generating the track, expressed in MeV/u.
        double particleEnergy;
        
        //! The type of particle generating the track (e.g. Chemical symbol for ions: H, He, Li, ...)
        /*!
            \sa Particle::type
         */
        std::string particleType;
        
        //! The track position (\c x coordinate) referred to the beam axis, expressed in mm.
        /*!
            \sa Particle::x
         */
        double x_track;
        
        //! The track position (\c y coordinate) referred to the beam axis, expressed in mm.
        /*!
            \sa Particle::y
         */
        double y_track;
        
        //! The weight in the beam of the particle generating the track. Useful in the case of "mixed fields".
        /*!
            \sa Particle::weight
         */
        double weight;
        
        //! The density of the medium expressed in \f$\frac{g}{cm^3}\f$.
        double density;
        
        //! The core radius expressed in um.
        /*!
            According to the LEM III formulation, it's proportional to the relativistic ion velocity:
            \f[
                \rho_{min}=\rho_c\cdot\beta
            \f]
         
            \sa R_C
         */
        double r_min;
        
        //! The radius of the track expressed in um.
        /*!
            According to the LEM III parametrization it's evaluated (and instantiated in the constructor) as:
            \f[
                r\_max=\gamma\,E^\delta
            \f]
            where \f$\gamma\f$ is #GAMMA, \f$\delta\f$ is #DELTA and \em E represents the specific energy of the ion.
         */
        double r_max;
        
        //! The LET in water of the particle generating the track expressed in MeV/um (according to the Bethe-Bloch formula).
        /*!
            \sa Particle::let
         */
        double let;
        
        //! The kinetic energy of the particle generating the track expressed in MeV.
        /*!
            \sa Particle::e_c
         */
        double e_c;
        
        //! A constant value required to evaluate #r_max.
        /*!
            It's defined as:
            \f[
                \lambda=\frac{1}{\rho\pi(1+\ln(\rho_{max}^2/\rho_{min}^2))}
            \f]
            where \f$\rho\f$ represents the density of the medium while \f$\rho_{max}\f$ and \f$\rho_{min}\f$ represent #r_max and #r_min respectively.
         
            It's expressed in \f$\frac{Gy\,\mu m^3}{MeV}\f$.
         */
        double lambda;
        
        //! The effective radius of the track (expressed in um) corresponding to the point where the local dose results equal to #doseCutoff.
        double r_eff;
        
        //! Minimum possible dose deposited evaluable, expressed in Gy.
        double doseCutoff;
        
        //! The integration step for the construction of the track profile expressed in um. It's created proportional to #r_min.
        double integrationStep;
        
        //! Limit between the calculation of Bessel function by means of series development and the calculation with the asymptotic exponential approximation
        double besselLimit;
        
        //! Half width of non-zero window in units of sigma.
        double numberOfSigma;
        
        //! The length of the dose profile expressed in um.
        int lengthDoseCurve;
        
        //! A pointer to the (logarithmic) radii corresponding to the calculated local dose of the track profile.
        double *logRhoCurve;
        
        //! A pointer to the values of the calculated (logarithmic) local dose of the track profile.
        double *logDoseCurve;
        
        //! The time associated to a particular event, expressed in hours.
        /*!
            The \a time isn't really a variable associated to the track, time flows during the irradiation process. For this reason, the initial value of this member represent the instant in which the track is generated, but then it will change during the execution assuming a value representing the interaction time between the track and a specific Nucleus.
         
            \note Since this track structure is used in the LEM model, that doesn't take into account (yet) the time structure of the irradiation, this data member is actually useless.
         */
        double time;
        
        //! Constants static variables and precalculated feature indices
        /*!
            It's the constant of conversion from \f$\frac{MeV\,dm^3}{Kg\mu m^3}\f$ to Gy. It's equal to \f$160.2177\,\frac{J\mu m^3}{MeV dm^3}\f$
         */
        static const double CONV;	//(J*um^3)/(MeV*dm^3) : Constant of conversion (MeV*dm^3)/(Kg*um^3) -> Gy
        
        //! Useful constant value.
        /*!
            It's necessary to evaluate #r_max; it's equal to 1.7, according to the LEM III parametrization.
         */
        static const double DELTA;	// exponent for calculation of maximum radius
        
        //! Useful constant value.
        /*!
            It's necessary to evaluate #r_max; it's equal to 0.062 \f$\frac{\mu m}{MeV^\delta}\f$, according to the LEM III parametrization.
         */
        static const double GAMMA;
        
        //! The maximum core radius, corresponding to ions traveling at the speed of light, expressed in um.
        static const double R_C;
        
        //! The radical diffusion length, expressed in um.
        /*!
            It's a constant value equal to 4 nm representing the spreading of the induced radical species taking place at longer time scales (a few microseconds).
         */
        static const double SIGMA;	// um : radical diffusion length
        
        //! Maximum length of the dose profile (i.e. maximum number of steps, equal to the length of the arrays #tmpLogRhoCurve and #tmpLogDoseCurve).
        static const int MAX_LENGTH_DOSE_CURVE = 1000;
        
        //! Array to temporary store the progressive radii corresponding to the profile of the curve.
        static double tmpLogRhoCurve[MAX_LENGTH_DOSE_CURVE];
        
        //! Array to temporary store the calculated values of local dose deposited, constituting the radial profile of the curve.
        static double tmpLogDoseCurve[MAX_LENGTH_DOSE_CURVE];
        
        //! The number of existing Track_Elsasser2008 objects.
        /*!
            It's incremented in the constructor and decremented in the destructor.
         
            \sa Track_Elsasser2008() and ~Track_Elsasser2008()
         */
        static int numberOfElsasser2008Tracks;
    };
    
}


#endif /* TRACK_ELSASSER2008_H */
