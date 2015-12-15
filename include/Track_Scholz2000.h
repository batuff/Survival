#ifndef TRACK_SCHOLZ_2000_H
#define TRACK_SCHOLZ_2000_H

#include "Track.h"

namespace Survival {
    
    class Particle;

    //! Inherited from the Track class, it implements the LEM I track model.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007
     
        The structure of the track is defined in a publication by Scholz and Kraft (\ref LEMI "1") but then (minimally) modified.
     
        According to this parametrization, the average local dose deposition pattern (amorphous track) \f$d(\rho)\f$ is divided in a \em core, with a radius \f$\rho_{min}=0.01\,um\f$, where the local dose is constantly equal to
        \f[
            d(\rho)=\lambda\frac{LET}{\rho_{min}^2}\qquad\rho<\rho_{min}
        \f]
        and an outer region which radius is indicated with \f$\rho_{max}\f$ and coincides with the radius of the track where the local dose is defined by
        \f[
            d(\rho)=\lambda\frac{LET}{\rho^2}\qquad\rho_{min}<\rho<\rho_{max}
        \f]
        while the local dose is equal to zero if the distance is greater than \f$\rho_{max}\f$. \f$\lambda\f$ is a parameter defined as:
        \f[
            \lambda=\frac{1}{\rho\pi(1+\ln(\rho_{max}^2/\rho_{min}^2))}
        \f]
        where \f$\rho\f$ in this formula represents the density of the medium.
     
        The radius \f$\rho_{max}\f$ of the track is evaluated as a function of the specific energy of the ion, as:
        \f[
            \rho_{max}=\gamma\,E^\delta
        \f]
        where \f$\gamma=0.062\,\frac{\mu m}{MeV^\delta}\f$ and \f$\delta=1.7\f$ are constant values.
     
        This class provides some methods to evaluate and get the local dose deposited by the ion along the track and the radial integral of the track.
     
        \anchor LEMI 1. M. Scholz and G. Kraft, "Track structure and the calculation of biological effects of heavy charged particles", \a Advances \a in \a Space \a Research \b 18, 5-14 (1996)
     
        \sa Track_Elsasser2007 and Track_Elsasser2008
     */
    class Track_Scholz2000 : public Track
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Converts a particle (object of class Particle), passed by reference, in a track according to the amorphous track model defined by Scholz \em et \em al.
         
            \param particle The particle generating the track in the medium, passed by reference.
            \param density The density of the medium expressed in \f$\frac{g}{cm^3}\f$. The default value is the density of water.
            \param t The time corresponding to the generation of the track in the target. The default value is 0. (See also the documentation of the data member #time)
         
            Also #r_max and #lambda are instantiated getting information from particle.
         */
        Track_Scholz2000(const Particle &particle,
                         const double density,   // kg/dm^3 or g/cm^3
                         double t = 0.0);
        
        //! Destructor.
        virtual ~Track_Scholz2000() {};
        
        //! Returns a pointer to a new Track_Scholz2000 object created as a copy of an existent one by means of the copy constructor.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         */
        virtual Track_Scholz2000* clone() const;
        
        //! Returns the distance from the center of the track, knowing the local dose deposited.
        /*!
            It's the inverse function of getLocalDose(). Since the maximum possible local dose deposited is: \f$d_{MAX}=\lambda\frac{LET}{\rho_{min}^2}\f$, if d is greater than \f$d_{MAX}\f$ it returns -1 (nonsense value). Else it returns:
            \f[
                \sqrt{\frac{\lambda\,LET}{d}}
            \f]
            where d is the local dose.
         
            \note If \f$d<d_{MIN}=\lambda\frac{LET}{\rho_{max}^2}\f$ it returns #r_max.
         
            \param localDose The local dose deposited expressed in Gy.
         
            \return The distance from the center of the track, expressed in um.
         
            \sa getLocalDose()
         */
        virtual double getDistance(const double localDose) const;  // Gy -> um
        
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
            According to the parametrization by Scholz \em et \em al., the function evaluates the local dose evaluating some possible cases:
                - If the distance \f$\rho\f$ is smaller than the core radius #R_MIN it returns \f$d(\rho)=\lambda\frac{LET}{R\_MIN^2}\f$
                - If the distance \f$\rho\f$ is greater than the radius of the track #r_max it returns 0
                - Else it returns \f$d(\rho)=\lambda\frac{LET}{\rho^2}\f$
         
            \param distance The distance from the track center expressed in um.
         
            \return The local dose (in Gy) deposited by the track at a certain distance from the track center (expressed in um).
         
            \sa #DELTA, #GAMMA and #lambda
         */
        virtual double getLocalDose(const double distance) const;  // distance from track center (um)
        
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
                                 double &returnY) const;
        
        //! Returns the dose (in Gy) deposited by the track evaluated as the radial integral of the track profile between \f$r_{begin}\f$ and \f$r_{end}\f$ (expressed in um).
        /*!
            The local dose varies along the track profile according to the structure defined by Scholz \em et \em al. (see Track_Scholz2000). This function evaluates the integral of the radial profile in \f$\left[r_{begin};\ r_{end}\right]\f$.
         
            \param r_begin Lower limit of integration, expressed in um.
            \param r_end Upper limit of integration, expressed in um.
         
            \return The dose (in Gy) deposited by the track evaluated as the radial integral of the track profile between \f$r_{begin}\f$ and \f$r_{end}\f$.
         
            \warning The execution of the program will be terminated if incorrect limits of integration are chosen, that is:
                - If \f$r_{begin}<0\f$
                - If \f$r_{end}<r_{begin}\f$
         
            \sa getLocalDose()
         */
        virtual double getRadialIntegral(const double r_begin,
                                         const double r_end) const;
        
        //! Returns the radius of the track expressed in um.
        /*!
            \return The radius of the track expressed in um.
         
            \sa #R_MIN
         */
        virtual double getRadius() const {return r_max;};
        
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
            The function divides the track radius (#r_max) into 300 logarithmically spaced distances from the track center, for each of these distances it evaluates the local dose deposited calling the function getLocalDose(); during the process it saves each distance and the corresponding dose in a new file. This allow to reconstruct the track profile.
         
            \return The name of the file created.
         
            \sa getLocalDose()
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
        virtual void setTime(double t) {time = t;};
        
    private:
        
        //! The type of particle generating the track (e.g. Chemical symbol for ions: H, He, Li, ...)
        /*!
            \sa Particle::type
         */
        std::string particleType;
        
        //! The specific energy of the particle generating the track, expressed in MeV/u.
        double particleEnergy;
        
        //! The radius of the track expressed in um.
        /*!
            According to the parametrization by Scholz and Kraft it's evaluated (and instantiated in the constructor) as:
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
            where \f$\rho\f$ represents the density of the medium while \f$\rho_{max}\f$ and \f$\rho_{min}\f$ represent #r_max and #R_MIN respectively.
         
            It's expressed in \f$\frac{Gy\,\mu m^3}{MeV}\f$.
         */
        double lambda;
        
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
        
        //! The time associated to a particular event, expressed in hours.
        /*!
            The \a time isn't really a variable associated to the track, time flows during the irradiation process. For this reason, the initial value of this member represent the instant in which the track is generated, but then it will change during the execution assuming a value representing the interaction time between the track and a specific Nucleus.
         
            \note Since this track structure is used in the LEM model, that doesn't take into account (yet) the time structure of the irradiation, this data member is actually useless.
         */
        double time;
        
        //! Useful constant value.
        /*!
            It's necessary to evaluate #r_max; it's equal to 1.7, according to the LEM I parametrization.
         */
        static const double DELTA;
        
        //! Useful constant value.
        /*!
            It's necessary to evaluate #r_max; it's equal to 0.062 \f$\frac{\mu m}{MeV^\delta}\f$, according to the LEM I parametrization.
         */
        static const double GAMMA;
        
        //! The core radius expressed in um.
        /*!
            It's taken equal to 0.01 um, according to the LEM I parametrization.
         */
        static const double R_MIN;
    };
    
}


#endif /* TRACK_SCHOLZ2000_H */
