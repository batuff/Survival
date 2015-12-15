#ifndef TRACK_KIEFER_CHATTERJEE_H
#define TRACK_KIEFER_CHATTERJEE_H

#include "Track.h"

namespace Survival {
    
    class Particle;

    //! Inherited from the Track pure virtual class, it implements the Kiefer-Chatterje amorphous track structure, used in the MKM model.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Lorenzo Marengo
        \author Germano Russo
        \date 2011-2015
     
        The structure of the track is defined in the publication by Kase \a et \a al. (\ref KASE "1") in which they suggest to use the Kiefer-Chatterjee amorphous model.
        The track profile is characterized by an inner region, the \a core, and an outer region, the \a penumbra.
        The radius of the core and the penumbra respectively are evaluated by means of the following relations
        \f[
            R_{core}=\eta\cdot\beta_{Ion}
        \f]
        \f[
            R_{penumbra}=\gamma\left(\frac{E}{A}\right)^\delta
        \f]
        Where \f$\gamma\f$, \f$\delta\f$ and \f$\eta\f$ are constants defined in the publication reference (see below), \f$\beta_{Ion}\f$ represents the ratio between the speed of the ion and the speed of light, E is the energy of the ion and A his mass number.
        The values implemented for \f$\gamma\f$, \f$\delta\f$ and \f$\eta\f$ are:
            - \f$\gamma=0.0616\,\frac{\mu m}{MeV^\delta}\f$
            - \f$\delta=1.7\f$
            - \f$\eta=0.0116\,\mu m\f$
     
        This class provides some methods to evaluate and get the local dose deposited by the ion along the track and the radial integral of the track.
     
        \anchor KASE 1. Y. Kase, T. Kanai, N. Matsufuji, Y. Furusawa, T. Elsasser, and M. Scholz, "Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation", \a Physics \a in \a Medicine \a and \a Biology \b 53, 37-59 (2008).
     
        \sa Track, Nucleus_MKM and Nucleus_tMKM
     */
    class Track_KieferChatterjee : public Track
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Converts a particle (object of class Particle), passed by reference, in a track according to the Kiefer-Chatterjee amorphous track model.
            All data members are instantiated on the basis of the informations stored in the Particle object. (For a more detailed description of the instantiation of each member respectively look at its specific documentation).
         
            \param particle The particle generating the track in the medium, passed by reference.
            \param density The density of the medium expressed in \f$\frac{g}{cm^3}\f$. The default value is the density of water.
            \param t The time corresponding to the generation of the track in the target. The default value is 0. (See also the documentation of the data member #time).
         */
        Track_KieferChatterjee(const Particle &particle,
                               const double density,    // kg/dm^3 or g/cm^3
                               double t = 0.0);
        
        //! Returns a pointer to a new Track_KieferChatterjee object created as a copy of an existent one by means of the copy constructor.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         */
        virtual Track_KieferChatterjee* clone() const;
        
        //! Destructor.
        virtual ~Track_KieferChatterjee() {};
        
        //! Returns the ratio between the speed of the ion generating the track and the speed of light.
        /*!
            \return The ratio between the speed of the ion generating the track and the speed of light.
         
            \sa #beta
         */
        double getBeta() const {return beta;};
        
        //! Returns the distance from the center of the track, knowing the local dose deposited.
        /*!
            It's the inverse function of getLocalDose(). Since the maximum possible local dose deposited is: \f$d_{MAX}=\frac{k_p}{R_c^2}\f$, if d is greater than \f$d_{MAX}\f$ it returns -1 (nonsense value). Else it returns:
            \f[
                \sqrt{\frac{k_p}{d}}
            \f]
            where d is the local dose, \f$R_c\f$ is #r_core and \f$k_p\f$ is #k_p.
         
            \note If \f$d<d_{MIN}=\frac{k_p}{R_p^2}\f$ it returns #r_penumbra (\f$R_p\f$).
         
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
        
        //! Returns the #k_p value used to evaluate the dose to the core and to the penumbra.
        /*!
            \return The #k_p value used to evaluate the dose to the core and to the penumbra.
         */
        double getKp() const {return k_p;};
        
        //! Returns the LET in water of the particle generating the track expressed in MeV/um.
        /*!
            \return The LET in water of the particle generating the track expressed in MeV/um (according to the Bethe-Bloch formula).
         
            \sa #let
         */
        virtual double getLet() const {return let;};
        
        //! Returns the local dose (in Gy) deposited by the track at a certain distance from the track center (expressed in um).
        /*!
            The function evaluates some possible cases:
                - If the distance is smaller than the core radius it returns #dose_core
                - If the distance is greater than the penumbra radius it returns 0
                - Else it returns the dose to the penumbra calculated as: \f$\frac{k_p}{d^2}\f$, where \a d is the "distance" parameter
         
            \param distance The distance from the track center expressed in um.
         
            \return The local dose (in Gy) deposited by the track at a certain distance from the track center (expressed in um).
         
            \sa #r_core and #r_penumbra
         */
        virtual double getLocalDose(const double distance) const;
        
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
        
        //! Returns the dose (in Gy) deposited by the track evaluated as the radial integral of the track profile between \f$r_{min}\f$ and \f$r_{max}\f$ (expressed in um).
        /*!
            \param r_min Lower limit of integration, expressed in um.
            \param r_max Upper limit of integration, expressed in um.
         
            \return The dose (in Gy) deposited by the track evaluated as the radial integral of the track profile between \f$r_{min}\f$ and \f$r_{max}\f$.
         
            \warning The execution of the program will be terminated if incorrect limits of integration are chosen, that is:
                - If \f$r_{min}<0\f$
                - If \f$r_{max}<r_{min}\f$
         
            \sa #dose_core and getLocalDose()
         */
        virtual double getRadialIntegral(const double r_min,
                                         const double r_max) const;
        
        //! Returns the penumbra radius of the track expressed in um.
        /*!
            \return The penumbra radius of the track expressed in um.
         
            \sa #r_penumbra and #r_core
         */
        virtual double getRadius() const {return r_penumbra;};
        
        //! Returns the core radius of the track expressed in um.
        /*!
            \return The core radius of the track expressed in um.
         
            \sa #r_core and #r_penumbra
         */
        double getRCore() const {return r_core;};
        
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
        
        //! Returns the Barkas effective charge.
        /*!
            \return The Barkas effective charge.
         
            \sa #z_eff
         */
        double getZBarkas() const {return z_eff;};
        
        //! Saves the local dose deposited by the track at different distances from its center to reconstruct the track profile.
        /*!
            The function divides the track penumbra radius into 300 logarithmically spaced distances from the track center, for each of these distances it evaluates the local dose deposited calling the function getLocalDose(); during the process it saves each distance and the corresponding dose in a new file. This allow to reconstruct the track profile.
         
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
        
        //! The type of particle generating the track (e.g. Chemical symbol for ions: H, He, Li, ...).
        /*!
            \sa Particle::type
         */
        std::string particleType;
        
        //! The specific energy of the particle generating the track, expressed in MeV/u.
        double particleEnergy;
        
        //! The core radius of the track expressed in \f$\mu m\f$
        /*!
            It is initialized by means of the following relation:
            \f[
                R_{core}=\eta\cdot\beta_{Ion}
            \f]
            Where \f$\eta=0.0116\,\mu m\f$ is a constant value and \f$\beta_{Ion}\f$ represents the ratio between the speed of the ion and the speed of light.
         
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double r_core;
        
        //! The penumbra radius of the track expressed in \f$\mu m\f$
        /*!
            It is initialized by means of the following relation:
            \f[
                R_{penumbra}=\gamma\left(\frac{E}{A}\right)^\delta
            \f]
            Where \f$\gamma=0.0616\,\frac{\mu m}{MeV^\delta}\f$ and \f$\delta=1.7\f$ are constant values, E represents the energy of the ion and A his mass number (Particle::A).
         
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double r_penumbra;
        
        //! Constant dose in the core of the track expressed in Gy.
        /*!
            Following the Kiefer-Chatterjee amorphous model, the dose in the core is a constant value evaluated by means of the relation:
            \f[
                D_c=\frac{1}{\pi R_c^2}\left(\frac{LET}{\rho}-2\pi k_p\ln\left(\frac{R_p}{R_c}\right)\right)
            \f]
            Where \f$R_c\f$ and \f$R_p\f$ represent respectively the radius of the core and the penumbra; \f$\rho\f$ represents the density of the medium (Tracks::density), LET is the unrestricted linear energy transfer for the incident ion, \f$k_p\f$ (#k_p) is a function of \f$\beta_{Ion}\f$ (#beta) and \f$Z_{eff}\f$ (#z_eff), that is the ratio between the speed of the ion and the speed of light and the Barkas effective charge respectively.
         
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double dose_core;
        
        //! It's the ratio between the speed of the ion and the speed of light.
        /*!
            It is evaluated by the information stored in the Particle object using the following relation:
            \f[
                \beta=\sqrt{1-\left(\frac{E_k}{E_0}+1\right)^{-2}}
            \f]
            Where \f$E_k\f$ represents the kinetic energy of the particle (#e_c; Particle::e_c) and \f$E_0\f$ his rest energy (Particle::restEnergy).
         
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double beta;
        
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
        
        //! Barkas effective charge
        /*!
            The Barkas effective charge evaluated by means of the following relation:
            \f[
                Z_{eff}=Z\cdot\left(1-\exp\left(-\frac{125\beta}{Z^{2/3}}\right)\right)
            \f]
            Where Z represents the particle charge (Particle::charge) and #beta the ratio between the speed of the ion and the speed of light.
         
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double z_eff;
        
        //! Value used to evaluate the dose to the core and to the penumbra. It's a function of the Barkas effective charge and of the \f$\beta_{Ion}\f$ (#beta).
        /*!
            The value is evaluated by means of the following relation:
            \f[
                k_p=1.25\cdot 10^{-4}\cdot\left(\frac{Z_{eff}}{\beta_{Ion}}\right)^2
            \f]
            It is instantiated in the Constructor.
         
            \sa Track_KieferChatterjee(const Particle&, const double, double)
         */
        double k_p; // fattore di normalizzazione della traccia
        
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
         */
        double time;
        
        //! Constant value defined in the Kiefer-Chatterjee amorphous model.
        /*!
            It's equal to:
            \f[
                \delta=1.7
            \f]
         */
        static const double DELTA;
        
        //! Constant value defined in the Kiefer-Chatterjee amorphous model.
        /*!
            It's equal to:
            \f[
                \gamma=0.0616\,\frac{\mu m}{MeV^\delta}
            \f]
         */
        static const double GAMMA;
        
        //! Constant value defined in the Kiefer-Chatterjee amorphous model.
        /*!
            It's equal to:
            \f[
                \eta=0.0116\,\mu m
            \f]
         */
        static const double ETA;
    };

}


#endif /* TRACK_KIEFERCHATTERJEE_H */
