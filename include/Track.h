#ifndef TRACK_H
#define TRACK_H

#include <string>

namespace Survival {

    //! Constructed on the base of an ion Particle object, this class represents the local dose distribution around that ion.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007--2015
     
        Constructed on the base of an ion Particle object, this class represents the local dose distribution around that ion. It is a pure virtual class, defined by the inherited Track_Scholz2000, Track_Elsasser2007, Track_Elsasser2008 and Track_KieferChatterjee classes, which implement the track models of LEM I, II, III and MKM respectively.
     */
    class Track
    {
    public:
        
        //! Constructor of a pure virtual class (empty).
        Track() {};
        
        //! Destructor of a pure virtual class (empty).
        virtual ~Track() {};
        
        //! Declaration of the pure virtual function clone (for a more detailed description see the derived classes).
        virtual Track* clone() const = 0;
        
        //! Declaration of the pure virtual function getDistance (for a more detailed description see the derived classes).
        virtual double getDistance(const double localDose) const = 0;  // Gy -> um
        
        //! Declaration of the pure virtual function getKineticEnergy (for a more detailed description see the derived classes).
        virtual double getKineticEnergy() const = 0;  // MeV
        
        //! Declaration of the pure virtual function getLet (for a more detailed description see the derived classes).
        virtual double getLet() const = 0;  // MeV/um
        
        //! Declaration of the pure virtual function getLocalDose (for a more detailed description see the derived classes).
        virtual double getLocalDose(const double distance) const = 0;  // distance from track center (um) -> Gy
        
        //! Declaration of the pure virtual function getParticleEnergy (for a more detailed description see the derived classes).
        virtual double getParticleEnergy() const = 0;
        
        //! Declaration of the pure virtual function getParticleType (for a more detailed description see the derived classes).
        virtual std::string getParticleType() const = 0;
        
        //! Declaration of the pure virtual function getPosition (for a more detailed description see the derived classes).
        virtual void getPosition(double &returnX,
                                 double &returnY) const = 0;  // track position referred to the beam axis (mm)
        
        //! Declaration of the pure virtual function getRadialIntegral (for a more detailed description see the derived classes).
        virtual double getRadialIntegral(const double r_min,
                                         const double r_max) const = 0;  // extremes of integration (um)
        
        //! Declaration of the pure virtual function getRadius (for a more detailed description see the derived classes).
        virtual double getRadius() const = 0;  // um
        
        //! Declaration of the pure virtual function getTime (for a more detailed description see the derived classes).
        virtual double getTime() const = 0;
        
        //! Declaration of the pure virtual function getWeight (for a more detailed description see the derived classes).
        virtual double getWeight() const = 0;
        
        //! Declaration of the pure virtual function saveTrack (for a more detailed description see the derived classes).
        virtual std::string saveTrack() const = 0;
        
        //! Declaration of the pure virtual function setPosition (for a more detailed description see the derived classes).
        virtual void setPosition(const double x,
                                 const double y) = 0;  // track position referred to the beam axis (mm)
        
        //! Declaration of the pure virtual function setTime (for a more detailed description see the derived classes).
        virtual void setTime(double t) = 0;
    };
    
}


#endif /* TRACK_H */
