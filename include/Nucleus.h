#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <string>
#include <vector>

namespace Survival {
    
    class Tracks;
    class Track;
    class CellLine;

    //! Linked to a specific CellLine object, this class defines the structure of the cellular nucleus to be irradiated and evaluates the dose absorbed during the interaction with a track.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007--2015
     
        The idea is that this class receives a set of tracks contained in a given Tracks object and corresponding to a certain spatial configuration of ion transversals, with the objective of:
            - superimposing the tracks in order to compute the composite local dose distribution;
            - transforming local doses in local number of lethal events by queries to the CellLine object;
            - integrating the local dose and the local number of lethal events
            - giving back the mean dose and mean survival estimation
     
        Since these tasks can be in principle accomplished in several ways, this class has been declared pure virtual. The present derived implementation are Nucleus_MKM, Nucleus_tMKM, Nucleus_Integral, Nucleus_Integral_t, Nucleus_Pixel and Nucleus_MonteCarlo classes which implements the structure of the nucleus defined in the LEM I, II, III and MKM models.
     
        For a more detailed description see the derived classes.
     
        \sa Nucleus_MKM, Nucleus_tMKM, Nucleus_Integral, Nucleus_Integral_t, Nucleus_Pixel and Nucleus_MonteCarlo
     */
    class Nucleus
    {
    public:
        
        //! Constructor of a pure virtual class (empty).
        Nucleus() {};
        
        //! Destructor of a pure virtual class (empty).
        virtual ~Nucleus() {};
        
        //! Declaration of the virtual function getNucleusDoses (for a more detailed description see the derived classes).
        virtual void addNucleusDoses(Nucleus &){};
        
        //! Declaration of the pure virtual function getCellType (for a more detailed description see the derived classes).
        virtual std::string getCellType() const = 0;
        
        //! Declaration of the pure virtual function cleanNucleus (for a more detailed description see the derived classes).
        virtual void cleanNucleus() = 0;
        
        //! Declaration of the pure virtual function clone (for a more detailed description see the derived classes).
        virtual Nucleus* clone(const CellLine&) = 0;
        
        //! Declaration of the pure virtual function distributeDose (for a more detailed description see the derived classes).
        virtual void distributeDose(const Track &track) = 0;
        
        //! Declaration of the pure virtual function distributeDose (for a more detailed description see the derived classes).
        virtual void distributeDose(const Tracks &tracks) = 0;
        
        //! Declaration of the virtual function getAC (for a more detailed description see the derived classes).
        //virtual double getAC() const {return 0.;};
        
        //! Declaration of the virtual function getDomainRadius (for a more detailed description see the derived classes).
        virtual double getDomainRadius() {return 0.;};
        
        //! Declaration of the virtual function getDosesandLethals (for a more detailed description see the derived classes).
        virtual void getDosesAndLethals(std::vector<double> &,std::vector<double> &,std::vector<double> &,std::vector<double> &){}; // funzione farlocca
        
        //! Declaration of the pure virtual function getDoseAndSurvival (for a more detailed description see the derived classes).
        virtual void getDoseAndSurvival(double &dose,
                                        double &doseUncertainty,
                                        double &survival,
                                        double &survivalUncertainty) const = 0;
        
        //! Declaration of the pure virtual function getInNucleusCount (for a more detailed description see the derived classes).
        virtual int getInNucleusCount() const = 0;
        
        //! Declaration of the pure virtual function getIntersectionCount (for a more detailed description see the derived classes).
        virtual int getIntersectionCount() const = 0;
        
        //! Declaration of the virtual function getNumberOfDomains (for a more detailed description see the derived classes).
        virtual int getNumberOfDomains() {return 1;};
        
        //! Declaration of the pure virtual function getPosition (for a more detailed description see the derived classes).
        virtual void getPosition(double &returnX,
                                 double &returnY) const = 0;
        
        //! Declaration of the pure virtual function getRadius (for a more detailed description see the derived classes).
        virtual double getRadius() const = 0;
    };
    
}


#endif /* NUCLEUS_H */
