#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>

namespace Survival {

    //! This class defines the object "particle".
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007
     
        This class defines the object "particle". It is used as a C++ struct to contain, for a certain particle in a given position in space, recorded characteristics like type, charge, mass number, kinetic energy, LET in water and position. It has no member functions; it has only data members that identify the characteristics of the particle. Note that all data members are defined \c public.
     */
    class Particle
    {
    public:
        
        //! The type of particle (e.g. Chemical symbol for ions: H, He, Li, ...).
        /*!
            This data member identify the particle type. Note that, if the particle is an ion, only ions with atomic number \f$Z<=10\f$ are supported.
         */
        std::string type;
        
        //! The charge of the particle expressed in elementary charge units.
        int charge;
        
        //! The rest energy of the particle expressed in MeV.
        double restEnergy;
        
        //! The mass number of the particle.
        int A;
        
        //! The kinetic energy of the particle expressed in MeV.
        double e_c;
        
        //! The LET in water of the particle expressed in MeV/um (according to the Bethe-Bloch formula).
        double let;
        
        //! The particle position (\c x coordinate) referred to the beam axis, expressed in mm.
        double x;
        
        //! The particle position (\c y coordinate) referred to the beam axis, expressed in mm.
        double y;
        
        //! The particle position (\c z coordinate) referred to the depth of penetration in matter.
        double z;
        
        //! The weight of this particular particle type in the beam. Useful in the case of "mixed fields".
        double weight;
    };
    
}


#endif /* PARTICLE_H */
