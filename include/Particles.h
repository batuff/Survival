#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <string>

namespace Survival {
    
    class Particle;

    //! This is a container class, used to group together Particle objects. It implements the structure and methods to manage a vector of particles.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2008
     
        This is a container class, used to group together Particle objects. It implements the structure and methods to manage a vector of particles, which is the only data member of the class. It provides also functionalities to select particles belonging to a specific region of space or corresponding to a certain category (e.g. Ions).
     
        \sa Particle and Tracks
     */
    class Particles
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            \param numberOfParticles The length of the vector or, likewise, the number of particles to be stored in the object.
         
            \sa Tracks
         */
        Particles(const int numberOfParticles = 0);
        
        //! Overload of the constructor. Instantiates and sets the object by loading a file where a spectrum of particles is defined.
        /*!
            The method reads the file "file_name" and, row by row, instantiates, sets and stores the particles in #particleVector.
         
            \param file_name A \c string identifying the file containing the spectrum of particles.
         
            \warning Almost no check will be done on the file.
         
            \sa loadSpectrum()
         */
        Particles(const std::string file_name);
        
        //! Destructor.
        ~Particles() {};
        
        //! Overload of the \c << operator to add a new particle at the end of the vector.
        /*!
            \param particle The particle to be added.
         */
        void operator<<(const Particle &particle);
        
        //! Overload of the \c << operator to add a vector of particles at the end of the vector.
        /*!
            \param particles The object containing the vector of particles to be added.
         */
        void operator<<(const Particles &particles);
        
        //! Overload of the \c [] operator to access at the n-th element of the vector.
        /*!
            \param index The position of the element in the vector.
         
            \return A reference to the element at the specified position in the vector.
         */
        Particle& operator[](const int index);
        
        //! Overload of the \c [] operator to access at the n-th element of the vector.
        /*!
            \param index The position of the element in the vector.
         
            \return A \c const reference to the element at the specified position in the vector.
         */
        const Particle& operator[](const int index) const;
        
        //! Returns the dose averaged LET of the vector of particles.
        /*!
            Evaluates and returns the dose averaged LET of the group of particles stored in the vector starting from the single LET and weight of each particle by means of the following relation:
            \f[
                LET_d=\frac{\sum w_i\cdot LET_i^2}{\sum w_i\cdot LET_i}
            \f]
         
            \return The dose averaged LET.
         */
        double getDoseAveragedLet() const;
        
        //---------------------
        //Particles* getIons();
        
        //! Selects and returns only the ions identified in the vector of particles.
        /*!
            The method check if both the charge and the mass number are greater than 0.
         
            \return A Particles object which is the subset of all ions stored in the original Particles object.
         
            \sa getIons(const int) and getIons(const int, const int)
         */
        Particles getIons();
        
        //! Selects and returns only the ions with a particular charge identified in the vector of particles.
        /*!
            Overload of the `getIons()` function to provide the possibility of selecting a particular charge value.
         
            \param charge The charge of the ions to be selected.
         
            \return A pointer to an object of the class Particles which is the subset of all ions, with a particular charge value, stored in the original Particles object.
         
            \sa getIons() and getIons(const int, const int)
         */
        Particles* getIons(const int charge);
        
        //! Selects and returns only the ions with a particular charge and mass number identified in the vector of particles.
        /*!
            Overload of the `getIons()` function to provide the possibility of selecting particular value of charge and mass.
         
            \param charge The charge of the ions to be selected.
            \param A The mass number of the ions to be selected.
         
            \return A pointer to an object of the class Particles which is the subset of all ions, with a particular value of charge and mass, stored in the original Particles object.
         
            \sa getIons() and getIons(const int, const int)
         */
        Particles* getIons(const int charge,
                           const int A);
        
        //! Returns the mean LET of the vector of particles, expressed in keV/um.
        /*!
            Evaluates and returns the mean LET of the group of particles stored in the vector starting from the single LET and weight of each particle by means of the following relation:
            \f[
                \langle LET\rangle=\frac{\sum w_i\cdot LET_i}{\sum w_i}
            \f]
         
            \return The mean LET, expressed in keV/um.
         
            \sa Particle::let
         */
        double getMeanLet() const;
        
        //! Returns a \c string identifying the file containing the spectrum of particles.
        /*!
            \return A \c string identifying the file containing the spectrum of particles (#spectrum_file).
         */
        std::string getSpectrumFile() const {return spectrum_file;};
        
        //! Returns the total LET of the vector of particles.
        /*!
            \return The total LET, evaluated as the sum of the LET of each particle.
         
            \sa Particle::let
         */
        double getTotalLet() const;
        
        //! Returns the total weight of the vector of particles.
        /*!
            \return The total weight, evaluated as the sum of the weight of each particle.
         
            \sa Particle::weight
         */
        double getTotalWeight() const;
        
        //! Selects and returns only the particles with coordinates between \f$\left[ x_{min},\,x_{max}\right]\f$ and \f$\left[ y_{min},\,y_{max}\right]\f$.
        /*!
            \param x_min The minimum value of the \c x coordinate, expressed in mm.
            \param x_max The maximum value of the \c y coordinate, expressed in mm.
            \param y_min The minimum value of the \c x coordinate, expressed in mm.
            \param y_max The maximum value of the \c y coordinate, expressed in mm.
         
            \return A pointer to an object of the class Particles which is a subset of the original Particles object.
         
            \sa getIons() and getWithDistanceBetween()
         */
        Particles* getWithCoordinatesBetween(const double x_min,
                                             const double x_max,
                                             const double y_min,
                                             const double y_max);
        
        //! Selects and returns only the particles with a distance from the origin between \f$\left[ distance_{min},\,distance_{max}\right]\f$.
        /*!
            Selects all particles in an annulus between \f$\left[ distance_{min},\,distance_{max}\right]\f$.
         
            \param distance_min The minimum distance from the origin, expressed in mm.
            \param distance_max The maximum distance from the origin, expressed in mm.
         
            \return A pointer to an object of the class Particles which is a subset of the original Particles object.
         
            \sa getIons() and getWithCoordinatesBetween()
         */
        Particles* getWithDistanceBetween(const double distance_min,
                                          const double distance_max);
        
        //! Loads a spectrum from a file and stores it in particleVector.
        /*!
            \param file_name A \c string identifying the file containing the spectrum of particles.
         */
        void loadSpectrum(const std::string file_name);
        
        //! For each particle of particle of particle vector, if its LET or energy is undefined the method tries to set it using the Bethe-Bloch formula, assuming it's an ion.
        /*!
            The actual LET of each particle is assumed to be expressed in keV/um, while the kinetic energy in MeV.
         
            \warning This function is NOT general and thought only for a precise particular purpose. Execution of the program will be terminated if an element of the vector isn't an ion, and there aren't any checks on other particle features.
         */
        void reconstructIonLETandEnergy();
        
        //! Sets the name of the file containing the spectrum of particles.
        /*!
            \param file_name A \c string identifying the file containing the spectrum of particles.
         */
        void setSpectrumFile(const std::string file_name){spectrum_file=file_name;};
        
        //! Returns the size of the vector of particles.
        /*!
            \return The size of the vector of particles: #particleVector
         */
        int size() const;
        
    private:
        
        //! The vector of particles.
        /*!
            \sa Particle
         */
        std::vector< Particle > particleVector;
        
        //! A \c string identifying the file containing the spectrum of particles.
        std::string spectrum_file;
    };
    
}


#endif /* PARTICLES_H */
