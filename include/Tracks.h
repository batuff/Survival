#ifndef TRACKS_H
#define TRACKS_H

#include <vector>
#include <string>

namespace Survival {
    
    class Track;
    class Particles;

    //! This is a container class for Track objects; it implements the structure and methods to manage a vector of tracks.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2008
     
        This is a container class for Track objects; it implements the structure and methods to manage a vector of tracks. It can be created directly from a Particles object, specifying a unique Track type, or it can be loaded with single Track objects, making use of polymorphism.
     
        \sa Track and Particles
     */
    class Tracks
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Converts a vector of particles (object of class Particles), passed by reference, in a vector of tracks, depending on the parametrization chosen for the Track model.
            The possible choices for the parametrization are:
                - Scholz2000
                - Elsasser2007
                - Elsasser2008
                - KieferChatterjee
         
            \warning The execution of the program will be terminated if an inexistent parametrization is chosen.
         
            \param particles The reference to a Particles object, that is the vector of particles generating the tracks.
            \param trackType The parametrization chosen for the Track model.
            \param massDensity The density of the medium expressed in \f$\frac{g}{cm^3}\f$. The default value is the density of water.
         
            \sa Tracks(const int, const double), Track_Scholz2000, Track_Elsasser2007, Track_Elsasser2008 and Track_KieferChatterjee
         */
        Tracks(const Particles &particles,
               const std::string trackType,
               const double massDensity = 1.0);
        
        //! Instantiates and sets the object. Overload of the constructor.
        /*!
            \param numberOfTracks The number of tracks to be stored in the vector
            \param massDensity The density of the medium expressed in \f$\frac{g}{cm^3}\f$. The default value is the density of water.
         
            \sa Tracks(const int, const double)
         */
        Tracks(const int numberOfTracks = 0,
               const double massDensity = 1.0);
        
        //! Destructor.
        /*!
            The destructor calls the function eraseAll() that deletes each element of the vector and the vector itself.
         
            \sa eraseAll()
         */
        ~Tracks();
        
        //! Overload of the \c << operator to add a new track at the end of the vector.
        /*!
            \param track The track to be added.
         */
        void operator<<(const Track &track);
        
        //! Overload of the \c << operator to add a vector of tracks at the end of the vector.
        /*!
            \param tracks The object containing the vector of tracks to be added.
         */
        void operator<<(const Tracks &tracks);
        
        //! Overload of the \c [] operator to access at the n-th element of the vector.
        /*!
            \param index The position of the element in the vector.
         
            \return A \c const reference to the element at the specified position in the vector.
         */
        const Track& operator[](const int index) const;
        
        //! Deletes each element of the vector and erases the vector itself.
        void eraseAll();
        
        //! Returns the #density of the medium in \f$\frac{g}{cm^3}\f$.
        /*!
            \return The #density of the medium in \f$\frac{g}{cm^3}\f$.
         */
        double getDensity() const {return density;};
        
        //! Evaluates and returns the dose averaged LET of the vector of tracks, expressed in keV/um.
        /*!
            It's evaluated as:
            \f[
                LET_d=\frac{\sum_i\,LET_i^2\,w_i}{\sum_i\,LET_i\,w_i}
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector.
         
            \return The dose averaged LET of the vector of tracks (#trackVector), expressed in keV/um.
         */
        double getDoseAveragedLet() const;
        
        //! Evaluates and returns the mean energy of the vector of tracks, expressed in MeV.
        /*!
            It's evaluated as:
            \f[
                \langle E\rangle=\frac{\sum_i\,E_i\,w_i}{\sum_i\,w_i}
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector.
         
            \return The mean energy of the vector of tracks (#trackVector), expressed in MeV.
         
            \sa Track::getKineticEnergy()
         */
        double getMeanEnergy() const;
        
        //! Evaluates and returns the mean LET of the vector of tracks, expressed in keV/um.
        /*!
            It's evaluated as:
            \f[
                \langle LET\rangle=\frac{\sum_i\,LET_i\,w_i}{\sum_i\,w_i}
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector.
         
            \return The mean LET of the vector of tracks (#trackVector), expressed in keV/um.
         
            \sa Particles::getMeanLet().
         */
        double getMeanLet() const;
        
        //! Evaluates and returns the standard deviation of the dose averaged LET of the vector of tracks, expressed in keV/um.
        /*!
            It's evaluated as:
            \f[
                \sigma_{LET_d}=\left(\frac{\sum_i(\,LET_i-LET_d)^2\,LET_i\,w_i}{\sum_i\,LET_i\,w_i}\right)^{1/2}
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector and \f$LET_d\f$ the dose averaged LET.
         
            \return The standard deviation of the dose averaged LET of the vector of tracks (#trackVector), expressed in keV/um.
         
            \sa getDoseAveragedLet()
         */
        double getSigmaDoseAveragedLet() const;
        
        //! Evaluate and returns standard deviation of the mean energy of the vector of tracks, expressed in MeV.
        /*!
            It's evaluated as:
            \f[
                \sigma_{\langle E\rangle}=\left(\frac{\sum_i (E_i-\langle E\rangle)^2\,w_i}{\sum_i\,w_i}\right)^2
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector and \f$\langle E\rangle\f$ the mean energy.
         
            \return The standard deviation of the mean energy of the vector of tracks (#trackVector), expressed in MeV.
         
            \sa getMeanEnergy()
         */
        double getSigmaMeanEnergy() const;
        
        //! Evaluate and returns standard deviation of the mean LET of the vector of tracks, expressed in keV/um.
        /*!
            It's evaluated as:
            \f[
                \sigma_{\langle LET\rangle}=\left(\frac{\sum_i (LET_i-\langle LET\rangle)^2\,w_i}{\sum_i\,w_i}\right)^2
            \f]
            where \f$w_i\f$ indicates the weight of the i-th track in the vector and \f$\langle LET\rangle\f$ the mean LET.
         
            \return The mean LET of the vector of tracks (#trackVector), expressed in keV/um.
         
            \sa Particles::getMeanLet
         */
        double getSigmaMeanLet() const;
        
        //! Returns a \c string identifying the file containing the spectrum of particles to be converted in tracks.
        /*!
            \return A \c string identifying the file containing the spectrum of particles to be converted in tracks.
         */
        std::string getSpectrumFile() const {return spectrum_file;};
        
        //! Returns the total weight of the group of particles generating the vector of tracks.
        /*!
            \return The total weight of the group of particles generating the vector of tracks, evaluated as the sum of the weight of each particle.
         
            \sa Particles::getTotalWeight
         */
        double getTotalWeight() const;
        
        //! Check if the vector is monoenergetic.
        /*!
            \return A boolean value identifying if #trackVector is monoenergetic, that is \f$size = 1\f$.
         
            \sa size()
         */
        bool isMonoenergetic() const;
        
        //! Sets the #density of the medium.
        /*!
            \param d The chosen value of #density of the medium to be set, expressed in \f$\frac{g}{cm^3}\f$
         */
        void setDensity(const double d){density = d;};
        
        //! Returns the size of the vector of tracks.
        /*!
            \return The size of the vector of tracks: #trackVector.
         */
        int size() const;
        
    private:
        
        //! A dynamic vector of Track pointers.
        std::vector< Track* > trackVector;  //for the use of polymorphism
        
        //! The density of the medium expressed in \f$\frac{g}{cm^3}\f$.
        double density;
        
        //! A \c string identifying the file containing the spectrum of particle to be converted in tracks.
        std::string spectrum_file;
    };
    
}


#endif /* TRACKS_H */
