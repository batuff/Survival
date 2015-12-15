#ifndef USEFUL_FUNCTIONS_H
#define USEFUL_FUNCTIONS_H

#include <vector>
#include <string>

namespace Survival {

    //! Portable wrapper for mkdir. Internally used by mkdir().
    /*!
        \param path The full path of the directory to create.
     
        \return 0 on success, otherwise -1.
     
        \sa folder_exists() and mkdir()
     */
    int _mkdir(const char* path);
    
    //! Returns the kinetic energy associated to the value of LET imposed for a specified ion.
    /*!
        The method opens and read the data file correspondent to the ion chosen until it finds in the table the nearest neighbors of the LET imposed, then it interpolates these values to get the correct kinetic energy required.
    
        \param ionType A \c string identifying the ion (chemical symbol).
        \param let_imposed The LET (in MeV/um) correspondent to the required kinetic energy.
     
        \return The kinetic energy in MeV/u associated to the value of LET imposed for a specified ion.
     
        \warning The execution of the program will be terminated if:
            - The data file of the ion chosen doesn't exist.
            - The environment hasn't been set correctly.
            - The LET imposed is out of the available range for the ion specified.
     
        \sa BetheBlochTable and betheBloch_Srim()
     */
    double betheBloch_inv_Srim(std::string ionType, double let_imposed);
    
    //! Returns the kinetic energy associated to the value of kinetic energy imposed for a specified ion.
    //!
    /*!
        The method makes use of a BetheBlochTable object, managed by an STL map (in pairing with a \c string identifying ion type chosen). The idea is to avoid reading the file every time the function BetheBlochTable::GetLet() is called. For this reason the map is defined \c static, and BetheBlochTable() (which reads the correspondent file) is called only when an ion is used for the first time.
     
        \param ionType A \c string identifying the ion (chemical symbol).
        \param e_c_imposed The kinetic energy (in MeV/u) correspondent to the required kinetic energy.
     
        \return The LET in MeV/um associated to the value of kinetic energy imposed for a specified ion.
     
        \sa betheBloch_inv_Srim()
     */
    double betheBloch_Srim(std::string ionType, double e_c_imposed);
    
    
    //! Perform a fit on a data set of doses and associated survival according to the linear quadratic formula.
    /*!
        The method perform a fit (making use of the GSL library) on the data set of doses (as \c x coordinate) and survivals (as \c y coordinate), weighted on the survival uncertainty, according to the linear quadratic formula:
        \f[
            S=\exp(-\alpha D-\beta D^2)
        \f]
        The method returns the LQ parameters \f$\alpha\f$ and \f$\beta\f$ and some other useful informations on the fit (such as the \f$\chi^2\f$ or the uncertainties on the estimated parameters) by overwriting the corresponding variables passed by reference.
     
        \param dose The vector of doses (in Gy) to be used in the fit.
        \param survival The vector of survivals to be used in the fit.
        \param survivalUncertainty The vector of uncertainties associated to the evaluated survival to be used as weights in the fit.
        \param alpha The linear parameter \f$\alpha\f$ estimated, expressed in \f$Gy^{-1}\f$.
        \param alphaUncertainty The uncertainty associated to the \f$\alpha\f$ parameter (in \f$Gy^{-1}\f$).
        \param beta The quadratic parameter \f$\beta\f$ estimated, expressed in \f$Gy^{-2}\f$.
        \param betaUncertainty The uncertainty associated to the \f$\beta\f$ parameter (in \f$Gy^{-2}\f$).
        \param chiSquared The value of the \f$\chi^2\f$ obtained from the fit.
        \param incompleteGammaQ The normalized incomplete Gamma Function \f$Q(a,x)=\frac{1}{\Gamma(a)}\int_x^{+\infty}t^{a-1}e^{-t}dt\f$ (\f$a>0\f$; \f$x>=0\f$).
     
        \sa Calculus::slow_meanDose_meanSurvival()
     */
    void fit_LQ(std::vector<double> dose,
                std::vector<double> survival,
                std::vector<double> survivalUncertainty,
                double &alpha, double &alphaUncertainty,
                double &beta, double &betaUncertainty,
                double &chiSquared, double &incompleteGammaQ);
    
    //! Checks if a folder exists.
    /*!
        \param foldername The path to the folder to check.
     
        \return \c true if the folder exists, \c false otherwise.
     
        \sa _mkdir() and mkdir()
     */
    bool folder_exists(std::string foldername);
    
    //! Recursive, portable wrapper for mkdir.
    /*!
        \param path The full path of the directory to create.
     
        \return 0 on success, otherwise -1.
     
        \sa _mkdir() and folder_exists()
     */
    int mkdir(const char* path);
    
    //! Parses the input arguments in the \c main to set the calculation parameters.
    /*!
        Parses the input arguments passed by the user in order to set the calculation parameters for the simulation. If an incorrect setting occurs, the program will be terminated and the method display an hint to the user to fix the problem.
     
        \note All parameters, with the exception of \c argc and \c argv, are passed by reference to be overwritten.
     
        \param argc The number of input arguments.
        \param argv A pointer to the array of \c chars, that is the input arguments set by the user.
        \param cellType A \c string identifying the name of the cell line.
        \param model A \c string identifying the model to be used in the simulation. "LEMI", "LEMII", "LEMIII", "MKM" and "tMKM" supported.
        \param trackType A \c string defining the type of Track to be used in the simulation. "KieferChatterjee", "Scholz2000", "Elsasser2007" and "Elsasser2008" supported.
        \param parametrizationType
        \param calculusType A \c string identifying the way to perform
        \param precision Fix the ending condition of the Monte Carlo simulation.
        \param parallelismType The number of threads needed to be used in the simulation. It has to be a ...
        \param doses A vector of \c double containing MIN, MAX and number of doses to be simulated in order to construct the survival curve and extrapolate the LQ parameters.
        \param parameter_name A vector of \c string identifying the name of the model parameters correspondent to the model chosen.
        \param MKM_alpha0 The MKM \f$\alpha_0\f$ parameter, expressed in \f$Gy^{-1}\f$.
        \param MKM_beta0 The MKM \f$\beta_0\f$ parameter, expressed in \f$Gy^{-2}\f$.
        \param MKM_rNucleus The radius of the MKM nucleus.
        \param MKM_rDomain The radius of domains in the MKM nucleus.
        \param tMKM_ac The time constant representing the cellular repair kinetics, expressed in \f$h^{-1}\f$.
        \param LEM_alpha0 The LEM \f$\alpha_0\f$ parameter, expressed in \f$Gy^{-1}\f$.
        \param LEM_beta0 The LEM \f$\beta_0\f$ parameter, expressed in \f$Gy^{-2}\f$.
        \param LEM_rNucleus The radius of the LEM nucleus.
        \param LEM_Dt The transition dose beyond which the standard linear quadratic parametrization is no more valid according to the LEM parametrization.
        \param ionType A \c string identifying the ion generating the track (Chemical symbol: H, He, Li, ...). Ions until "Ne" are supported.
        \param particleA The mass number of the ion chosen.
        \param particleZ The atomic number of the ion chosen.
        \param trackMode A \c string identifying the modality to pass the vector of particles in the mixed fields case. "histogram" or "random" are supported.
        \param energyType A \c string identifying if input values are energies or LETs. "LET" or "energy" supported.
        \param energies A vector of \c double to store energies ot LETs to be used.
        \param nFraction A vector of \c double containing the number of fractions to be delivered in a fractionated treatment (MIN, MAX and STEP). Compatible with tMKM model.
        \param timeSpacing A vector of \c double containing the time spacing between fractions of a fractionated treatment (MIN, MAX and STEP). Compatible with tMKM model.
        \param fracDeliveryTime A vector of \c double containing the fraction delivery time of a fractionated treatment (MIN, MAX and STEP). Compatible with tMKM model.
        \param saveAlphaBeta A boolean value identifying if LQ data are to be saved.
        \param saveMeans A boolean value identifying if mean doses and survivals are to be saved.
        \param saveCell A boolean value identifying if the data of dose absorbed and survival observed by each single cell irradiated in the Monte Carlo simulation are to be saved.
        \param projectName The name of the project, chosen by the user.
        \param mono A boolean value identifying if the radiation is monoenergetic or not.
        \param spectrum_file A \c string identifying the name of the file containing the spectrum of energies to be used in the simulation (complete with its absolute or relative path).
     */
    void parse(int argc,
               char* argv[],
               std::string& cellType,
               std::string& model,
               std::string& trackType,
               std::string& parametrizationType,
               std::string& calculusType,
               double& precision,
               int& parallelismType,
               std::vector<double>& doses,
               std::vector<std::string>& parameter_name,
               double& MKM_alpha0,
               double& MKM_beta0,
               double& MKM_rNucleus,
               double& MKM_rDomain,
               double& tMKM_ac,
               double& LEM_alpha0,
               double& LEM_beta0,
               double& LEM_rNucleus,
               double& LEM_Dt,
               std::string& ionType,
               int& particleA,
               int& particleZ,
               std::string& trackMode,
               std::string& energyType,
               std::vector<double>& energies,
               int& nFraction,
               double& timeSpacing,
               double& fracDeliveryTime,
               bool& saveAlphaBeta,
               bool& saveMeans,
               bool& saveCell,
               std::string& projectName,
               bool& mono,
               std::string& spectrum_file);
    
    //! Display an hint to the user to correctly use the executable.
    void Usage();
    
}

#endif /* USEFUL_FUNCTIONS_H */
