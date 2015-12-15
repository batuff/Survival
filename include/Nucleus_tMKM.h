#ifndef NUCLEUS_TMKM_H
#define NUCLEUS_TMKM_H

#include "Nucleus.h"

namespace Survival {
    
    class Nucleus_Pixel;
    class Nucleus_Integral_t;

    //! Inherited from the Nucleus pure virtual class, it implements the nucleus defined in the Monte Carlo temporal reformulation of the MKM model (MCt-MKM).
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2015
     
        Similar in its structure to the Nucleus_MKM class, it provides some method to manage the temporal structure of the irradiation to support the MonteCarlo temporal-Microdosimetric Kinetic Model (MCt-MKM, \ref MCt-MKM "1"). It keeps track of the history of the irradiation, associating to each dose deposited a precise temporal instant. The total number of lethal events observed and the associated cellular survival are evaluated considering also the repaired kinetics of the cell.
     
        \anchor MCt-MKM 1. L. Manganaro, G. Russo, A. Attili, "Advanced modeling of the temporal effect in particle therapy: from radiobiological evaluation to treatment planning", \a Medical \a Physics, (Submitted).
     */
    class Nucleus_tMKM : public Nucleus
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            When the constructor is called, it instantiates the object creating a hexagonal structure of circular domains, using the informations stored in the CellLine reference object, such as the radius of the nucleus and the radius of the single domain. This is made by calling the createDomains() function.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa createDomains(), cleanNucleus(), Nucleus_tMKM(const CellLine&, double, int, const double, const double) and Nucleus_MKM
         */
        Nucleus_tMKM(const CellLine &cellLineRef,
                     const double xPosition = 0.0,
                     const double yPosition = 0.0);
        
        //! Instantiates and sets the object. Overload of the constructor.
        /*!
            Provide the possibility to specify also the total number of domains and the radius of each domain.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param domainRadius The radius of the single domain in the tMKM nucleus, expressed in um.
            \param numberOfDomains The number of domains that constitute the tMKM nucleus.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa createDomains(), Nucleus_Integral_t and Nucleus_tMKM(const CellLine&, const double, const double)
         */
        Nucleus_tMKM(const CellLine &cellLineRef,
                     double domainRadius,
                     int numberOfDomains,
                     const double xPosition = 0.0,
                     const double yPosition = 0.0);
        
        //! Destructor.
        /*!
            The destructor deletes the domainCell, domains and each domains[i-th] pointers created when the object is instantiated.
         
            \sa createDomains()
         */
        virtual ~Nucleus_tMKM();
        
        //! Adds a constant value of dose absorbed in each domain of the nucleus in a specific instant.
        /*!
            The method calls systematically for each domain the function Nucleus_Integral_t::addBackgroundDose(), that is the override of this method itself in the Nucleus_Integral_t class. The result is to add a constant value of dose absorbed in each domain (at a specific instant) and consequently in the whole nucleus.
         
            \param dose The dose to be added expressed in Gy.
            \param t The time associated to the dose added, expressed in hours.
         */
        void addBackgroundDose(const double dose,
                               const double t);
        
        //! Resets to zero all counters (#inNucleusCount and #intersectionCount) and doses in the nucleus and his domain.
        /*!
            Resets to zero #inNucleusCount and #intersectionCount and calls systematically for each domain the function Nucleus_Integral_t::cleanNucleus(), that this the override of this method itself in the Nucleus_Integral_t class.
         */
        virtual void cleanNucleus();
        
        //! Returns a pointer to a new Nucleus_tMKM object. It not really a clone but a new clean object.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         
            \note To create a real clone of another nucleus, a better implementation of the copy constructor is needed.
         */
        virtual Nucleus_tMKM* clone(const CellLine&);
        
        //! When an interaction between a Particle and the tMKM nucleus occurs, the method increases #inNucleusCount and #intersectionCount counters and proceeds to distribute the dose in each domain.
        /*!
            This function is the first step to evaluate the dose deposited by the radiation in the nucleus or, better, in each domain of the tMKM nucleus. It checks if any interaction exists between Nucleus and Particle (that is the Track generated by the particle) looking at their positions and radius. If it's true, it increases the respective counters (#inNucleusCount and #intersectionCount) and calls the method Nucleus_Integral_t::distributeDose() in a \c for loop over the total number of domains.
         
            \param track A reference to the Track generated by the particle in the nucleus.
         */
        virtual void distributeDose(const Track &track);
        
        //! Overload of distributeDose(const Track &track) to manage a Tracks object, it simply calls distributeDose(const Track &track) for each track of the container.
        /*!
            Since the Tracks class is a container for Track objects, this method calls distributeDose(const Track &track) in a \c for loop over each track contained in the Tracks object.
         */
        virtual void distributeDose(const Tracks &tracks);
        
        //----------------------------------------------------
        //double evaluateG(double alpha_ion, double beta_ion);
        
        //! Returns the name of the cell line to which the nucleus belongs.
        /*!
            \return A \c string corresponding to the name of the cell line to which the nucleus belongs, getting the information by #cellLine.
         */
        virtual std::string getCellType() const;
        
        //! Returns the radius of the domain corresponding to the CellLine to which the tMKM nucleus belongs.
        /*!
            \return #domainRadius That is the radius of the domain corresponding to the CellLine to which the tMKM nucleus belongs.
         */
        virtual double getDomainRadius() {return domainRadius;};
        
        //! Returns (overwriting parameters passed by reference) the total dose absorbed by the nucleus and the associated survival, evaluated taking into account the time structure of the irradiation, with respective uncertainties.
        /*!
            The function calls, in a \c for loop over the total number of domains, the methods Nucleus_Integral_t::getDoses() and Nucleus_Integral_t::getTimes() to get the complete history of doses deposited in the nucleus by the radiation, each dose associated to a specific instant. Then it evaluates for each domain the total number of lethal events observed, taking into account the time structure of the irradiation, by means of the following relation:
            \f[
                L_d=-\alpha_d\left(\sum_{i=1}^N z_i\right) -\beta_d\left(\sum_{i=1}^N\right)^2 -2\beta\sum_{i=1}^{N-1}\sum_{j=i+1}^{N}\left[1-\exp\left(-(a+c)(t_j-t_i)\right)\right]z_i\,z_j
            \f]
            where N represents the total number of interaction events in the domain, the sum (a+c) represents the time constant characteristics of the cellular repair kinetics and \f$z_i\f$ and \f$t_i\f$ represent the i-th element of the vectors of times and doses absorbed respectively.
         
            The total number of lethal events in the nucleus is evaluated as the sum of the ones observed in each domain and the cellular survival is evaluated as a negative exponential function of the total number of lethal events (according to the poissonian statistics):
            \f[
                S=\exp(-L_{TOT})
            \f]
         
            Finally the method overwrite the total dose absorbed by the nucleus and the cellular survival with respective uncertainties.
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival associated to the dose absorbed by the nucleus (with its specific time structure), passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and survival, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         */
        virtual void getDoseAndSurvival(double &dose,
                                        double &doseUncertainty,
                                        double &survival,
                                        double &survivalUncertainty) const;
        
        //! Return the dose absorbed in the indexOfDomain-th domain, expressed in Gy.
        /*!
            The function calls the method Nucleus_Integral_t::getDose() from the indexOfDomain-th domain.
         
            If an incorrect index is selected (e.g. greater than the total number of domains) the dose absorbed is set to -1.
         
            \param indexOfDomain The index associated to the domain.
         
            \return The dose absorbed in the indexOfDomain-th domain, expressed in Gy.
         */
        double getDoseForDomain(int indexOfDomain) const;
        
        //! Returns the number of times that the nucleus has been crossed through by a Particle.
        /*!
            \return #inNucleusCount The number of times that the nucleus has been crossed through by a Particle.
         
            \sa distributeDose(const Track&)
         */
        virtual int getInNucleusCount() const {return inNucleusCount;};
        
        //! Returns the number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
        /*!
            \return #intersectionCount The number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
         
            \sa distributeDose(const Track&)
         */
        virtual int getIntersectionCount() const {return intersectionCount;};
        
        //! Returns the number of domains composing the tMKM nucleus.
        /*!
            \return The number of domains composing the tMKM nucleus.
         */
        virtual int getNumberOfDomains() {return numberOfDomains;};
        
        //! Returns the nucleus position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the nucleus referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c y coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
         
            \sa Track_KieferChatterjee::getPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const;
        
        //! Returns the effective radius of the Nucleus_tMKM object.
        /*!
            \return #r_nucleus, the effective radius of the Nucleus_tMKM object. Since the structure of the final tMKM nucleus is "hexagon-like" this radius is different from the radius stored in the #cellLine reference. It's the distance between the center of the nucleus and the farthest point of the nucleus itself, expressed in um.
         */
        virtual double getRadius() const {return r_nucleus;};
        
        //! Save data corresponding to the dose absorbed by each domain and the number of lethal events observed, useful to debug.
        /*!
            \warning This method hasn't been implemented yet.
         
            \param fileName The name of the file where to save data.
         
            \warning The execution of the program will be terminated if fileName refers to an inexistent file.
         */
        void saveLocalDose(const std::string fileName) const;
        
    private:
        
        //! Create the domains as pointers to Nucleus_Integral_t objects, placed to form a hexagonal shape, spiraling from (0,0).
        /*!
            This function is called by the constructor every times a Nucleus_tMKM is created, and it's responsible to instantiate and place the domains in the right position.
         
            The structure is created in such a way that the center of each domain is placed on the vertex of a regular hexagon and the distance between two nearest neighbors is exactly equal to twice the radius of the domain (#domainRadius). Some concentric hexagons are created to places all the domains defined (#numberOfDomains). The center of each hexagon coincides with the position of the first domain created, that is also the center of the tMKM nucleus. This structure is created by means of the rotate() method.
         */
        void createDomains();
        
        //! Performs a 60 degrees clockwise rotation.
        /*!
            The rotation matrix is defined by:
            \f[
                \left( \begin{array}{cc}
                \cos\theta & \sin\theta \\
                -\sin\theta & \cos\theta \end{array} \right)
            \f]
         
            \param xTranslation The \c x coordinate of the point where to start the 60 degrees clockwise rotation.
            \param yTranslation The \c y coordinate of the point where to start the 60 degrees clockwise rotation.
         */
        void rotate(double &xTranslation,
                    double &yTranslation);
        
        //! A reference to a CellLine object where the characteristics of the cell line to which the tMKM nucleus belongs are stored.
        const CellLine &cellLine;
        
        //! The radius of the domain corresponding to the CellLine to which the tMKM nucleus belongs.
        /*!
            This information is stored in the #cellLine reference and then copied to this variable in the constructor. It's expressed in um.
         
            \sa Nucleus_tMKM()
         */
        double domainRadius;
        
        //! The number of domains composing the tMKM nucleus.
        /*!
            It's evaluated in the constructor as the ratio between the areas of the tMKM nucleus, whose radius \f$R_N\f$ is stored in the #cellLine reference, and the single domain, characterized by a radius \f$R_d\f$ (#domainRadius):
            \f[
                N_d=\frac{R_N^2}{R_d^2}
            \f]
         
            \sa Nucleus_tMKM()
         */
        int numberOfDomains;
        
        //! The position of the center of the nucleus (\c x coordinate) referred to the beam axis, expressed in mm.
        const double x_nucleus;
        
        //! The position of the center of the nucleus (\c y coordinate) referred to the beam axis, expressed in mm.
        const double y_nucleus;
        
        //! The linear-quadratic parameter \f$\alpha\f$ associated to each of the domains composing the tMKM nucleus.
        /*!
            It's instantiated in the constructor dividing the \f$\alpha\f$ parameter stored in the #cellLine reference by the total number of domains (#numberOfDomains): \f$\frac{\alpha}{N_d}\f$. Note that \f$\alpha\f$ is a parameter of the model, that ideally represents the value of the linear-quadratic \f$alpha_X\f$ parameter identified in the case of irradiation with a photon beam.
         
            \sa CellLine, CellLine::getParameters_LQ_noDt() and Nucleus_tMKM()
         */
        double alpha_d;
        
        //! The linear-quadratic parameter \f$\beta\f$ associated to each of the domains composing the tMKM nucleus.
        /*!
            It's instantiated in the constructor dividing the \f$\beta\f$ parameter stored in the #cellLine reference by the total number of domains (#numberOfDomains): \f$\frac{\beta}{N_d}\f$. Note that \f$\beta\f$ is a parameter of the model, that ideally represents the value of the linear-quadratic \f$beta_X\f$ parameter identified in the case of irradiation with a photon beam.
         
            \sa CellLine, CellLine::getParameters_LQ_noDt() and Nucleus_tMKM()
         */
        double beta_d;
        
        //! It's the effective radius of the Nucleus_tMKM object.
        /*!
            Since the structure of the final tMKM nucleus is "hexagon-like" this radius is different from the radius stored in the #cellLine reference. It's the distance between the center of the nucleus and the farthest point of the nucleus itself.
         
            It's is defined in the createDomains() function, called in the constructor, and it's expressed in um.
         
            \sa Nucleus_tMKM()
         */
        double r_nucleus;
        
        //! A pointer to a CellLine object, storing the information about the cell line to which the tMKM nucleus belongs.
        /*!
            It's defined in the createDomains() method and deleted in the destructor ~Nucleus_tMKM.
         */
        CellLine *domainCell;
        
        //! A pointer to pointers, where the objects finally pointed are Nucleus_Integral_t objects corresponding to the domains composing the tMKM nucleus.
        Nucleus_Integral_t* *domains;
        
        //! The number of times that the nucleus has been crossed through by a Particle.
        /*!
            It's incremented by means of the distributeDose(const Track&) method.
         */
        int inNucleusCount;
        
        //! The number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
        /*!
            It's incremented by means of the distributeDose(const Track&) method.
         */
        int intersectionCount;
    };
    
}


#endif /* NUCLEUS_TMKM_H */
