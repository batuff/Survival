#ifndef NUCLEUS_MKM_H
#define NUCLEUS_MKM_H

#include "Nucleus.h"

namespace Survival {
    
    class Nucleus_Pixel;
    class Nucleus_Integral;
    
    //! Inherited from the Nucleus pure virtual class, it implements the cellular nucleus as defined and used in the MKM model.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2011--2015
     
        The microdosimetric-kinetic model (MKM, \ref MKM "2") is based on a cellular nucleus divided into subcellular structures referred to as \a domains similar to the \a sites defined in the theory of dual radiation action (TDRA) by Kellerer and Rossi (\ref TDRA "1"). This class implements the nucleus structure characteristics of the MKM model in all his formulations and provides some methods to evaluate the dose absorbed in the interaction with particles and get informations about lethal events observed and the associated cellular survival.
     
        \anchor TDRA 1. A. Kellerer and H. Rossi, "A generalized formulation of dual radiation action", \a Radiation \a Research \b 75, 471-488 (1978)
     
        \anchor MKM 2. The MKM model was formulated by Hawkins in 1994, then modified over subsequent years and recently reformulated (and here implemented) following a Monte Carlo approach. The original published reference for the MKM is:
            - R.B. Hawkins, "A Statistical Theory of Cell Killing by Radiation of Varying Linear Energy Transfer", \a Radiation \a Research \b 140, 366-374 (1994).
        While for the recent Monte Carlo reformulation:
            - L. Manganaro, G. Russo, A. Attili, "Advanced modeling of the temporal effect in particle therapy: from radiobiological evaluation to treatment planning", \a Medical \a Physics, (Submitted)
     */
    class Nucleus_MKM : public Nucleus
    {
        // tell the compiler we are implicitly overriding and then overloading -- suppress [-Woverloaded-virtual] warnings
        using Survival::Nucleus::addNucleusDoses;
        using Survival::Nucleus::getDosesAndLethals;
        
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            When the constructor is called, it instantiates the object creating a hexagonal structure of circular domains, using the informations stored in the CellLine reference object, such as the radius of the nucleus and the radius of the single domain. This is made by calling the createDomains() function.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa createDomains(), cleanNucleus(), Nucleus_MKM(const CellLine&, double, int, const double, const double) and Nucleus_tMKM
         */
        Nucleus_MKM(const CellLine &cellLineRef,
                    const double xPosition = 0.0,
                    const double yPosition = 0.0);
        
        //! Instantiates and sets the object. Overload of the constructor.
        /*!
            Provide the possibility to specify also the total number of domains and the radius of each domain.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param domainRadius The radius of the single domain in the MKM nucleus, expressed in um.
            \param numberOfDomains The number of domains that constitute the MKM nucleus.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa createDomains(), Nucleus_Integral and Nucleus_MKM(const CellLine&, const double, const double)
         */
        Nucleus_MKM(const CellLine &cellLineRef,
                    double domainRadius,
                    int numberOfDomains,
                    const double xPosition = 0.0,
                    const double yPosition = 0.0);
        
        //! Destructor.
        /*!
            The destructor deletes the domainCell, domains and each domains[i-th] pointers created when the object is instantiated.
         
            \sa createDomains()
         */
        virtual ~Nucleus_MKM();
        
        //! Adds a constant value of dose absorbed in each domain of the nucleus.
        /*!
            The method calls systematically for each domain the function Nucleus_Integral::addBackgroundDose(), that is the override of this method itself in the Nucleus_Integral class. The result is to add a constant value of dose absorbed in each domain and consequently in the whole nucleus.
         
            \param dose The dose to be added expressed in Gy.
         */
        void addBackgroundDose(const double dose);
        
        //! Performs the sum of the dose absorbed by two different nucleus domain to domain (adds one to the other).
        /*!
            This method calls systematically for each domain the function Nucleus_Integral::getDoseAndLethalForDomain() to get the respectively domain dose and then adds it to the corresponding domain in the current nucleus via the Nucleus_Integral::addBackgroundDose() function.
         
            \warning The execution of the program will be terminated if the other nucleus is not geometrically similar to this one (e.g. different #r_nucleus, #domainRadius or #numberOfDomains).
         
            \param nucleus A reference to another Nucleus_MKM object to evaluate the sum of the doses.
         */
        virtual void addNucleusDoses(Nucleus_MKM &nucleus);
        
        //! Resets to zero all counters (#inNucleusCount and #intersectionCount) and doses in the nucleus and his domain.
        /*!
            Resets to zero #inNucleusCount and #intersectionCount and calls systematically for each domain the function Nucleus_Integral::cleanNucleus(), that this the override of this method itself in the Nucleus_Integral class.
         */
        virtual void cleanNucleus();
        
        //! Returns a pointer to a new Nucleus_MKM object. It not really a clone but a new clean object.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         
            \note To create a real clone of another nucleus, a better implementation of the copy constructor is needed.
         */
        virtual Nucleus_MKM* clone(const CellLine&);
        
        //! When an interaction between a Particle and the MKM nucleus occurs, the method increases #inNucleusCount and #intersectionCount counters and proceeds to distribute the dose in each domain.
        /*!
            This function is the first step to evaluate the dose deposited by the radiation in the nucleus or, better, in each domain of the MKM nucleus. It checks if any interaction exists between Nucleus and Particle (that is the Track generated by the particle) looking at their positions and radius. If it's true, it increases the respective counters (#inNucleusCount and #intersectionCount) and calls the method Nucleus_Integral::distributeDose() in a \c for loop over the total number of domains.
         
            \param track A reference to the Track generated by the particle in the nucleus.
         */
        virtual void distributeDose(const Track &track);
        
        //! Overload of distributeDose(const Track &track) to manage a Tracks object, it simply calls distributeDose(const Track &track) for each track of the container.
        /*!
            Since the Tracks class is a container for Track objects, this method calls distributeDose(const Track &track) in a \c for loop over each track contained in the Tracks object.
         */
        virtual void distributeDose(const Tracks &tracks);
        
        //! Returns the name of the cell line to which the nucleus belongs.
        /*!
            \return A \c string corresponding to the name of the cell line to which the nucleus belongs, getting the information by #cellLine.
         */
        virtual std::string getCellType() const;
        
        //! Returns the radius of the domain corresponding to the CellLine to which the MKM nucleus belongs.
        /*!
            \return #domainRadius That is the radius of the domain corresponding to the CellLine to which the MKM nucleus belongs.
         */
        virtual double getDomainRadius() {return domainRadius;};
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed and the number of lethal events observed in a specified domain identified by domainIndex, with respective uncertainties.
        /*!
            The function calls Nucleus_Integral::getDoseAndSurvival() to get the dose absorbed by the domainIndex-th domain and then evaluates the number of lethal events by means of the linear-quadratic relation:
            \f[
                L=\alpha_d\,D+\beta_d\,D^2
            \f]
         
            \param domainIndex The index referred to the domain, passed by reference to be overwritten.
            \param dose The dose absorbed by the domainIndex-th domain, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param lethal The lethal events observed in the domainIndex-th domain, passed by reference to be overwritten.
            \param lethalUncertainty The uncertainty associated to the lethal events observed, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and lethals, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \warning If an incorrect domainIndex is chosen (e.g. greater than the real number of domains) -1 will be assigned also to dose and lethals, no exceptions are thrown.
         
            \sa getDosesAndLethals(), getDoseAndSurvival() and getDoseForDomain()
         */
        void getDoseAndLethalForDomain(int domainIndex,
                                       double &dose,
                                       double &doseUncertainty,
                                       double &lethal,
                                       double &lethalUncertainty) const;
        
        //! Returns (overwriting parameters passed by reference) the vectors containing doses absorbed and number of lethal events observed with respective uncertainties, each element of the vector refers to one of the domains.
        /*!
            The function calls, in a \c for loop over the total number of domains, the method getDoseAndLethalForDomain() to get the informations desired.
         
            \param doses The vector of doses absorbed (in Gy), each element refers to a specific domain, passed by reference to be overwritten.
            \param dosesUncertainty The vector of uncertainties associated to doses absorbed (in Gy), each element refers to a specific domain, passed by reference to be overwritten.
            \param lethals The vector of lethal events observed, each element refers to a specific domain, passed by reference to be overwritten.
            \param lethalsUncertainty The vector of uncertainties associated to the number of lethal events observed, each element refers to a specific domain, passed by reference to be overwritten.
         
            \sa getDoseAndSurvival() and getDoseForDomain()
         */
        virtual void getDosesAndLethals(std::vector<double> &doses,
                                        std::vector<double> &dosesUncertainty,
                                        std::vector<double> &lethals,
                                        std::vector<double> &lethalsUncertainty) const;
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and the associated survival, with respective uncertainties.
        /*!
            The function calls, in a \c for loop over the total number of domains, the method Nucleus_Integral::getDoseAndSurvival() to get the dose absorbed by each domain and for each one evaluates the number of lethal events by means of the linear-quadratic relation:
            \f[
                L_i=\alpha_d\,D+\beta_d\,D^2
            \f]
            (An equivalent solution could be reached calling directly getDoseAndLethalForDomain()).
            Then the total dose absorbed by the nucleus is obtained by the average of the doses absorbed by domains while the cellular survival is evaluated as a negative exponential function of the total number of lethal events (according to the poissonian statistics):
            \f[
                S=\exp(-L_{TOT})
            \f]
            where \f$L_{TOT}=\sum L_i\f$
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival associated to the dose absorbed by the nucleus, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and survival, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \sa getDosesAndLethals(), getDoseAndLethalForDomain() and getDoseForDomain()
         */
        virtual void getDoseAndSurvival(double &dose,
                                        double &doseUncertainty,
                                        double &survival,
                                        double &survivalUncertainty) const;
        
        //! Return the dose absorbed in the indexOfDomain-th domain, expressed in Gy.
        /*!
            The function calls the method Nucleus_Integral::getDoseAndSurvival() from the indexOfDomain-th domain.
         
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
        
        //! Returns the number of domains composing the MKM nucleus.
        /*!
            \return The number of domains composing the MKM nucleus.
         */
        virtual int getNumberOfDomains() {return numberOfDomains;};
        
        //! Returns the nucleus position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the nucleus referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c y coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
         
            \sa Nucleus_Integral::getPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const;
        
        //! Returns the effective radius of the Nucleus_MKM object.
        /*!
            \return #r_nucleus, the effective radius of the Nucleus_MKM object. Since the structure of the final MKM nucleus is "hexagon-like" this radius is different from the radius stored in the #cellLine reference. It's the distance between the center of the nucleus and the farthest point of the nucleus itself, expressed in um.
         */
        virtual double getRadius() const {return r_nucleus;};
        
        //! Save data corresponding to the dose absorbed by each domain and the number of lethal events observed, useful to debug.
        /*!
            This method write on a file, for each domain of the nucleus:
                - The index of the domain;
                - The radius of the domain;
                - The \c x and \c y coordinates (in um) identifying its position referred to the beam axis;
                - The dose absorbed by the domain;
                - The number of lethal events observed.
        
            \param fileName The name of the file where to save data.
         
            \warning The execution of the program will be terminated if fileName refers to an inexistent file.
         */
        void saveLocalDose(const std::string fileName) const;
        
    private:
        
        //! Create the domains as pointers to Nucleus_Integral objects, placed to form a hexagonal shape, spiraling from (0,0).
        /*!
            This function is called by the constructor every times a Nucleus_MKM is created, and it's responsible to instantiate and place the domains in the right position.
         
            The structure is created in such a way that the center of each domain is placed on the vertex of a regular hexagon and the distance between two nearest neighbors is exactly equal to twice the radius of the domain (#domainRadius). Some concentric hexagons are created to places all the domains defined (#numberOfDomains). The center of each hexagon coincides with the position of the first domain created, that is also the center of the MKM nucleus. This structure is created by means of the rotate() method.
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
        
        //! A reference to a CellLine object where the characteristics of the cell line to which the MKM nucleus belongs are stored.
        const CellLine &cellLine;
        
        //! The radius of the domain corresponding to the CellLine to which the MKM nucleus belongs.
        /*!
            This information is stored in the #cellLine reference and then copied to this variable in the constructor. It's expressed in um.
         
            \sa Nucleus_MKM()
         */
        double domainRadius;
        
        //! The number of domains composing the MKM nucleus.
        /*!
            It's evaluated in the constructor as the ratio between the areas of the MKM nucleus, whose radius \f$R_N\f$ is stored in the #cellLine reference, and the single domain, characterized by a radius \f$R_d\f$ (#domainRadius):
            \f[
                N_d=\frac{R_N^2}{R_d^2}
            \f]
         
            \sa Nucleus_MKM()
         */
        int numberOfDomains;
        
        //! The position of the center of the nucleus (\c x coordinate) referred to the beam axis, expressed in mm.
        const double x_nucleus;
        
        //! The position of the center of the nucleus (\c y coordinate) referred to the beam axis, expressed in mm.
        const double y_nucleus;
        
        //! The linear-quadratic parameter \f$\alpha\f$ associated to each of the domains composing the MKM nucleus.
        /*!
            It's instantiated in the constructor dividing the \f$\alpha\f$ parameter stored in the #cellLine reference by the total number of domains (#numberOfDomains): \f$\frac{\alpha}{N_d}\f$. Note that \f$\alpha\f$ is a parameter of the model, that ideally represents the value of the linear-quadratic \f$alpha_X\f$ parameter identified in the case of irradiation with a photon beam.
         
            \sa CellLine, CellLine::getParameters_LQ_noDt() and Nucleus_MKM()
         */
        double alpha_d;
        
        //! The linear-quadratic parameter \f$\beta\f$ associated to each of the domains composing the MKM nucleus.
        /*!
            It's instantiated in the constructor dividing the \f$\beta\f$ parameter stored in the #cellLine reference by the total number of domains (#numberOfDomains): \f$\frac{\beta}{N_d}\f$. Note that \f$\beta\f$ is a parameter of the model, that ideally represents the value of the linear-quadratic \f$beta_X\f$ parameter identified in the case of irradiation with a photon beam.
         
            \sa CellLine, CellLine::getParameters_LQ_noDt() and Nucleus_MKM()
         */
        double beta_d;
        
        //! It's the effective radius of the Nucleus_MKM object.
        /*!
            Since the structure of the final MKM nucleus is "hexagon-like" this radius is different from the radius stored in the #cellLine reference. It's the distance between the center of the nucleus and the farthest point of the nucleus itself.
         
            It's is defined in the createDomains() function, called in the constructor, and it's expressed in um.
         
            \sa Nucleus_MKM()
         */
        double r_nucleus;
        
        //! A pointer to a CellLine object, storing the information about the cell line to which the MKM nucleus belongs.
        /*!
            It's defined in the createDomains() method and deleted in the destructor ~Nucleus_MKM.
         */
        CellLine *domainCell;
        
        //! A pointer to pointers, where the objects finally pointed are Nucleus_Integral objects corresponding to the domains composing the MKM nucleus.
        Nucleus_Integral* *domains;
        
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


#endif /* NUCLEUS_MKM_H */
