#ifndef CELL_LINE_H
#define CELL_LINE_H

#include <string>
#include <vector>

namespace Survival {

    //! Represents the cell line to which the nucleus belongs and hosts the radiobiological and structural characteristics of the cell and the different parametrizations of the models.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2011--2015
     
        Besides its hosting function, this class computes the local number of lethal events corresponding to a local dose deposition by means of the LEM I, LEM II, LEM III and LQ survival parametrizations for X-rays irradiation.
        Different parametrizations can be contemporary loaded in a CellLine object; the parametrization in use is specified via the public method  setParametrization(). The evaluation of the clustering damage enhancement of LEM II and LEM III can be performed in several ways: it can be fully generated on the fly via Monte Carlo, loaded from an external file or generated using an approximated analytical expression.
     */
    class CellLine
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            Each numeric data member is set to zero by default, with the exception of the parameters of the constructor. The selected parametrization is set to "noParametrization".
         
            \param cell_type A \c string identifying the name of the cell line.
            \param nucleus_radius The radius of the nucleus characteristic for the cell line, expressed in um (default 10 um).
            \param domain_radius The radius of the domain associated to the MKM parametrization of the nucleus, expressed in um (default 10 um).
         */
        CellLine(const std::string cell_type,
                 const double nucleus_radius = 10.0,
                 const double domain_radius = 10.0);
        
        //! Destructor.
        ~CellLine() {};
        
        //! Adds the LQ_noDt parametrization to the cell line (used in the MKM model) for the evaluation of the cellular survival.
        /*!
            Sets #alpha_X and #beta_X to the passed values.
         
            \param alphaX The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param betaX The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
         
            \sa parametrization_LQ_noDt()
         */
        void addParametrization_LQ_noDt(const double alphaX,
                                        const double betaX);
        
        //! Adds the LQ_noDt_T parametrization to the cell line (used in the MCt-MKM model) for the evaluation of the cellular survival.
        /*!
            Sets #alpha_X, #beta_X and #ac to the passed values.
         
            \param alphaX The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param betaX The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
            \param ac_ The time constant associated to the repair kinetics of the nucleus, expressed in \f$h^{-1}\f$.
         
            \sa parametrization_LQ_noDt_T()
         */
        void addParametrization_LQ_noDt_T(const double alphaX,
                                          const double betaX,
                                          const double ac_ = 2.187);
        
        //! Adds the LEM I parametrization of the survival to the cell line.
        /*!
            \param alphaX The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param betaX The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
            \param Dt The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
         
            \sa parametrization_LQ()
         */
        void addParametrization_LQ(const double alphaX,
                                   const double betaX,
                                   const double Dt);
        
        //! Adds the LEM II parametrization of the survival to the cell line.
        /*!
            \param alphaX The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param betaX The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
            \param Dt2 The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
            \param genome_Length The genome length expressed in unit of base pairs (#genomeLength).
            \param alphaSSB The number of SSBs directly produced by the irradiation in the whole genome per unit of dose absorbed (#alpha_SSB).
            \param alphaDSB The number of DSBs directly produced by the irradiation in the whole genome per unit of dose absorbed (#alpha_DSB).
            \param basePairs The distance (in number of based pairs) between two SSBs resulting in a DSB (#base_Pairs).
         
            \sa parametrization_LQ2()
         */
        void addParametrization_LQ2(const double alphaX,
                                    const double betaX,
                                    const double Dt2,
                                    const double genome_Length,
                                    const double alphaSSB = 1250.0,
                                    const double alphaDSB = 30.0,
                                    long int basePairs = 25);
        
        //! Adds the LEM III parametrization of the survival to the cell line.
        /*!
            \param alphaX The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param betaX The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
            \param Dt3 The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
            \param genome_Length The genome length expressed in unit of base pairs (#genomeLength).
            \param alphaSSB The number of SSBs directly produced by the irradiation in the whole genome per unit of dose absorbed (#alpha_SSB).
            \param alphaDSB The number of DSBs directly produced by the irradiation in the whole genome per unit of dose absorbed (#alpha_DSB).
            \param basePairs The distance (in number of based pairs) between two SSBs resulting in a DSB (#base_Pairs).
         
            \sa parametrization_LQ3()
         */
        void addParametrization_LQ3(const double alphaX,
                                    const double betaX,
                                    const double Dt3,
                                    const double genome_Length,
                                    const double alphaSSB = 1250.0,
                                    const double alphaDSB = 30.0,
                                    long int basePairs = 25);
        
        //! Evaluate the damage enhancement factor by means of an analytic approximated expression.
        /*!
            The damage enhancement factor is evaluated by means of the expression:
            \f[
                \left\langle\eta(D)\right\rangle=1+\sum_{n=1}^{h}\frac{\exp(-(\tilde{\alpha}_{SSB}+\tilde{\alpha}_{DSB})D)}{\tilde{\alpha}_{DSB}\,D}\frac{(\tilde{\alpha}_{SSB}\,D)^n(1-2^{(1-n)})}{n!}
            \f]
            derived from statistical consideration on the probability to observed near SSB on the DNA. In the formula, \em h represents #base_Pairs and \f$\tilde{\alpha}_{SSB}\f$ and \f$\tilde{\alpha}_{DSB}\f$ represent #alpha_SSB #alpha_DSB respectively multiplied by the ratio between #base_Pairs and #genomeLength.
         
            \note The analytic formula underestimates the real value of the damage enhancement factor. The best way to evaluate it is therefore the Monte Carlo simulation performed via the damageEnhancement() factor. The problem is the time necessary to the evaluation; hence for rapid estimates it could be used this method.
         */
        double analyticDamageEnhancement(const double dose) const;
        
        //! Evaluate the damage enhancement factor corresponding to a certain dose absorbed via a Monte Carlo simulation.
        /*!
            For the fixed dose:
                - It directly generates a number of DSB given by: \f$N_{DSB}(D)=\alpha_{DSB}\,L_{Genome}\,D\f$; placed in random position in the genome.
                - It directly generates a number of SSB given by: \f$N_{SSB}(D)=\alpha_{SSB}\,L_{Genome}\,D\f$; placed in random position on the two strands. SSB near to a DSB (in a window of width #base_Pairs centered on the DSB) are excluded from the computation.
                - The DNA is read base by base, when a SSB is identified, if it can be found another SSB inside a window of width #base_Pairs, a counter is incremented (\f$N_{2SSB}(D)\f$).
         
            The value of the resulting damage enhancement factor is evaluated by means of the relation:
            \f[
                \eta=\frac{N_{DSB}(D)+N_{2SSB}(D)}{N_{DSB}(D)}
            \f]
         
            \param dose The dose absorbed by the cell, expressed in Gy.
         
            \return The value of the damage enhancement factor correspondent to the dose absorbed.
         */
        double damageEnhancement(const double dose) const;
        
        //! Returns the time constant associated to the repair kinetics of the nucleus, expressed in \f$h^{-1}\f$.
        /*!
            \return The time constant associated to the repair kinetics of the nucleus, expressed in \f$h^{-1}\f$.
         */
        double getAC() const {return ac;};
        
        //! Returns a string identifying the name of the cell line.
        /*!
            \return #cellType The name of the cell line to which the nucleus belongs.
         */
        std::string getCellType() const {return cellType;};
        
        //! Returns the radius of the domain expressed in um.
        /*!
            \return #domainRadius The radius of the domain relative to the MKM parametrization of the nucleus, expressed in um.
         
            \sa Nucleus_MKM
         */
        double getDomainRadius() const {return domainRadius;};
        
        //! Returns the natural logarithm of the cellular survival evaluated on the basis of the selected parametrization.
        /*!
            The logarithmic survival is evaluated by calling the selected parametrization.
         
            \param dose The dose absorbed by the nucleus, needed to evaluate the cellular survival, expressed in Gy.
         
            \return The cellular survival associated to the dose absorbed.
         
            \sa parametrization_LQ(), parametrization_LQ_noDt() and CellLine::getLogSurvival_X(const vector<double>, const vector<double>)
         */
        double getLogSurvival_X(const double dose) const;
        
        //! Overload. Returns the natural logarithm of the cellular survival taking into account the time structure of the irradiation.
        /*!
            The logarithmic survival is evaluated by calling getParameters_LQ_noDt_T() if the corresponding parametrization is selected.
         
            \param doses A vector containing the sequence of doses deposited in the nucleus, expressed in Gy, each elements is associated to one interaction.
            \param times A vector containing the sequence of interaction times (expressed in hours), each elements is associated to one interaction.
         
            \return The cellular survival associated to the specific structure of doses absorbed.
         
            \warning The execution of the program will be terminated if the parametrization selected isn't "parametrization_LQ_noDt_T".
         
            \sa parametrization_LQ(), parametrization_LQ_noDt() and getLogSurvival_X()
         */
        double getLogSurvival_X(const std::vector<double>doses,
                                const std::vector<double>times) const;
        
        //! Returns the radius of the nucleus expressed in um.
        /*!
            \return #nucleusRadius The radius of the nucleus expressed in um.
         */
        double getNucleusRadius() const {return nucleusRadius;};
        
        //! Returns the linear quadratic \f$\alpha\f$ and \f$\beta\f$ parameters and the transition dose, corresponding to the selected parametrization, by overwriting three \c double variables passed by reference.
        /*!
            In the case of the "noDt-parametrization", D_t is set to -1.
         
            \param returnAlpha_X The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
            \param returnBeta_X The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
            \param returnD_t The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
         
            \sa setParametrization()
         */
        void getParameters(double &returnAlpha_X,
                           double &returnBeta_X,
                           double &returnD_t) const;
        
        //! Returns the linear quadratic \f$\alpha\f$ and \f$\beta\f$ parameters corresponding to the "LQ_noDt" parametrization by overwriting two \c double variables passed by reference.
        /*!
            \param returnAlpha_X The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param returnBeta_X The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if an incorrect parametrization is selected.
         
            \sa setParametrization() and parametrization_LQ_noDt()
         */
        void getParameters_LQ_noDt(double &returnAlpha_X,
                                   double &returnBeta_X) const;
        
        //! Returns the linear quadratic \f$\alpha\f$ and \f$\beta\f$ parameters corresponding to the "LQ_noDt" parametrization and the time constant by overwriting three \c double variables passed by reference.
        /*!
            \param returnAlpha_X The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param returnBeta_X The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
            \param ac_ The time constant associated to the repair kinetics of the nucleus, expressed in \f$h^{-1}\f$, passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if an incorrect parametrization is selected.
         
            \sa setParametrization() and parametrization_LQ_noDt_T()
         */
        void getParameters_LQ_noDt_T(double &returnAlpha_X,
                                     double &returnBeta_X,
                                     double &ac_) const;
        
        //void getParameters_LQ(double &returnAlpha_X1,
        //                      double &returnBeta_X1,
        //                      double &returnD_t) const;
        
        //! Returns the parameters characteristic of the LQ2 parametrization (used in the LEM II formulation) overwriting some variables passed by reference.
        /*!
            \param returnAlpha_X2 The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param returnBeta_X2 The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
            \param returnD_t2 The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy, passed by reference to be overwritten.
            \param returnGenomeLength The genome length expressed in unit of base pairs, passed by reference to be overwritten.
            \param returnAlpha_SSB The number of SSBs directly produced by the irradiation in the whole genome per unit of dose absorbed, passed by reference to be overwritten.
            \param returnAlpha_DSB The number of DSBs directly produced by the irradiation in the whole genome per unit of dose absorbed, passed by reference to be overwritten.
            \param returnBase_Pairs The distance (in number of based pairs) between two SSBs resulting in a DSB, passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if an incorrect parametrization is selected.
         
            \sa setParametrization, addParametrization_LQ2 and parametrization_LQ2
         */
        void getParameters_LQ2(double &returnAlpha_X2,
                               double &returnBeta_X2,
                               double &returnD_t2,
                               double &returnGenomeLength,
                               double &returnAlpha_SSB,
                               double &returnAlpha_DSB,
                               long int &returnBase_Pairs) const;
        
        //! Returns the parameters characteristic of the LQ3 parametrization (used in the LEM III formulation) overwriting some variables passed by reference.
        /*!
            \param returnAlpha_X3 The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param returnBeta_X3 The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
            \param returnD_t3 The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy, passed by reference to be overwritten.
            \param returnGenomeLength The genome length expressed in unit of base pairs, passed by reference to be overwritten.
            \param returnAlpha_SSB The number of SSBs directly produced by the irradiation in the whole genome per unit of dose absorbed, passed by reference to be overwritten.
            \param returnAlpha_DSB The number of DSBs directly produced by the irradiation in the whole genome per unit of dose absorbed, passed by reference to be overwritten.
            \param returnBase_Pairs The distance (in number of based pairs) between two SSBs resulting in a DSB, passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if an incorrect parametrization is selected.
         
            \sa setParametrization, addParametrization_LQ3 and parametrization_LQ3
         */
        void getParameters_LQ3(double &returnAlpha_X3,
                               double &returnBeta_X3,
                               double &returnD_t3,
                               double &returnGenomeLength,
                               double &returnAlpha_SSB,
                               double &returnAlpha_DSB,
                               long int &returnBase_Pairs) const;
        
        //! Get the value of the damage enhancement factor from the precalculated curve (#doseForEta, #etaPre) interpolating the nearest neighbors of the dose imposed.
        /*!
            \param dose The dose absorbed by the cell, expressed in Gy.
         
            \return The value of the damage enhancement factor correspondent to the dose absorbed.
         */
        double interpolatedDamageEnhancement(const double dose) const;
        
        //! If called it interrupts the execution of the program, because no parametrizazion is selected.
        double noParametrization(const double dummy) const;
        
        //! If called it interrupts the execution of the program, because no parametrizazion is selected.
        double noParametrization(const std::vector<double>v1,
                                 const std::vector<double>v2) const;
        
        //! Returns the logarithmic cellular survival corresponding to a particular dose absorbed. Implements the parametrization used in the MKM formulation.
        /*!
            The survival is evaluated by means of the standard linear quadratic relation:
            \f[
                S=\exp(-\alpha\,D-beta\,D^2)
            \f]
            
            \note The function returns the natural logarithm of the survival, NOT the survival.
        
            \param dose The dose absorbed expressed in Gy.
         
            \return The logarithmic cellular survival corresponding to a particular dose absorbed.
         
            \sa setParametrization(), parametrization_LQ_noDt_T() and parametrization_LQ()
         */
        double parametrization_LQ_noDt(const double dose) const;
        
        //! Returns the logarithmic cellular survival associated to a sequence of doses absorbed with a specific time structure. Implements the parametrization used in the MCt-MKM formulation.
        /*!
            At present, this parametrization is associated to the tMKM model (\ref tMKM "1").
            The logarithmic survival (L) is evaluated by means of the following relation:
            \f[
                L=-\alpha_d\left(\sum_{i=1}^N z_i\right) -\beta_d\left(\sum_{i=1}^N\right)^2 -2\beta\sum_{i=1}^{N-1}\sum_{j=i+1}^{N}\left[1-\exp\left(-(a+c)(t_j-t_i)\right)\right]z_i\,z_j
            \f]
            where N represent the length of the vector of doses (or times), the sum (a+c) represents #ac and \f$z_i\f$ and \f$t_i\f$ represent the i-th element of the vectors of times and doses absorbed respectively.
         
            \param doses The vector representing the sequence of doses absorbed, expressed in Gy.
            \param times The vector representing the sequence of interaction times, expressed in hours.
        
            \return The logarithmic cellular survival associated to a sequence of doses absorbed with a specific time structure.
         
            \sa setParametrization(), parametrization_LQ_noDt_T(), parametrization_LQ() and Nucleus_tMKM
         
            \anchor tMKM 1. L. Manganaro, ..., A. Attili, ...
         */
        double parametrization_LQ_noDt_T(const std::vector<double> doses,
                                         const std::vector<double> times) const;
        
        //! Returns the logarithmic cellular survival corresponding to a particular dose absorbed. Implements the parametrization used in the LEM I formulation.
        /*!
            It implements the LEM I parametrization of the cellular survival based on the assumption that beyond a dose \f$D_t\f$ the standard parametrization of the LQ model is no more valid as it was observed "an exponential tail" (function of the dose absorbed).
            \f[
                S=\exp(-\alpha\,D-beta\,D^2)\qquad D<D_t
            \f]
            \f[
                S=\exp(-\alpha\,D_t-beta\,D_t^2)\exp(-s(D-D_t))=S_t\exp(-s(D-D_t))\qquad D>=D_t
            \f]
            where \f$s=\alpha+2\,\beta D_t\f$.
         
            \param dose The dose absorbed expressed in Gy.
         
            \return The logarithmic cellular survival corresponding to a particular dose absorbed.
         
            \sa setParametrization(), parametrization_LQ_noDt_T() and parametrization_LQ()
         */
        double parametrization_LQ(const double dose) const;
        
        //! Returns the logarithmic cellular survival corresponding to a particular dose absorbed. Implements the parametrization used in the LEM II formulation.
        /*!
            The LEM II formulation adopts the same parametrization used in the LEM I (see parametrization_LQ()) but, to accounts for the so called \em clustering \em effect, an enhancement factor \f$\eta\f$ is added which multiplies the dose absorbed when \f$D>D_t\f$.
            The resulting parametrization can be written as:
            \f[
                S=\exp(-\alpha\,D-beta\,D^2)\qquad D<D_t
            \f]
            \f[
                S=\exp(-\alpha\,D_t-beta\,D_t^2)\exp\left[-s(\eta (D)D-D_t)\right]=S_t\exp\left[-s(\eta (D)D-D_t)\right]\qquad D>=D_t
            \f]
            where \f$s=\alpha+2\,\beta D_t\f$.
         
            \f$\eta\f$ is a function of the dose absorbed, and there are several ways to generate it:
                - It can be generated via a Monte Carlo Simulation (damageEnhancement())
                - It can be generated via an analytic approximation (analyticDamageEnhancement())
                - It can be read from an external file (readDamageEnhancement())
         
            \param dose The dose absorbed expressed in Gy.
         
            \return The logarithmic cellular survival corresponding to a particular dose absorbed.
         
            \sa setParametrization(), parametrization_LQ_noDt_T() and parametrization_LQ()
         */
        double parametrization_LQ2(const double dose) const;
        
        //! Returns the logarithmic cellular survival corresponding to a particular dose absorbed. Implements the parametrization used in the LEM II formulation.
        /*!
            The parametrization used is identical to the one defined in the LEM II formulation (see parametrization_LQ2()). Briefly, the survival is evaluated as:
            \f[
                S=\exp(-\alpha\,D-beta\,D^2)\qquad D<D_t
            \f]
            \f[
                S=\exp(-\alpha\,D_t-beta\,D_t^2)\exp\left[-s(\eta (D)D-D_t)\right]=S_t\exp\left[-s(\eta (D)D-D_t)\right]\qquad D>=D_t
            \f]
         
            \param dose The dose absorbed expressed in Gy.
         
            \return The logarithmic cellular survival corresponding to a particular dose absorbed.
         
            \sa setParametrization(), parametrization_LQ_noDt_T() and parametrization_LQ()
         */
        double parametrization_LQ3(const double dose) const;
        
        //! Reads and returns the value of the damage enhancement factor correspondent to the required dose from an external file.
        /*!
            \param dose The dose absorbed, expressed in Gy, to calculate the correspondent damaga enhancement factor.
         
            \return The value of the damage enhancement factor.
         
            \note This method was thought to import the exact curve published in the LEM description.
         
            \warning The execution of the program will be terminated if the file doesn't exist.
         */
        double readDamageEnhancement(const double dose) const;
        
        //! Sets the radius of the domain relative to the MKM parametrization of the nucleus.
        /*!
            \param domainRadius_ The radius of the domain expressed in um.
         
            \sa #domainRadius, Nucleus_MKM
         */
        void setDomainRadius(double domainRadius_) {domainRadius = domainRadius_;};
        
        //! Sets the radius of the nucleus.
        /*!
            \param nucleusRadius_ The radius of the nucleus expressed in um.
         
            \sa #nucleusRadius
         */
        void setNucleusRadius(double nucleusRadius_) {nucleusRadius = nucleusRadius_;};
        
        //! Sets an X-ray parametrization for the evaluation of the cellular survival.
        /*!
            \param parametrization_type A string indicating the name of the desired parametrization.
         
            The possible choices for the parametrization and the correspondents options set by the function are listed in the following table:
         
            | Parametrization           | Model   | selectedParametrization     | selectedDamageEnhancement       | selectedEtaGeneration       |
            | :-----------------------: | :-----: | :-------------------------: | :-----------------------------: | :-------------------------: |
            | LQ                        | LEM I   | parametrization_LQ()        | None                            | None                        |
            | LQ2                       | LEM II  | parametrization_LQ2()       | interpolatedDamageEnhancement() | readDamageEnhancement()     |
            | LQ3                       | LEM III | parametrization_LQ3()       | interpolatedDamageEnhancement() | readDamageEnhancement()     |
            | LQ_noDt                   | MKM     | parametrization_LQ_noDt()   | None                            | None                        |
            | LQ_noDt_T                 | MCt-MKM | parametrization_LQ_noDt_T() | None                            | None                        |
            | LQ2_readfile              | LEM II  | parametrization_LQ2()       | interpolatedDamageEnhancement() | readDamageEnhancement()     |
            | LQ2_interpolated_MC       | LEM II  | parametrization_LQ2()       | interpolatedDamageEnhancement() | damageEnhancement()         |
            | LQ2_interpolated_analytic | LEM II  | parametrization_LQ2()       | interpolatedDamageEnhancement() | analyticDamageEnhancement() |
            | LQ2_punctual_analytic     | LEM II  | parametrization_LQ2()       | analyticDamageEnhancement()     | analyticDamageEnhancement() |
            | LQ2_punctual_MC           | LEM II  | parametrization_LQ2()       | damageEnhancement()             | damageEnhancement()         |
         
            If the parametrization selected needs an eta generated, then the function generates it by (recursively) using the selectedEtaGeneration pointer, storing the calculated values of dose and \f$\eta\f$ in #doseForEta and #etaPre respectively.
         
            \note \c #selectedParametrization, \c #selectedDamageEnhancement and \c #selectedEtaGeneration are pointers to functions.
         
            \warning The execution of the program will be terminated if an inexistent parametrization is selected.
         */
        void setParametrization(const std::string parametrization_type);
        
        
    private:
        
        //! A \c string identifying the name of the cell line.
        std::string cellType;
        
        //! The radius of the nucleus characteristic for the cell line, expressed in um.
        double nucleusRadius;
        
        //! The radius of the domain associated to the MKM parametrization of the nucleus, expressed in um.
        double domainRadius;
        
        
        //! The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
        double alpha_X;
        
        //! The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
        double beta_X;
        
        //! A boolean value identifying if the "LQ_noDt" parametrization is selected.
        bool isLQ_noDtLoaded;
        
        
        //! The time constant associated to the repair kinetics of the nucleus, expressed in \f$h^{-1}\f$.
        double ac;
        
        //! A boolean value identifying if the "LQ_noDt_T" parametrization is selected.
        bool isLQ_noDt_TLoaded;
        
        
        //! The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
        double alpha_X1;
        
        //! The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
        double beta_X1;
        
        //! The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
        double D_t;
        
        //! The coefficient of the exponential tail: \f$s=\alpha+2\beta D_t\f$
        double s;
        
        //! The logarithmic survival associated to a dose absorbed #D_t, evaluated according to the standard linear quadratic parametrization.
        double logS_t;
        
        //! A boolean value identifying if the "LQ" parametrization is selected.
        bool isLQloaded;
        
        
        //! The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
        double alpha_X2;
        
        //! The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
        double beta_X2;
        
        //! The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
        double D_t2;
        
        //! The coefficient of the exponential tail: \f$s=\alpha+2\beta D_t\f$
        double s2;
        
        //! The logarithmic survival associated to a dose absorbed #D_t, evaluated according to the standard linear quadratic parametrization.
        double logS_t2;
        
        //! The genome length expressed in number of base pairs.
        double genomeLength;
        
        //! The number of SSBs directly produced by the irradiation in the whole genome per unit of dose absorbed.
        double alpha_SSB;
        
        //! The number of DSBs directly produced by the irradiation in the whole genome per unit of dose absorbed.
        double alpha_DSB;
        
        //! The distance (in unit of based pairs) between two SSBs resulting in a DSB.
        long int base_Pairs;
        
        //! A boolean value identifying if one the "LQ2" parametrizations is selected.
        /*!
            Possible cases are:
                - "LQ2"
                - "LQ2_interpolated_analytic"
                - "LQ2_interpolated_MC"
                - "LQ2_interpolated_readfile"
                - "LQ2_punctual_analytic"
                - "LQ2_punctual_MC"
         */
        bool isLQ2loaded;
        
        
        //! The linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$.
        double alpha_X3;
        
        //! The linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$.
        double beta_X3;
        
        //! The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy.
        double D_t3;
        
        //! The coefficient of the exponential tail: \f$s=\alpha+2\beta D_t\f$
        double s3;
        
        //! The logarithmic survival associated to a dose absorbed #D_t, evaluated according to the standard linear quadratic parametrization.
        double logS_t3;
        
        //! A boolean value identifying if the "LQ3" parametrization is selected.
        bool isLQ3loaded;
        
        
        //! A pointer to functions that identifies the selected way to evaluate the enhancement factor \f$\eta (D)\f$ in LEM II and III formulations.
        /*!
            Possible choices are:
                - interpolatedDamageEnhancement()
                - analyticDamageEnhancement()
                - damageEnhancement()
         
            It's used in parametrization_LQ2() and parametrization_LQ3() methods.
         */
        double (CellLine::*selectedDamageEnhancement)(const double dose) const;
        
        //! A pointer to functions that identifies the selected way to evaluate the enhancement factor \f$\eta (D)\f$ in LEM II and III formulations.
        /*!
            Possible choices are:
                - readDamageEnhancement()
                - analyticDamageEnhancement()
                - damageEnhancement()
         
            It's used in setParametrization() to calculate a curve for \f$\eta\f$ as a function of the dose absorbed to be stored in #doseForEta and #etaPre arrays.
         
            \sa parametrization_LQ2() and parametrization_LQ3()
         */
        double (CellLine::*selectedEtaGeneration)(const double dose) const;
        
        //! A pointer to functions that identifies the selected parametrization.
        /*!
            Possible parametrizations are:
                - noParametrization()
                - parametrization_LQ_noDt() (MKM)
                - parametrization_LQ() (LEM I)
                - parametrization_LQ2() (LEM II)
                - parametrization_LQ3() (LEM III)
         */
        double (CellLine::*selectedParametrization)(const double dose) const;
        
        //! A pointer to functions that identifies the selected parametrization when the temporal effects of the irradiation are taken into account.
        /*!
            Possible parametrizations are:
                - noParametrization(const std::vector<double>, const std::vector<double>)
                - parametrization_LQ_noDt_T() (MCt-MKM)
         */
        double (CellLine::*selectedParametrizationT)(const std::vector<double>doses,
                                                     const std::vector<double>times) const;
        
        
        //! An array used to store the values of dose for to precalculate the enhancement factor curve.
        /*!
            It is constituted by 200 values logarithmically spaced in \f$[100,5\cdot 10^6]\f$.
         */
        double doseForEta[200];
        
        //! An array containing the precalculated values of the enhancement factor as a function of the dose absorbed.
        /*!
            It is constituted by 200 values indicated the enhancement factor for each value of the #doseForEta array.
         */
        double etaPre[200];
        
        //! A boolean data member that indicates if the selected parametrization requires the generation of the enhancement factor.
        static bool needEtaGenerated;
        
        //! An array representing the genome of the cell.
        static int DNA[10000000];
        
        //! A boolean array representing the single strand breaks (SSB) in the first strand.
        static bool SSB1[10000000];
        
        //! A boolean array representing the single strand breaks (SSB) in the second strand.
        static bool SSB2[10000000];
        
        //! A boolean array representing the double strand breaks (DSB) in the genome.
        static bool DSB[10000000];
    };
    
}


#endif /* CELLLINE_H */
