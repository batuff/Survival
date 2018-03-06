#ifndef CALCULUS_H
#define CALCULUS_H

#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

namespace Survival {

    class Nucleus;
    class CellLine;
    class Tracks;
    
    //! It implements some methods to simulate the irradiation process evaluating the cellular survival on a population of cells and extrapolating the radiobiological LQ parameters.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2011--2018
     
        The class makes use of Track and Nucleus objects to simulate the irradiation process. Different methods are implemented in order to perform Monte Carlo simulations or approximated analytic evaluations of the process.
     */
    class Calculus
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            It needs references to Tracks, CellLine and Nucleus objects to be used in the simulation. It needs also the prefix of the output file and provides the possibility to fix the seed of the pseudorandom numbers generator and to manage the number of threads (only necessary for methods that implement a multithread system).
            The class makes use of the GSL library to manage the generation of pseudorandom numbers, hence an instance of the Tausworthe generator is created and the pseudorandom numbers generator is seeded.
         
            \param tracksRef A reference to a Tracks object to manage the radiation.
            \param cellLineRef A reference to a CellLine object containing the structural and radiobiological informations.
            \param nucleusRef A reference to a Nucleus object to be irradiated.
            \param savePrefix A \c string indicating the prefix of the output file.
            \param model The model used in the simulation.
            \param p_type The type of parallelism desired (see #nThreads).
            \param randomSeed The seed to be used in the initialization of the pseudorandom numbers generator.
         */
        Calculus(const Tracks &tracksRef,
                 const CellLine &cellLineRef,
                 Nucleus &nucleusRef,
                 std::string savePrefix,
                 std::string model,
                 int p_type = 1,
                 long randomSeed = 0);
        
        //! Destructor.
        /*!
            Deletes the Calculus object and the gsl_rng pointer.
         */
        ~Calculus();
        
        //! Returns the number of threads used in the simulation.
        /*!
            \return The number of threads used in the simulation.
         
            \sa #nThreads
         */
        int getNThreads(void) const {return nThreads;};
        
        //! Generate an ordered fixed-size sequence of pseudorandom numbers in [begin, begin+duration].
        /*!
            \param nEv Size of the sequence, that is the number of pseudorandom numbers generated.
            \param begin The minimum number extractable.
            \param duration The width of the interval.
         
            \return A STL vector containing an ordered sequence of pseudorandom numbers.
         
            \warning The execution of the program will be terminated if a negative "duration" is chosen.
         */
        std::vector<double> generateSequence(int nEv,
                                             double begin,
                                             double duration);
        
        //! Perform a Monte Carlo simulation delivering a fixed total dose with a defined time structure to a cell population, evaluating the resulting cell survival in order to get the value of the G factor that defines the temporal effect of the irradiation.
        /*!
            The user has to fix:
            1. The precision required for the simulation, that is the statistics to be reached to terminate the simulation. Two possibilities are supported, so the user can:
                - Fix the number of iterations, hence the \c precision has to be an integer value greater (or at least equal) to 1.
                - Define a constraint on the precision to reach in the simulation in the evaluation of the cell survival (precisely the relative error on the survival), hence the \c precision has to be a \c double in (0, 1).
            2. The structure of the irradiation, which is defined by:
                - The total dose to be delivered.
                - The number of fractions of the treatment.
                - The time spacing between two consecutive fractions.
                - The time needed to deliver a single fraction.
            Once set these features, the irradiation process of the cell population is simulated by calling the histogram_dose_survival_t() method in a parallel loop. Each thread performs the irradiation of a single cell (with its complete time structure) and at the end of each irradiation the method updates the total mean dose and survival observed, with associated uncertainties, and determines if the precision set by the user is reached or not. All the values of dose and survival obtained, cell by cell, are saved in an output file.
            Once reached the precision required, the method evaluates the G factor using the LQ \f$\alpha\f$ and \f$\beta\f$ parameters obtained in the case of acute irradiation (parameters of the function) by means of the relation:
            \f[
                G=\frac{-\ln(S)-\alpha D}{\beta D^2}
            \f]
            where D is the total dose delivered and S the resulting survival.
            It evaluates also the G factor following the analytical expression (valid only in case of fractionated treatment with neglecting the intra-fraction time structure):
            \f[
                G=1-\frac{2}{D^2}\sum_{i=0}^{N-1}\sum_{j=i+1}^{N-2}\left(1-\exp[-(a+c)(t_j-t_i)]\right)d^2
            \f]
            where d represents the dose delivered per fraction, N the total number of fractions, \f$(t_j-t_i)\f$ represents the temporal distance between two fractions and (a+c) is the time constant (see Nucleus_tMKM::ac).
            The calculated values, together with the informations on the time structure are saved in another output file.
         
            \note Survival uncertainty: the estimation of the survival uncertainty is knowingly biased, but extremely efficient and memory-friendly. Note that this quantity is used only to decide when it is possible to stop the simulation hence it doesn't affect the result (and that the bias decrease with increasing iteration).
         
            \param trackMode Defined for completeness. Only "histogram" is supported.
            \param totalDose The total dose to be delivered, expressed in Gy.
            \param nFrac The total number of fractions constituting the treatment.
            \param timeSpacing The time spacing between two consecutive fractions, expressed in hours.
            \param fracDeliveryTime The time needed to deliver a single fraction, expressed in hours.
            \param precision The precision to be reached in the simulation (could be a fixed number of iterations or a constraint on the survival precision).
            \param alpha LQ \f$\alpha\f$ parameter obtained in acute conditions.
            \param beta LQ \f$\beta\f$ parameter obtained in acute conditions.
         
            \warning The execution of the program will be terminated if the precision is not set correctly.
         
            \sa histogram_dose_survival_t() and slow_meanDose_meanSurvival()
         */
        void evaluateG(const std::string trackMode,
                       double totalDose,
                       int nFrac,
                       double timeSpacing,
                       double fracDeliveryTime,
                       double precision,
                       double alpha,
                       double beta);
        
        //! Returns the model used in the simulation.
        /*!
            \return A string indicating the model used in the simulation.
         */
        std::string getModel() const {return model;};
        
        //! Simulates, via the Monte Carlo method, the irradiation process of a nucleus which is part of the cellular population.
        /*!
            First of all the function calculate the fluence of the beam on the base of the dose imposed, by means of the following relation:
            \f[
                \Phi=\rho\frac{D}{K\cdot \langle LET\rangle}
            \f]
            where D is the dose imposed, K is a constant to accounts for the different units of measure used, \f$\rho\f$ is the density of the target and \f$\langle LET\rangle\f$ is the mean value of LET of the radiation (that coincides with the LET of the particle only in the monoenergetic case).
            The number of particles interacting with the nucleus is randomly extracted with a poisson distribution whose mean value is determined on the basis of the fluence, opportunely normalized accounting also for the weight of the particle used. Then, in two nested \c for loops over the number of particles extracted and over the #tracks vector, the dose deposited in the nucleus is evaluated by means of the Nucleus::distributeDose() method. At the end of the loop, dose and survival (with respective uncertainties) are updated by calling the Nucleus::getDoseAndSurvival() method.
         
            \note Parallelism: If the parallelism is disabled, the same nucleus could be irradiated and cleaned recursively in order to save memory; but if the parallelism is used, each thread must work on a different nucleus, hence this method provide the possibility to indicate also the nucleus to be irradiated with the \c nuc_cp parameter. This is the only difference between this method and histogram_dose_survival() which is actually useless. It was not deleted only for a chronological reason, because the parallelism was implemented later, but at present it is not used.
         
            \param doseImposed The dose imposed to be delivered to the nucleus, expressed in Gy.
            \param dose The dose absorbed by the nucleus in the irradiation, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival observed, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
            \param nuc_cp A reference to the nucleus to be irradiated.
            \param clean A boolean value indicating if the nucleus has to be cleaned at the end of the evaluation.
         
            \sa random_dose_survival_p() and histogram_dose_survival_t()
         */
        void histogram_dose_survival_p(const double doseImposed,
                                       double &dose,
                                       double &doseUncertainty,
                                       double &survival,
                                       double &survivalUncertainty,
                                       Nucleus &nuc_cp,
                                       bool clean = true);
        
        //! Simulates via the Monte Carlo method the irradiation process of a nucleus, which is part of the cellular population, in a single fraction with a defined duration.
        /*!
            First of all the function calculate the fluence of the beam on the base of the dose imposed, by means of the following relation:
            \f[
                \Phi=\rho\frac{D}{K\cdot \langle LET\rangle}
            \f]
            where D is the dose imposed, K is a constant to accounts for the different units of measure used, \f$\rho\f$ is the density of the target and \f$\langle LET\rangle\f$ is the mean value of LET of the radiation (that coincides with the LET of the particle only in the monoenergetic case).
            The number of particles interacting with the nucleus is randomly extracted with a poisson distribution whose mean value is determined on the basis of the fluence, opportunely normalized accounting also for the weight of the particle used. Before irradiating, the time structure of the irradiation has to be defined: starting from an instant indicated as parameter of the function (\f$t_0\f$), a sequence of times is randomly extracted with uniform probability in \f$[t_0,\;t_0+\Delta t]\f$ (where \f$\Delta t\f$ indicates the fraction delivery time); each time extracted is associated to a particle as its temporal index.
            Then, in two nested \c for loops over the number of particles extracted and over the #tracks vector, the dose deposited in the nucleus is evaluated by means of the Nucleus::distributeDose() method.
         
            \param nuc_cp A reference to the nucleus to be irradiated.
            \param doseImposed The dose to be delivered in the irradiation in a single fraction, expressed in Gy.
            \param time The starting point of the temporal sequence, expressed in hours.
            \param fracDeliveryTime The time needed to deliver a single fraction, expressed in hours.
         
            \sa evaluateG() and histogram_dose_survival_p()
         */
        void histogram_dose_survival_t(Nucleus &nuc_cp,
                                       const double doseImposed,
                                       const double time,
                                       const double fracDeliveryTime);
        
        //! Simulates, via the Monte Carlo method, the irradiation process of a nucleus which is part of the cellular population. This function was thought to control the dose delivery inside the MKM nucleus. It's similar to histogram_dose_survival() but provide also informations on the microscopic dose deposition pattern in each domain of the nucleus.
        /*!
            First of all the function calculate the fluence of the beam on the base of the dose imposed, by means of the following relation:
            \f[
                \Phi=\rho\frac{D}{K\cdot \langle LET\rangle}
            \f]
            where D is the dose imposed, K is a constant to accounts for the different units of measure used, \f$\rho\f$ is the density of the target and \f$\langle LET\rangle\f$ is the mean value of LET of the radiation (that coincides with the LET of the particle only in the monoenergetic case).
            The number of particles interacting with the nucleus is randomly extracted with a poisson distribution whose mean value is determined on the basis of the fluence, opportunely normalized accounting also for the weight of the particle used. Then, in two nested \c for loops over the number of particles extracted and over the #tracks vector, the dose deposited in the nucleus is evaluated by means of the Nucleus::distributeDose() method. At the end of the loop, dose and survival (with respective uncertainties) are updated by calling the Nucleus::getDoseAndSurvival() method. It returns also the microscopical informations on doses deposited and lethal events observed in each domain by overwriting the correspondent parameters.
         
            \param doseImposed The dose imposed to be delivered to the nucleus, expressed in Gy.
            \param dose The dose absorbed by the nucleus in the irradiation, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival observed, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
            \param doses The vector of doses absorbed by each domain, expressed in Gy, passed by reference to be overwritten.
            \param lethals The vector of lethal events observed in each domain, passed by reference to be overwritten.
            \param dosesUncertainty The uncertainties associated to the doses absorbed, expressed in Gy, passed by reference to be overwritten.
            \param lethalsUncertainty The uncertainties associated to the lethal events observed, passed by reference to be overwritten.
            \param clean A boolean value indicating if the nucleus has to be cleaned at the end of the evaluation.
         
            \sa histogram_dose_survival_p()
         */
        void histogram_dose_survival_with_domains(const double doseImposed,
                                                  double &dose,
                                                  double &doseUncertainty,
                                                  double &survival,
                                                  double &survivalUncertainty,
                                                  std::vector<double> &doses,
                                                  std::vector<double> &lethals,
                                                  std::vector<double> &dosesUncertainty,
                                                  std::vector<double> &lethalsUncertainty,
                                                  bool clean = true);
        
        //! Simulates, via the Monte Carlo method, the irradiation process of a nucleus which is part of the cellular population.
        /*!
            First of all the function calculate the fluence of the beam on the base of the dose imposed, by means of the following relation:
            \f[
                \Phi=\rho\frac{D}{K\cdot \langle LET\rangle}
            \f]
            where D is the dose imposed, K is a constant to accounts for the different units of measure used, \f$\rho\f$ is the density of the target and \f$\langle LET\rangle\f$ is the mean value of LET of the radiation (that coincides with the LET of the particle only in the monoenergetic case).
            The number of particles interacting with the nucleus is randomly extracted with a poisson distribution whose mean value is determined on the basis of the fluence, opportunely normalized.
            Then, in a \c for loops over the total number of particles extracted, a different particle from the #tracks vector is randomly chosen (iteration by iteration) with uniform probability and the dose deposited in the nucleus is evaluated by means of the Nucleus::distributeDose() method. At the end of the loop, dose and survival (with respective uncertainties) are updated by calling the Nucleus::getDoseAndSurvival() method.
         
            \note Parallelism: If the parallelism is disabled, the same nucleus could be irradiated and cleaned recursively in order to save memory; but if the parallelism is used, each thread must work on a different nucleus, hence this method provide the possibility to indicate also the nucleus to be irradiated with the \c nuc_cp parameter. This is the only difference between this method and random_dose_survival() which is actually useless. It was not deleted only for a chronological reason, because the parallelism was implemented later, but at present it is not used.
         
            \param doseImposed The dose imposed to be delivered to the nucleus, expressed in Gy.
            \param dose The dose absorbed by the nucleus in the irradiation, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival observed, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
            \param nuc_cp A reference to the nucleus to be irradiated.
            \param clean A boolean value indicating if the nucleus has to be cleaned at the end of the evaluation.
         
            \sa random_dose_survival(), histogram_dose_survival_p() and slow_meanDose_meanSurvival()
         */
        void random_dose_survival_p(const double doseImposed,
                                    double &dose,
                                    double &doseUncertainty,
                                    double &survival,
                                    double &survivalUncertainty,
                                    Nucleus &nuc_cp,
                                    bool clean = true);
        
        //! This method was developed to improve both the \f$\alpha\f$ and \f$\beta\f$ estimation with respect to the approach of Scholz. It is currently applicable to the LEM I version only (the possibility of extending the approach to the subsequent versions is under investigation).
        /*!
            This method is based on a quasi-analytic solution of the LEM for the \f$\alpha_P\f$ and \f$\beta_P\f$ parameters.
         
            \warning It hasn't been implemented yet.
         */
        void rapidINFN_alphaIon_betaIon(double &alphaIon,
                                        double &betaIon);
        
        //! Fast implementation of the MKM as described in (\ref Kase_2008)
        /*!
         Within this approach, the estimate for \f$\alpha\f$ is obtained as:
         \f[
        \alpha=\frac{1-\exp(-\alpha_P\,\gamma_{nucleus})}{\gamma_{nucleus}}
        \f]
        where \f$\gamma_{nucleus}\f$ is the dose-weighted average of the specific energy deposited in the nucleus
        by a single track, evaluated by means of the approximated formula:
        \f[
        \gamma_{nucleus}=\frac{LET}{\rho\,\sigma}
        \f]
        while
        \f[
        \alpha_P=\alpha_0+\beta_0\,\gamma
        \f]
        indicating with \f$\alpha_0\f$ and \f$\beta_0\f$ the input LQ parameters of the MKM, which can be identified
         with the LQ parameters for X-ray reference irradiation,
        and with \f$\gamma\f$ the dose-weighted average of the dose deposited by a single event in the domain, or:
        \f[
        \gamma=\frac{\left\langle z_d^2\right\rangle}{\left\langle z_d\right\rangle}
        \f]
        where \f$z_d\f$ is the specific energy deposited in the single domain in each interaction event.
        For \f$\beta\f$ no recipe is available. In this implementation it is assumed to be constant and
        equal to \f$\beta_0\f$, even if this contrasts with most of the experimental data.
        To accounts also for the mixed fields, the method estimates \f$\alpha\f$ and \f$\beta\f$
        for each tracks of the #tracks vector accounting also for the particle weight (Particle::weight).
        Then \f$\alpha\f$ and \f$\beta\f$ are evaluating according to the TDRA:
        \f[
        \alpha=\frac{\sum_i\,\alpha_i\,LET_i}{\sum_i\,LET_i}
        \f]
        \f[
        \sqrt{\beta}=\frac{\sum_i\,\sqrt{\beta_i}\,LET_i}{\sum_i\,LET_i}
        \f]
        
        \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
        \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
        
        \sa slow_alphaIon_betaIon() and rapidMKM_Attili2013()
        
        \anchor Kase_2008 Kase, Y., Kanai, T., Matsufuji, N., Furusawa, Y., Elsässer, T., & Scholz, M. (2008).
        Biophysical calculation of cell survival probabilities using amorphous track structure models
        for heavy-ion irradiation. \a Physics \a in \a Medicine \a and \a Biology, 53(1), 37–59
         */
        void rapidMKM_Kase2008(double &alphaIon, double &betaIon);
        
        //! Extension of the rapidMKM_Kase2008() method with \f$\beta=\beta(\mathrm{LET})\f$
        /*!
         This implementation correspond to the implementation of rapidMKM_Kase2008() with the LET-dependent non-poissonian correction factor
         (\f$alpha/alpha_P\f$) applied to the quadratic term (\f$\beta_0\f$) too:
         \f[
         \beta=\left(\frac{\alpha}{\alpha_P}\right)^2\beta_0
         \f]
         */
        void rapidMKM_Kase2008_corrected_beta(double &alphaIon, double &betaIon);
        
        //! This method provide a fast original implementation of the MKM model, combining the methods described in (\ref Hawkins_2003) and (\ref Kase_2008).
        /*!
            Whithin this method the non-Poissonian corrective factor introduced in (\ref Hawkins_2003) is exactly evaluated using
            the track model adopted by (\ref Kase_2008) and performing an explicit integration of the track averaged dose in the cell nucleus.
         
            Within this approach, the estimate for \f$\alpha\f$ is obtained as:
            \f[
                \alpha=\frac{1-\exp(-\alpha_P\,\gamma_{nucleus})}{\gamma_{nucleus}}
            \f]
            where \f$\gamma_{nucleus}\f$ is the dose-weighted average of the specific energy deposited in the nucleus
            by a single track, evaluated by integrating:
            \f[
                \gamma_{nucleus}=\frac{\left\langle z_n^2\right\rangle}{\left\langle z_n\right\rangle}
            \f]
            where \f$z_n\f$ is the specific energy deposited in the single domain in each interaction event, while
            \f[
                \alpha_P=\alpha_0+\beta_0\,\gamma
            \f]
            indicating with \f$\alpha_0\f$ and \f$\beta_0\f$ the input LQ parameters of the model, which can be identified with
            the LQ parameter of the reference X-ray irradiation.
            \f$\gamma\f$ is the dose-weighted average of the dose deposited by a single event in the domain, or:
            \f[
                \gamma=\frac{\left\langle z_d^2\right\rangle}{\left\langle z_d\right\rangle}
            \f]
            where \f$z_d\f$ is the specific energy deposited in the single domain in each interaction event.
            For \f$\beta\f$ no recipe is available. In most of the applications it is assumed to be constant and
            equal to \f$\beta_X\f$, even if this contrasts with most of the experimental data.
            To accounts also for the mixed fields, the method estimates \f$\alpha\f$ and \f$\beta\f$
            for each tracks of the #tracks vector accounting also for the particle weight (Particle::weight).
            Then \f$\alpha\f$ and \f$\beta\f$ are evaluating according to the TDRA:
            \f[
                \alpha=\frac{\sum_i\,\alpha_i\,LET_i}{\sum_i\,LET_i}
            \f]
            \f[
                \sqrt{\beta}=\frac{\sum_i\,\sqrt{\beta_i}\,LET_i}{\sum_i\,LET_i}
            \f]
         
            \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
         
            \sa slow_alphaIon_betaIon() and rapidMKM_Kase2008()
         
            \anchor Hawkins_2003 Hawkins, R. B. (2003). A microdosimetric-kinetic model for the effect of non-Poisson distribution of lethal
                    lesions on the variation of RBE with LET. \a Radiation \a Research, 160(1), 61–69.
         
            \anchor Kase_2008 Kase, Y., Kanai, T., Matsufuji, N., Furusawa, Y., Elsässer, T., & Scholz, M. (2008).
                    Biophysical calculation of cell survival probabilities using amorphous track structure models
                    for heavy-ion irradiation. \a Physics \a in \a Medicine \a and \a Biology, 53(1), 37–59
         */
        void rapidMKM_Attili2013(double &alphaIon, double &betaIon);
        
        //! Extension of the rapidMKM_Attili2013() method with \f$\beta=\beta(\mathrm{LET})\f$
        /*!
         This implementation correspond to the implementation of rapidMKM_Attili2013() with the LET-dependent non-poissonian correction factor
         (\f$alpha/alpha_P\f$) applied to the quadratic term (\f$\beta_0\f$) too:
         \f[
           \beta=\left(\frac{\alpha}{\alpha_P}\right)^2\beta_0
         \f]
         */
        void rapidMKM_Attili2013_corrected_beta(double &alphaIon, double &betaIon);
        
        //! This method is based on the approach of Scholz (rapidScholz_alphaIon_betaIon()) but provides a more precise estimation of the \f$\alpha\f$ parameter.
        /*!
            This method was proposed by the INFN in 2011, with the work of Russo (\ref rapidRusso "1").
            To improve the rapidScholz_alphaIon_betaIon() method it is necessary to include the effect of ions that, passing near the nucleus, irradiate it only through their track penumbra.
         
            Here, the rigorous derivations of the formulas used is omitted, the reader interested look at the work of Russo (\ref rapidRusso "1").
         
            The LQ \f$\alpha\f$ parameter could be calculated in the monoenergetic case as:
            \f[
                \alpha=\alpha_X\left(1-\langle z\rangle_{dir}\frac{\rho A_{nucl}}{LET}\right)+(1-\langle S\rangle_{dir})\frac{\rho A_{nucl}}{LET}
            \f]
            where \f$\rho\f$ represents the density of the medium, \f$A_{nucl}\f$ the area of the nucleus, \f$\langle z\rangle_{dir}\f$ and \f$\langle S\rangle_{dir}\f$ are the single-event dose and survival corresponding to the ion traversing the nucleus at its center and \f$\alpha_X\f$ is the linear quadratic \f$\alpha\f$-parameter characteristic for X-rays.
            The \f$\beta\f$ parameter, as in the case of the approach of Scholz, is evaluated as
            \f[
                \beta=\left(\frac{\alpha}{\alpha_P}\right)^2\beta_P
            \f]
            where, in this case, diversely from the approach of Scholz:
            \f[
                \alpha_P=\alpha_X\left(1-\langle z\rangle_{dir}\frac{\rho A_{nucl}}{LET}\right)-\ln(\langle S\rangle_{dir})\frac{\rho A_{nucl}}{LET}
            \f]
            and
            \f[
                \beta_P=\frac{s-\alpha_P}{2D_t}
            \f]
            (see the LEM II parametrization for \em s and \f$D_t\f$, CellLine::parametrization_LQ2()).
            To accounts also for the mixed fields, the method evaluates \f$\langle S\rangle_{dir}\f$ for each tracks of the #tracks vector estimating its \f$\alpha\f$ accounting also for the particle weight (Particle::weight).
            Then \f$\alpha\f$ and \f$\beta\f$ are evaluating according to the theory of dual radiation action (\ref TDRA "2"):
            \f[
                \alpha=\frac{\sum_i\,\alpha_i\,LET_i}{\sum_i\,LET_i}
            \f]
            \f[
                \sqrt{\beta}=\frac{\sum_i\,\sqrt{\beta_i}\,LET_i}{\sum_i\,LET_i}
            \f]
            The function provide also the possibility to plot the results.
         
            \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
         
            \sa rapidScholz_alphaIon_betaIon() and rapidINFN_alphaIon_betaIon()
         
            \anchor rapidRusso 1. G. Russo, "Develpment of a radiobiological database for carbon ion Treatment Planning Systems - Modelling and simulating the irradiation process", \a PhD \a Thesis, Università degli studi di Torino (2011).
         
            \anchor TDRA 2. M. Zaider and H.H. Rossi, "The synergistic effets of different radiations", \a Radiation \a Research \b 160, 61-69 (2003).
         */
        void rapidLEM_Russo2011(double &alphaIon, double &betaIon);
        
        //! This method provide a faster approximate implementation of the LEM model that avoids the Monte Carlo simulation.
        /*!
            This method was proposed by Krämer and Scholz in 2006 (\ref rapidLEM "1") and it is applicable only to the case of monoenergetic irradiation, but the estimation is extended to the mixed-field case exploiting the Theory of Dual Radiation Action (TDRA) by Zaider and Rossi (\ref TDRA "2").
            The method assumes that only direct impacts of ions on the cell nucleus are relevant for the evaluation of the \f$\alpha\f$ parameter that could be calculated (in the monoenergetic case) as:
            \f[
                \alpha=\frac{\rho A_{nucl}}{LET}(1-\langle S\rangle_{dir})
            \f]
            where \f$\rho\f$ represents the density of the medium, \f$A_{nucl}\f$ the area of the nucleus and \f$\langle S\rangle_{dir}\f$ is the single-event survival corresponding to the ion traversing the nucleus at its center, disregarding the small dependence over the ion impact parameter.
            The \f$\beta\f$ parameter is estimated as:
            \f[
                \beta=\left(\frac{\alpha}{\alpha_P}\right)^2\beta_P
            \f]
            where
            \f[
                \alpha_P=-\ln(\langle S\rangle_{dir})\frac{\rho A_{nucl}}{LET}
            \f]
            and
            \f[
                \beta_P=\frac{s-\alpha_P}{2D_t}
            \f]
            (see the LEM II parametrization for \em s and \f$D_t\f$, CellLine::parametrization_LQ2()).
            To accounts also for the mixed fields, the method evaluates \f$\langle S\rangle_{dir}\f$ for each tracks of the #tracks vector estimating its \f$\alpha\f$ accounting also for the particle weight (Particle::weight).
            Then \f$\alpha\f$ and \f$\beta\f$ are evaluating according to the TDRA:
            \f[
                \alpha=\frac{\sum_i\,\alpha_i\,LET_i}{\sum_i\,LET_i}
            \f]
            \f[
                \sqrt{\beta}=\frac{\sum_i\,\sqrt{\beta_i}\,LET_i}{\sum_i\,LET_i}
            \f]
            The function provide also the possibility to plot the results.
         
            \note For a rigorous derivation of the formulas used see the published reference.
         
            \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
         
            \sa rapidRusso_alphaIon_betaIon() and rapidINFN_alphaIon_betaIon()
         
            \anchor rapidLEM 1. M. Krämer and M. Scholz, "Rapid calculation of biological effects in ion radiotherapy", \a Physics \a in \a medicine \a and \a biology \b 51, 1959-1970 (2006).
         
            \anchor TDRA 2. M. Zaider and H.H. Rossi, "The synergistic effets of different radiations", \a Radiation \a Research \b 160, 61-69 (2003).
         */
        void rapidLEM_Scholz2006(double &alphaIon, double &betaIon);
        
        //! Sets the number of threads.
        /*!
            \param nTh The number of threads to be set.
         
            \sa #nThreads
         */
        void setNThreads(int nTh){nThreads=nTh;};
        
        //! Sets the prefix of the output file name.
        /*!
            \param save_prefix The prefix of the output file name.
         
            \sa #savePrefix
         */
        void setSavePrefix(std::string save_prefix){savePrefix=save_prefix;};
        
        //! Method called to perform a Monte Carlo simulation to reproduce the irradiation process getting the LQ parameters \f$\alpha\f$ and \f$\beta\f$.
        /*!
            The idea is to simulate the irradiation process of an entire cell population and to repeat the simulation for different values of dose imposed in order to obtain a complete survival curve. For each value of dose imposed the mean value of dose absorbed and cellular survival observed in the population are obtained by calling the slow_meanDose_meanSurvival() method. Then the survival curve is fitted by means of the fit_LQ() method in order to get the LQ parameters. Finally the function returns these parameters with the associated uncertainties by overwriting the correspondent variables passed by reference.
         
            \param trackMode A \c string defining the modality to pass the vector of particles in the mixed fields case. The possibilities are "histogram" or "random".
         
            \param parameters The vector containing the model parameters to be used in the simulation (\ref pars "1").
            \param dosesImposed A vector containing the values of nominal dose to be simulated, expressed in Gy.
            \param precision Fix the ending condition of the Monte Carlo simulation.
            \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param alphaIonUncertainty The uncertainty associated to the \f$\alpha\f$ parameter (in \f$Gy^{-1}\f$), passed by reference to be overwritten.
            \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
            \param betaIonUncertainty The uncertainty associated to the \f$\beta\f$ parameter (in \f$Gy^{-2}\f$), passed by reference to be overwritten.
            \param nFraction The total number of fraction, in case of fractionated treatment.
            \param timeSpacing The time spacing between fractions, expressed in hours.
            \param fracDeliveryTime The delivery time of each fraction, expressed in hours.
            \param saveAlphaBeta A boolean parameter indicating if the extrapolated \f$\alpha\f$ and \f$\beta\f$ parameters are to be saved.
            \param saveMeans A boolean parameter indicating if the informations on mean dose and survival observed are to be saved.
            \param saveCell A boolean parameter indicating if the dose-survival data of each cell irradiated are to be saved.
            \param title_means The title of the file where the method will save the informations on mean dose and survival resultin from the simulation.
         
            \warning The execution of the program will be terminated if the minimum dose imposed is greater than the maximum one or if an inexistent track mode is selected.
         
            \sa slow_alphaIon_betaIon_with_Domains() and histogram_dose_survival()
         
            \anchor pars 1. These parameters are stored in the CellLine object but the way to get them is a little bit tricky, hence this is an easier and not aestethically perfect way to get these informations. This has to be fixed in the next versions of the program.
         */
        void slow_alphaIon_betaIon(const std::string trackMode,
                                   const std::vector<double> parameters,
                                   const std::vector<double> dosesImposed,
                                   const double precision,
                                   double &alphaIon,
                                   double &alphaIonUncertainty,
                                   double &betaIon,
                                   double &betaIonUncertainty,
                                   const int nFraction,
                                   const double timeSpacing,
                                   const double fracDeliveryTime,
                                   const bool saveAlphaBeta,
                                   const bool saveMeans,
                                   const bool saveCell,
                                   const std::string title_means);
        
        //! This function was thought to control the dose delivery inside the MKM nucleus. It's similar to slow_alphaIon_betaIon() but provide also informations on the microscopic dose deposition pattern in each domain of the nucleus.
        /*!
            The idea is to simulate the irradiation process of an entire cell population and to repeat the simulation for different values of dose imposed in order to obtain a complete survival curve. For each value of dose imposed the mean value of dose absorbed and cellular survival observed in the population are obtained by calling the slow_meanDose_meanSurvival_with_Domains() method. Then the survival curve is fitted by means of the fit_LQ() method in order to get the LQ parameters. Finally the function returns these parameters with the associated uncertainties by overwriting the correspondent variables passed by reference.
         
            \param trackMode Defined for completeness. Only "histogram" is supported.
            \param minDose The minimum value of dose imposed expressed in Gy.
            \param maxDose The maximum value of dose imposed expressed in Gy.
            \param numberOfDoses The total number of doses imposed (together with \c minDose and \c maxDose this allow to define the sequence).
            \param precision Fix the ending condition of the Monte Carlo simulation.
            \param alphaIon The LQ \f$\alpha\f$ parameter expressed in \f$Gy^{-1}\f$, passed by reference to be overwritten.
            \param alphaIonUncertainty The uncertainty associated to the \f$\alpha\f$ parameter (in \f$Gy^{-1}\f$), passed by reference to be  overwritten.
            \param betaIon The LQ \f$\beta\f$ parameter expressed in \f$Gy^{-2}\f$, passed by reference to be overwritten.
            \param betaIonUncertainty The uncertainty associated to the \f$\beta\f$ parameter (in \f$Gy^{-2}\f$), passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if the minimum dose imposed is greater than the maximum one or if an inexistent track mode is selected.
         
            \warning Even if not explicitly required, the nucleus has to be a Nucleus_MKM object.
         
            \sa slow_alphaIon_betaIon() and histogram_dose_survival_with_domains()
         */
        void slow_alphaIon_betaIon_with_Domains(const std::string trackMode,
                                                const double minDose,
                                                const double maxDose,
                                                const int numberOfDoses,
                                                const double precision,
                                                double &alphaIon,
                                                double &alphaIonUncertainty,
                                                double &betaIon,
                                                double &betaIonUncertainty);
        
        //! Perform the Monte Carlo simulation relatively to a single value of dose imposed and returns the mean values of dose absorbed and survival observed in the cell population.
        /*!
            The user has to fix the precision required for the simulation, that is the statistics to be reached to terminate the simulation. Two possibilities are supported, so the user can:
                - Fix the number of iterations, hence the \c precision has to be an integer value greater (or at least equal) to 1.
                - Define a constraint on the precision to reach in the simulation in the evaluation of the cell survival (precisely the relative error on the survival), hence the \c precision has to be a \c double in (0, 1).
            Once set the precision desired, the irradiation process of the cell population is simulated by calling the histogram_dose_survival_p() or the random_dose_survival_p() methods in a parallel loop (or histogram_dose_survival_t() in the case of temporal studies), depending on the trackMode selected (see below). Each thread performs the irradiation of a single cell and at the end of each irradiation the method updates the total mean dose and survival observed, with associated uncertainties, and determines if the precision set by the user is reached or not. All the values of dose and survival obtained, cell by cell, are saved in an output file.
         
            \note Survival uncertainty: the estimation of the survival uncertainty is knowingly biased, but extremely efficient and memory-friendly. Note that this quantity is used only to decide when it is possible to stop the simulation hence it doesn't affect the result (and that the bias decrease with increasing iteration).
         
            \note track mode: this is useful in the case of mixed fields, that is when the field is constituted by different particles. In that case the there are different ways to extract from the #tracks vector the particle to be used event by event:
                - If "histogram" is selected then the #tracks vector is interpreted as an histogram where each particle has a specific weight, defined by its frequency in the histogram. Different particle, in this case, are not necessary equiprobable.
                - If "random" is selected then the different tracks in the #tracks vector are considered equiprobable and a random extraction between them is performed event by event.
         
            \param trackMode A \c string defining the modality to pass the vector of particles in the mixed fields case. The possibilities are "histogram" or "random".
            \param doseImposed The value of dose to be delivered in the irradiation, expressed in Gy.
            \param precision The precision to be reached in the simulation (could be a fixed number of iterations or a constraint on the survival precision).
            \param meanDose The mean value of dose absorbed expressed in Gy, passed by reference to be overwritten.
            \param meanDoseUncertainty The uncertainty on meanDose, expressed in Gy, passed by reference to be overwritten.
            \param meanSurvival The mean value of cellular survival observed, passed by reference to be overwritten.
            \param meanSurvivalUncertainty The uncertainty associated to the mean survival, passed by reference to be overwritten.
            \param nFraction The number of fraction in the case of fractionated treatment.
            \param timeSpacing The time spacing between fractions, expressed in hours.
            \param fracDeliveryTime The delivery time of each fraction, expressed in hours.
            \param saveCell A boolean parameter indicating if the dose-survival data of each cell irradiated are to be saved or not.
         
            \warning The execution of the program will be terminated if the precision is not set correctly or if an inexistent track mode is selected.
         
            \sa slow_alphaIon_betaIon()
         */
        void slow_meanDose_meanSurvival(const std::string trackMode,
                                        const double doseImposed,
                                        const double precision,
                                        double &meanDose,
                                        double &meanDoseUncertainty,
                                        double &meanSurvival,
                                        double &meanSurvivalUncertainty,
                                        const int nFraction,
                                        const double timeSpacing,
                                        const double fracDeliveryTime,
                                        const bool saveCell);
        
        //! Perform the Monte Carlo simulation relatively to a single value of dose imposed and returns the mean values of dose absorbed and survival observed in the cell population. This function was thought to control the dose delivery inside the MKM nucleus. It's similar to slow_meanDose_meanSurvival() but provide also informations on the microscopic dose deposition pattern in each domain of the nucleus.
        /*!
            The user has to fix the precision required for the simulation, that is the statistics to be reached to terminate the simulation. Two possibilities are supported, so the user can:
                - Fix the number of iterations, hence the \c precision has to be an integer value greater (or at least equal) to 1.
                - Define a constraint on the precision to reach in the simulation in the evaluation of the cell survival (precisely the relative error on the survival), hence the \c precision has to be a \c double in (0, 1).
            Once set the precision desired, the irradiation process of the cell population is simulated by calling the histogram_dose_survival_with_domains() method in a loop that, in each iteration, performs the irradiation of a single cell. At the end of the irradiation, the method updates the total mean dose and survival observed, with associated uncertainties, and determines if the precision set by the user is reached or not. All the values of dose and survival obtained, cell by cell, are saved in an output file together with the total number of lethal events observed in the nucleus.
         
            \note Survival uncertainty: the estimation of the survival uncertainty is knowingly biased, but extremely efficient and memory-friendly. Note that this quantity is used only to decide when it is possible to stop the simulation hence it doesn't affect the result (and that the bias decrease with increasing iteration).
         
            \param trackMode Defined for completeness. Only "histogram" is supported.
            \param doseImposed The value of dose to be delivered in the irradiation, expressed in Gy.
            \param precision The precision to be reached in the simulation (could be a fixed number of iterations or a constraint on the survival precision).
            \param meanDose The mean value of dose absorbed expressed in Gy, passed by reference to be overwritten.
            \param meanDoseUncertainty The uncertainty on meanDose, expressed in Gy, passed by reference to be overwritten.
            \param meanSurvival The mean value of cellular survival observed, passed by reference to be overwritten.
            \param meanSurvivalUncertainty The uncertainty associated to the mean survival, passed by reference to be overwritten.
         
            \warning The execution of the program will be terminated if the precision is not set correctly or if an inexistent track mode is selected.
         
            \sa slow_alphaIon_betaIon() and slow_alphaIon_betaIon_with_Domains()
         */
        void slow_meanDose_meanSurvival_with_Domains(const std::string trackMode,
                                                     const double doseImposed,
                                                     const double precision,
                                                     double &meanDose,
                                                     double &meanDoseUncertainty,
                                                     double &meanSurvival,
                                                     double &meanSurvivalUncertainty);
        
        //! Evaluates the dose deposited in the nucleus using directly the #tracks vector without modifying it and without random numbers extractions.
        /*!
            It simply calls the Nucleus::distributeDose(const Tracks) method passing the #tracks vector and then it gets informations on dose deposited and survival observed, with associated uncertainties, via the Nucleus::getDoseAndSurvival() method.
         
            \param dose The dose absorbed by the nucleus in the irradiation, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival observed, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
            \param clean A boolean value indicating if the nucleus has to be cleaned at the end of the evaluation.
         
            \note It's actually unused.
         
            \sa histogram_dose_survival_p() and random_dose_survival_p()
         */
        void verbatim_dose_survival(double &dose,
                                    double &doseUncertainty,
                                    double &survival,
                                    double &survivalUncertainty,
                                    bool clean = true);
        
    private:
        
        //! A \c const reference to a Track object corresponding to the Particle interacting with #nucleus.
        const Tracks &tracks;
        
        //! A \c const reference to a CellLine object corresponding to the cell line to which the nucleus belongs.
        const CellLine &cellLine;
        
        //! A reference to the cellular nucleus.
        Nucleus &nucleus;
        
        //! The number of threads needed to be used in the simulation (if parallelism is supported).
        /*!
            Possible cases are:
                - 0: Uses a number of threads corresponding to the number of core of the machine executing this program.
                - 1: Uses 1 threads, i.e. Disabled multithread
                - A number greater than 1: Specifies the exact number of threads.
         */
        int nThreads;
        
        //! A pointer to a gsl_rng object, useful in the generation of pseudorandom numbers in the Monte Carlo simulation.
        gsl_rng *randomGenerator;
        
        //! The prefix of the output file.
        std::string savePrefix;
        
        //! The model used in the simulation.
        std::string model;
    };
    
}


#endif /* CALCULUS_H */
