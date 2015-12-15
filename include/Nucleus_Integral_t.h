#ifndef NUCLEUS_INTEGRAL_T_H_
#define NUCLEUS_INTEGRAL_T_H_

#include "Nucleus.h"

namespace Survival {

    //! Implements a nucleus as a 2D circular object and provides methods to evaluate the number of lethal events observed taking into account the time structure of the irradiation.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2015
     
        This class inherits from the Nucleus pure virtual class and it has the same structure of the Nucleus_Integral class. It defines the nucleus used in the MCt-MKM model and its peculiarity is to provide some methods to evaluate the temporal effect of the irradiation keeping track of its time structure. Every time an interaction between the nucleus and a particle occurs, the time of the event and the dose deposited are stored in dedicated vectors to be passed at the Nucleus_tMKM class.
     */
    class Nucleus_Integral_t : public Nucleus
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            \param cellLineRef A reference to the corresponding CellLine.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa cleanNucleus() and Nucleus_tMKM::createDomains()
         */
        Nucleus_Integral_t(const CellLine &cellLineRef,
                           const double xPosition = 0.0,
                           const double yPosition = 0.0);
        
        //Nucleus_Integral_t(const Nucleus_Integral_t& nn, //Questo è scritto bene e compila
        //                   const CellLine& cellLineRef);
        
        //! Destructor.
        virtual ~Nucleus_Integral_t() {};
        
        //! Adds a constant value of dose absorbed by the nucleus in a specific instant.
        /*!
            The method added an element to #doses and #times vectors, representing a constant dose absorbed by the nucleus in a specific instant.
            The method updates also the #totalNucleusDose value, adding a constant value chosen by the user.
         
            \param dose The value of dose absorbed to be added, expressed in Gy.
            \param t The time associated to the dose added, expressed in hours.
         
            \sa Nucleus_tMKM::addBackgroundDose()
         */
        void addBackgroundDose(const double dose,
                               const double t);
        
        //! Resets to zero #inNucleusCount and #intersectionCount counters, the total dose absorbed (#totalNucleusDose) and the vectors containing the history of times and doses deposited.
        /*!
            \sa Nucleus_tMKM::cleanNucleus()
         */
        virtual void cleanNucleus();
        
        //! Returns a pointer to a new Nucleus_Integral_t object. It not really a clone but a new clean object.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         
            \note To create a real clone of another nucleus, a better implementation of the copy constructor is needed.
         */
        virtual Nucleus_Integral_t* clone(const CellLine&);
        
        //! Integrates the radial profile of the track in the intersection area with the nucleus to evaluate the dose deposited.
        /*!
            It has the same implementation of Nucleus_Integral::distributeDose(), but once the dose deposited is calculated this method adds an element to #doses and #times vectors (through the \c push.back() method defined in the STL): the dose is the one calculated by the function while the time is got by means of the Track::getTime() method.
         
            \param track The Track of the particle interacting with the nucleus.
         */
        virtual void distributeDose(const Track &track);
        
        //! Overload of distributeDose(const Track &track) to manage a Tracks object, it simply calls distributeDose(const Track &track) for every track of the container.
        /*!
            Since the Tracks class is a container for Track objects, this method calls distributeDose(const Track &track) in a \c for loop over each track contained in the Tracks object.
         */
        virtual void distributeDose(const Tracks &tracks);
        
        //! Returns the name of the cell line to which the nucleus belongs.
        /*!
            \return A \c string corresponding to the name of the cell line to which the nucleus belongs, getting the information by #cellLine.
         */
        virtual std::string getCellType() const;
        
        //! Returns the total dose absorbed by the nucleus, expressed in Gy, overwriting a double variable passed by reference.
        /*!
            \sa #totalNucleusDose
         */
        void getDose(double &dose) const {dose = totalNucleusDose;};
        
        //! Returns the vector representing each dose deposited in the nucleus by a Particle when an interaction occurs.
        /*!
            \return A vector representing each dose deposited in the nucleus by a Particle when an interaction occurs, expressed in Gy.
         
            \sa #doses
         */
        std::vector<double> getDoses() const {return doses;};
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and the associated survival, with respective uncertainties.
        /*!
            The dose absorbed coincides with #totalNucleusDose, the survival is evaluated by means of the CellLine::getLogSurvival_X() method.
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival associated to the dose absorbed by the nucleus, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and survival, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \warning The implementation of this method is identical to the one defined in Nucleus_Integral::getDoseAndSurvival(); it isn't considered the time structure of the irradiation. The reason is that this class was thought to implement the MCt-MKM model, therefore the time structure is considered directly in the Nucleus_tMKM class.
         */
        virtual void getDoseAndSurvival(double &dose,
                                        double &doseUncertainty,
                                        double &survival,
                                        double &survivalUncertainty) const;
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and associated lethal events, with respective uncertainties.
        /*!
            The dose absorbed coincides with #totalNucleusDose, the number of lethal events is evaluated by means of the CellLine::getLogSurvival_X() method.
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param lethals The lethal events observed in the nucleus, passed by reference to be overwritten.
            \param lethalsUncertainty The uncertainty associated to the lethal events observed, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and lethal events, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \warning This method doesn't consider the time structure of the irradiation. The reason is that this class was thought to implement the MCt-MKM model, therefore the time structure is considered directly in the Nucleus_tMKM class.
         */
        void getDoseAndLethals(double &dose,
                               double &doseUncertainty,
                               double &lethals,
                               double &lethalsUncertainty) const;
        
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
        
        //! Returns the nucleus position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the nucleus referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c y coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
         
            \sa Nucleus_tMKM::getPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const;
        
        //! Returns the radius of the nucleus expressed in um.
        /*!
            \return #r_nucleus The radius of the nucleus expressed in um.
         */
        virtual double getRadius() const {return r_nucleus;};
        
        //! Returns the vector representing each instant in which an interaction with a Particle occurs. Each value is expressed in hours.
        /*!
            \return A vector representing each instant in which an interaction with a Particle occurs, expressed in hours.
         
            \sa times
         */
        std::vector<double> getTimes() const {return times;};
        
    private:
        
        //! Evaluate the length of the arc of circumference (expressed in radians) derived from the intersection between the nucleus and a circumference whose radius is equal to \c r and whose center is far \c b from the center of the nucleus.
        /*!
            The calculus is performed by considering all possible cases and relative subcases:
                -# if \c b is smaller than the radius of the nucleus (\f$R_N\f$), i.e. center inside the nucleus, and:
                    - if \f$r<R_N-b\f$ returns \f$2\pi\f$ (full circle);
                    - else returns the length of the arc, that is \f$2\arccos\left(\frac{b}{2r}+\frac{r}{2b}-\frac{R_N^2}{2br}\right)\f$.
                -# if the center is outside the nucleus and:
                    - if \f$r<b-R_N\f$ no intersection occurs, therefore it returns 0;
                    - else if \f$r<b+R_N\f$ it returns the length of the arc, that is \f$2\arccos\left(\frac{b}{2r}+\frac{r}{2b}-\frac{R_N^2}{2br}\right)\f$.
            Else no intersection occurs and it returns 0.
         
            \param r The radius of the circle, expressed in um.
            \param b The impact parameter of the track, expressed in um.
         
            \return The length of the arc of circumference, expressed in radians.
         
            \sa IntegrateWeightedRadialTrack() and distributeDose()
         */
        double ArcIntersectionWeight(double r,
                                     double b);
        
        //! Performs the integral of the radial profile of the track in the intersection area with the nucleus.
        /*!
            The problem is to integrate a function over a complex area originated by the random intersection of two circles.
            The way this method performs the task is to evaluate it numerically, dividing the area (or the radius to be covered) in a number of finite (small) step and evaluating for each step the length of the arc of circumference by means of the ArcIntersectionWeight() method. The sum of all these lengths is equal to the intersection area and it's overwritten in the correspondent parameter. For each step, defined by a specific radius, the method gets the local dose from the track (Track::getLocalDose()) and the integral is evaluated considering a step function constructed in this way. Finally the value of the integral is normalized over the intersection area.
         
            \param track A reference to the Track of the particle interacting with the nucleus.
            \param rMin Minimum radius, lower limit of integration, expressed in um.
            \param rMax Maximum radius, upper limit of integration, expressed in um.
            \param b The impact parameters of the track, expressed in um.
            \param area The intersection area between Track and Nucleus, passed by reference to be overwritten with its correct value.
            \param step The length of the radial step of integration.
         
            \return The value of the integral normalized over the intersection area, expressed in Gy.
         
            \sa distributeDose()
         */
        double IntegrateWeightedRadialTrack(const Track &track,
                                            double rMin,
                                            double rMax,
                                            double b,
                                            double &area,
                                            double step);
        
        //! Vector containing the sequence of interaction times (expressed in hours), each elements is associated to one interaction.
        /*!
            It's updated by distributeDose() every times an interaction occurs. This vector has the same length of #doses: they are strongly coupled because they represent the same event.
         
            \sa getTimes()
         */
        std::vector<double> times;
        
        //! Vector containing the sequence of doses deposited in the nucleus, expressed in Gy, each elements is associated to one interaction.
        /*!
            It's updated by distributeDose() every times an interaction occurs. This vector has the same length of #times: they are strongly coupled because they represent the same event. The sum of the elements of the vector is equal to #totalNucleusDose, representing the total dose absorbed by the nucleus.
         
            \sa getDoses()
         */
        std::vector<double> doses;
        
    protected:
        
        //! The total dose absorbed by the nucleus, expressed in Gy.
        /*!
            The value is initially set to zero in the constructor (through the cleanNucleus() method) and then updated by distributeDose() when an interaction with a Particle occurs.
         */
        double totalNucleusDose;
        
        //! A reference to a CellLine object where the characteristics of the cell line to which the nucleus belongs are stored.
        const CellLine &cellLine;
        
        //! The radius of the nucleus, expressed in um.
        /*!
            It's instantiated in the constructor getting the value from the CellLine object representing the cell line to which the nucleus belongs to.
         
            \sa Nucleus_Integral()
         */
        double r_nucleus;
        
        //! The position of the center of the nucleus (\c x coordinate) referred to the beam axis, expressed in mm.
        const double x_nucleus;
        
        //! The position of the center of the nucleus (\c y coordinate) referred to the beam axis, expressed in mm.
        const double y_nucleus;
        
        //! The number of times that the nucleus has been crossed through by a Particle.
        /*!
            It's incremented by means of the distributeDose(const Track&) method.
         */
        int inNucleusCount;
        
        //! The number of times that the nucleus interacted with a Particle.
        /*!
            It's incremented by means of the distributeDose(const Track&) method.
         */
        int intersectionCount;
    };
    
}


#endif /* NUCLEUS_INTEGRAL_T_H_ */
