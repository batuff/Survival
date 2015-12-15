#ifndef NUCLEUS_INTEGRAL_H_
#define NUCLEUS_INTEGRAL_H_

#include "Nucleus.h"

namespace Survival {

    //! Implements a nucleus as a 2D circular object and provides methods to evaluate the number of lethal events observed.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Lorenzo Marengo
        \author Germano Russo
        \date 2011--2015
     
        This class implements a nucleus as a 2D circular object that represents the cross section shown by the cell to the particle. It contains a reference to the cell line to which the nucleus belongs, which is used to get informations such as the radius. It provides method to evaluate the dose deposited by the radiation in the interaction with the nucleus itself and methods to get the dose absorbed and the associated cellular survival.
     
        \sa Nucleus_MKM and Nucleus_Integral_t
     */
    class Nucleus_Integral : public Nucleus
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            \param cellLineRef A reference to the corresponding CellLine.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
         
            \sa cleanNucleus() and Nucleus_MKM::createDomains()
         */
        Nucleus_Integral(const CellLine &cellLineRef,
                         const double xPosition = 0.0,  // nucleus position referred to the beam axis (mm)
                         const double yPosition = 0.0);  // nucleus position referred to the beam axis (mm)
        
        //! Destructor.
        virtual ~Nucleus_Integral() {};
        
        //! Adds a constant value of dose absorbed by the nucleus.
        /*!
            The method updates the #totalNucleusDose value, adding a constant value chosen by the user.
         
            \param dose The value of dose absorbed to be added, expressed in Gy.
         
            \sa Nucleus_MKM::addBackgroundDose()
         */
        void addBackgroundDose(const double dose);
        
        //! Resets to zero #inNucleusCount and #intersectionCount counters and the total dose absorbed (#totalNucleusDose).
        /*!
            \sa Nucleus_MKM::cleanNucleus()
         */
        virtual void cleanNucleus();
        
        //! Returns a pointer to a new Nucleus_Integral object. It not really a clone but a new clean object.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         
            \note To create a real clone of another nucleus, a better implementation of the copy constructor is needed.
         */
        virtual Nucleus_Integral* clone(const CellLine&);
        
        //! Integrates the radial profile of the track in the intersection area with the nucleus to evaluate the dose deposited.
        /*!
            It's a very tricky and clever method that performs the integral of the radial profile of the track in its common area with the nucleus.
            First of all it centers the origin of the reference system (RS) on the position of the nucleus, identifying the position of the track in the new RS. Then it examines all possible cases:
                - Track inside the nucleus;
                - Track outside the nucleus but interacting with it;
                - Track non interacting with the nucleus (in this case it trivially returns nothing).
            In the first two cases the respective counter is incremented (#inNucleusCount or #intersectionCount), then the method identifies the intersection area and calls IntegrateWeightedRadialTrack() that performs the integral in the area defined and evaluates the area itself (unless the Track is all contained in the nucleus: in that case the function Track::getRadialIntegral() is called, and the value of the area is calculated \em manually as \f$\pi R^2\f$).
         
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
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and the associated survival, with respective uncertainties.
        /*!
            The dose absorbed coincides with #totalNucleusDose, the survival is evaluated by means of an exponential function of the lethal events observed, evaluated through the CellLine::getLogSurvival_X() method.
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival associated to the dose absorbed by the nucleus, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and survival, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         */
        virtual void getDoseAndSurvival(double &dose,
                                        double &doseUncertainty,
                                        double &survival,
                                        double &survivalUncertainty) const;
        
        //! Returns the number of times that the nucleus has been crossed through by a Particle.
        /*!
            \return #inNucleusCount The number of times that the nucleus has been crossed through by a Particle.
         
            \sa distributeDose(const Track&)
         */
        virtual int getInNucleusCount() const {return inNucleusCount;};
        
        //! Returns the number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
        /*!
            \return #intersectionCount The number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
         
            \sa distributeDose()
         */
        virtual int getIntersectionCount() const {return intersectionCount;};
        
        //! Returns the nucleus position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the nucleus referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c y coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
         
            \sa Nucleus_MKM::getPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const;
        
        //! Returns the radius of the nucleus expressed in um.
        /*!
            \return #r_nucleus The radius of the nucleus expressed in um.
         */
        virtual double getRadius() const {return r_nucleus;};
        
        
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
        double ArcIntersectionWeight(double r, double b);
        
        //double IntegrateRadialTrack(const Track &track, double rMin, double rMax, double step = 0.03);
        
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
        double IntegrateWeightedRadialTrack(const Track &track, double rMin, double rMax, double b, double &area, double step);
        
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
            It's instantiated in the constructor getting the value from the CellLine object representing the cell line to which the nucleus belongs.
         
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
        
        //! The number of times that the nucleus interacted with a Particle that doesn't pass through the nucleus itself.
        /*!
            It's incremented by means of the distributeDose(const Track&) method.
         */
        int intersectionCount;
    };
    
}


#endif /* NUCLEUS_INTEGRAL_H_ */
