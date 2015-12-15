#ifndef NUCLEUS_PIXEL_H
#define NUCLEUS_PIXEL_H

#include "Nucleus.h"

namespace Survival {

    //! Implements the Pixel features to be used in the Nucleus_Pixel class.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007
     
        The object has simply a few public data members, like the pixel position or the dose deposited. The most important data members are #v and #numberOfSubPixels that identify the real structure. The Nucleus_Pixel is divided into four grid of pixels of decreasing dimension (this is useful to sample the track interacting with the nucleus with a higher frequency only when needed, that is near the position of ion transversals, where the local dose is rapidly varying). Hence a "tree structure" is created and the data member #v is used to dynamically allocate the memory necessary to contain the grid of pixel just smaller in which each pixel is divided. #numberOfSubPixels stores the number of pixels pointed by #v.
     
        \sa Nucleus_Pixel::createPixels()
     */
    class Pixel
    {
    public:
        
        //! The position of the pixel (\c x coordinate) referred to the nucleus center, expressed in um.
        double x;
        
        //! The position of the pixel (\c y coordinate) referred to the nucleus center, expressed in um.
        double y;
        
        //! The local dose deposited in the pixel, expressed in Gy.
        double dose;
        
        //! The number of pixel constituting the first inner grid.
        /*!
            \sa Nucleus_Pixel::createPixels()
         */
        int numberOfSubPixels;
        
        //! A pointer to the sub-pixels in which the pixel itself is divided (i.e. The pixels of the first subgrid).
        /*!
            \sa Nucleus_Pixel::createPixels()
         */
        Pixel *v;
    };
    
    
    //! Inherited from the Nucleus pure virtual class, it implements the nucleus structure used in the LEM.
    /*!
        \author Andrea Attili
        \author Lorenzo Manganaro
        \author Germano Russo
        \date 2007
     
        It performs the integration on several grids of pixels of varying resolution. Those grids are used to sample with a higher spatial frequency only when needed (i.e. near the position of ion transversals, where the local dose is rapidly varying). Thanks to this approach the single-event survival evaluation is both fast and accurate.
     */
    class Nucleus_Pixel : public Nucleus
    {
    public:
        
        //! Constructor. Instantiates and sets the object.
        /*!
            It divided the nucleus creating the pixel-structure by means of the createPixels() method. It creates four grid of pixels of decreasing dimension.
         
            \param cellLineRef A reference to the corresponding CellLine.
            \param xPosition The nucleus position (\c x coordinate of the center) referred to the beam axis, expressed in mm.
            \param yPosition The nucleus position (\c y coordinate of the center) referred to the beam axis, expressed in mm.
            \param pixelSide1 The side of the smallest (or the third) sub-grid of pixels (#pixelSide_1), expressed in um.
            \param scale1 The scale factor between the second and the third subgrid of pixels (#scale_1).
            \param radius1 The radius of the smallest circumference that defines the sampling of the track, expressed in um.
            \param scale2 The scale factor between the first and the second subgrid of pixels (#scale_2).
            \param radius2 The radius of the second circumference that defines the sampling of the track, expressed in um.
            \param scale3 The scale factor between the biggest grid of pixel and the first sub-grid (#scale_3).
            \param radius3 The radius of the biggest circumference that defines the sampling of the track, expressed in um.
         
            \sa createPixels() and distributeDose()
         */
        Nucleus_Pixel(const CellLine &cellLineRef,
                      const double xPosition = 0.0,
                      const double yPosition = 0.0,
                      const double pixelSide1 = 0.005,
                      const int scale1 = 2,
                      const double radius1 = 0.1,
                      const int scale2 = 10,
                      const double radius2 = 1.0,
                      const int scale3 = 10,
                      const double radius3 = 10.0);
        
        //! Destructor.
        /*!
            Cyclically deletes all pixels created, starting from the smallest grid and ending with the biggest one.
         
            \sa Pixel and createPixels()
         */
        virtual ~Nucleus_Pixel();
        
        //! Adds a constant value of dose absorbed in each pixel of the nucleus.
        /*!
            The method adds systematically a constant value of dose in each pixel my means of a \c for loop over #pixelVector.
         
            \param dose The dose to be added expressed in Gy.
         */
        void addBackgroundDose(const double dose);
        
        //! Performs the sum of the dose absorbed by two different nucleus pixel by pixel, starting from the smallest grid.
        /*!
            This method cycles in some nested \c for loops over the pixels of each subgrid, starting from the smallest one, and adds the doses of the correspondent pixel in the other nucleus.
         
            \warning The function doesn't do anything if the other nucleus is not geometrically similar to this one (e.g. different #r_nucleus of number of pixels - #numberOfSmallestPixels, #numberOfBiggestPixels), but the execution of the program is not terminated.
         
            \param nucleus A reference to another Nucleus_MKM object to evaluate the sum of the doses.
         */
        virtual void addNucleusDoses(Nucleus_Pixel &nucleus);
        
        //! Resets to zero all counters (#inNucleusCount and #intersectionCount) and doses absorbed in the nucleus, pixel by pixel.
        virtual void cleanNucleus();
        
        //! Returns a pointer to a new Nucleus_Pixel object. It not really a clone but a new clean object.
        /*!
            \warning It dynamically allocates memory to be deleted (somewhere) by the user.
         
            \note To create a real clone of another nucleus, a better implementation of the copy constructor is needed.
         */
        virtual Nucleus_Pixel* clone(const CellLine&);
        
        //! Distributes the dose deposited by a Track in the pixels constituting the nucleus.
        /*!
            The sampling of the track is divided into four zone delimited by three concentric circumferences of decreasing radius.
            Starting from the biggest grid, for each pixel it determines if there is intersection with the biggest circle defined for the sampling by calling the intersection(const double, const double, const double, const double, const double, const double, double) method; if the answer is:
                - no: then it calls Track::getLocalDose() to evaluate the dose deposited in that pixel (obviously only if its distance from the track is smaller than the track radius).
                - yes: then it recursively does the same process on the inner grid of pixels.
            In this way, the track is sampled with a higher spatial frequency only when needed (i.e. near the position of ion transversals, where the local dose is rapidly varying); thanks to this approach the single-event survival evaluation is both fast and accurate.
         
            \param track A reference to the track interacting with the nucleus.
         
            \sa createPixels()
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
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and the associated lethal events, with respective uncertainties.
        /*!
            Through four nested \c for loops over the grid of pixels it evaluates the total dose absorbed by the nucleus, considering the dose absorbed by each pixel opportunely normalized, parallel evaluating the total number of lethal events observed, calculated by means of the selected parametrization (CellLine::getLogSurvival_X()).
         
            \param doses The vector of doses absorbed (in Gy), each element refers to a specific pixel, passed by reference to be overwritten.
            \param dosesUncertainty The vector of uncertainties associated to doses absorbed (in Gy), each element refers to a specific pixel, passed by reference to be overwritten.
            \param lethals The vector of lethal events observed, each element refers to a specific pixel, passed by reference to be overwritten.
            \param lethalsUncertainty The vector of uncertainties associated to the number of lethal events observed, each element refers to a specific pixel, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to doses and lethals, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \sa getDoseAndSurvival()
         */
        void getDosesAndLethals(std::vector<double> &doses,
                                std::vector<double> &dosesUncertainty,
                                std::vector<double> &lethals,
                                std::vector<double> &lethalsUncertainty) const;
        
        //! Returns (overwriting parameters passed by reference) the dose absorbed by the nucleus and the associated survival, with respective uncertainties.
        /*!
            Through four nested \c for loops over the grid of pixels it evaluates the total dose absorbed by the nucleus, considering the dose absorbed by each pixel opportunely normalized, parallel evaluating the total number of lethal events observed, calculated by means of the selected parametrization (CellLine::getLogSurvival_X()).
            The cellular survival is then evaluated by means of the following relation:
            \f[
                S=\exp(-L_{TOT})
            \f]
            where \f$L_{TOT}\f$ is the sum of all the lethal events observed.
         
            \param dose The total dose absorbed by the nucleus, expressed in Gy, passed by reference to be overwritten.
            \param doseUncertainty The uncertainty associated to the dose absorbed, expressed in Gy, passed by reference to be overwritten.
            \param survival The cellular survival associated to the dose absorbed by the nucleus, passed by reference to be overwritten.
            \param survivalUncertainty The uncertainty associated to the cellular survival, passed by reference to be overwritten.
         
            \note The method was thought to associate also an uncertainty to dose and survival, but this possibility hasn't been implemented yet, therefore actually -1 is assigned to those values.
         
            \sa getDosesAndLethals()
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
         
            \sa distributeDose(const Track&)
         */
        virtual int getIntersectionCount() const {return intersectionCount;};
        
        //! Returns the number of pixels constituting the biggest grid.
        /*!
            \return The number of pixels constituting the biggest grid (#numberOfBiggestPixels).
         */
        int getNumberOfBiggestPixels() {return numberOfBiggestPixels;};
        
        //! The number of pixels constituting the smallest grid.
        /*!
            \return The number of pixels constituting the biggest grid (#numberOfSmallestPixels).
         */
        int getNumberOfSmallestPixels() {return numberOfSmallestPixels;};
        
        //! Returns the nucleus position (\c x and \c y coordinates) referred to the beam axis and expressed in mm overwriting two \c double variables passed by reference.
        /*!
            This is an unusual getter which needs two \c double variables, passed by reference, that will be overwritten with the \c x and \c y coordinates of the nucleus referred to the beam axis, expressed in mm.
         
            \param returnX The variable to be overwritten with the \c x coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
            \param returnY The variable to be overwritten with the \c y coordinate of the nucleus, expressed in mm, passed by reference to be overwritten.
         
            \sa Track::getPosition()
         */
        virtual void getPosition(double &returnX,
                                 double &returnY) const;
        
        //! Returns the radius of the nucleus, expressed in um.
        /*!
            \return The radius of the nucleus expressed in um.
         */
        virtual double getRadius() const {return r_nucleus;};
        
        //! Saves (writing on a file) the local dose deposited in the smallest pixels and their coordinates.
        /*!
            This method write on a file, for each of the smallest pixel of the nucleus:
                - The \c x and \c y coordinates (in um) identifying its position referred to the beam axis;
                - The dose absorbed (in Gy);
         
            \param fileName The name of the file where to save data.
         */
        void saveLocalDose(const std::string fileName) const;
        
        //! Sets a value of dose in each of the smallest pixels constituting the nucleus, getting it from an external vector.
        /*!
            It runs in a nested \c for loops over the smallest pixels, for each one setting its Pixel::dose (getting it from an external vector)
         
            \param doses A vector containing the doses (in Gy) to be assigned at each of the smallest pixels constituting the nucleus.
         */
        void writeDoses(std::vector<double> &doses);
        
    private:
        
        //! Called by the constructor to create the pixels structure.
        /*!
            The nucleus is divided into four grid of pixels of decreasing dimension. Starting from the biggest grid, the function perform a loop over the pixels and for each pixel it defines the \c x and \c y coordinates with respect to the beam axis (in um). In this way the first grid is created. Then, pixel by pixel, the method verify if there is intersection with the nucleus by means of the intersection(const double, const double, const double) function and, if there is, it proceeds recursively creating a subgrid of pixels inside it in the same way, and so on till the smallest grid.
         */
        void createPixels();
        
        //! Determines if there is intersection between a pixel of the grid and the nucleus.
        /*!
            Evaluates if there is at least one point of the pixel intersecting the circular nucleus looking at the pixel coordinates first, then at the innermost vertex and finally at the innermost point of the edge with respect to the center of the nucleus.
         
            \param x_pixel The \c x coordinate of the pixel referred to the beam axis, expressed in um.
            \param y_pixel The \c x coordinate of the pixel referred to the beam axis, expressed in um.
            \param pixel_side The length of the pixel side, expressed in um.
         
            \return A boolean value to indicate if an intersection occurs.
         */
        inline bool intersection(const double x_pixel,
                                 const double y_pixel,
                                 const double pixel_side) const;
        
        //! Determines if there is intersection between a pixel and a circle with a specified radius centered in the track position.
        /*!
            First identifies the distance between pixel and track, overwriting the distance-parameter. The it evaluates if there is at least one point of the pixel intersecting the circle looking at the distance first, then at the innermost vertex of the pixel and finally at the innermost point of the edge with respect to the position of the track.
         
            \param x_pixel The position of the pixel (\c x coordinate) referred to the beam axis, expressed in um.
            \param y_pixel The position of the pixel (\c y coordinate) referred to the beam axis, expressed in um.
            \param pixel_side The length of the pixel edge, expressed in um.
            \param x_track The position of the track (\c x coordinate) referred to the beam axis, expressed in um.
            \param y_track The position of the track (\c y coordinate) referred to the beam axis, expressed in um.
            \param radius The radius of the circle expressed in um.
            \param distance The distance between pixel and track, expressed in um.
         */
        inline bool intersection(const double x_pixel,
                                 const double y_pixel,
                                 const double pixel_side,
                                 const double x_track,
                                 const double y_track,
                                 const double radius,
                                 double &distance) const;
        
    protected:
        
        //! A reference to a CellLine object where the characteristics of the cell line to which the nucleus belongs are stored.
        const CellLine &cellLine;
        
        //! The radius of the nucleus, expressed in um.
        /*!
            It's instantiated in the constructor getting the value from the CellLine object representing the cell line to which the nucleus belongs.
         
            \sa Nucleus_Pixel()
         */
        double r_nucleus;
        
        //! A pointer to the pixels of the biggest sub-grid.
        Pixel *pixelVector;
        
        //! The number of pixels constituting the biggest grid.
        /*!
            \sa createPixels()
         */
        int numberOfBiggestPixels;
        
        //! The number of pixels constituting the smallest grid.
        /*!
            \sa createPixels()
         */
        int numberOfSmallestPixels;
        
        //! The position of the center of the nucleus (\c x coordinate) referred to the beam axis, expressed in mm.
        const double x_nucleus;
        
        //! The position of the center of the nucleus (\c y coordinate) referred to the beam axis, expressed in mm.
        const double y_nucleus;
        
        //! The side of the smallest (or the third) sub-grid of pixels, expressed in um. Default: 0.005 um.
        const double pixelSide_1;
        
        //! The scale factor between the second and the third subgrid of pixels. Default: 2.
        const int scale_1;
        
        //! The radius of the smallest circumference that defines the sampling of the track, expressed in um.
        /*!
            \sa distributeDose()
         */
        const double radius_1;
        
        //! The side of the second sub-grid of pixels, expressed in um. Default: 0.01 um.
        const double pixelSide_2;
        
        //! The scale factor between the first and the second subgrid of pixels. Default: 10.
        const int scale_2;
        
        //! The radius of the second circumference that defines the sampling of the track, expressed in um.
        /*!
            \sa distributeDose()
         */
        const double radius_2;
        
        //! The side of the first sub-grid of pixels, expressed in um. Default: 0.1 um.
        const double pixelSide_3;
        
        //! The scale factor between the biggest grid of pixel and the first sub-grid. Default: 10.
        const int scale_3;
        
        //! The radius of the biggest circumference that defines the sampling of the track, expressed in um.
        /*!
            \sa distributeDose()
         */
        const double radius_3;
        
        //! The side of the biggest grid of pixels, expressed in um. Default: 1 um.
        const double pixelSide_4;
        
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


#endif /* NUCLEUS_PIXEL_H */
