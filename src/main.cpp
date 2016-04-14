#include "Calculus.h"
#include "Nucleus_MKM.h"
#include "Nucleus_tMKM.h"
#include "Nucleus_Pixel.h"
#include "CellLine.h"
#include "Tracks.h"
#include "Particles.h"
#include "Particle.h"
#include "usefulFunctions.h"

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::ostringstream;

#include <iomanip>
using std::setprecision;
using std::scientific;

#include <cmath>
using std::exp;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::clog;
using std::fixed;
using std::ios;

using std::size_t;

using std::string;
using std::vector;

using namespace Survival;

/*!
    \mainpage
 
    \author Andrea Attili
    \author Lorenzo Manganaro
    \author Germano Russo
 
    \date 2015
 
    \section usage Usage
    \subsection Linux_MAC Unix and Unix-like systems.
    To execute the program the user has to call from the command line:
    \code{.sh}
    $ source setenv.sh
    $ ./survival -SIMULATION_OPTION CHOSEN_VALUES ...
    \endcode
 
    Then the user has the possibility to set a number of physical (and not only physical) parameters by using the syntax: \c -PARAMETER_NAME \c PARAMETER_VALUE
 
    Here is the complete list of parameters and their meaning:
        - \c -projectName It's a string representing the prefix to give at any file and directories that will be created in the simulation. The default value is "NewProject".
        - \c -output The user has the possibility to choose between three kinds of output (and all possible combination between them):
            -# "LQ_pars" Then a file will be created, named "PROJECTNAME_LQparameters_MKM.csv", containing the information about the parameters chosen for the simulation and the values of the simulated LQ \f$\alpha\f$ and \f$\beta\f$ parameters (a new line for each energy evaluated).
            -# "meanValues" Then a file will be created, named "PROJECTNAME_survival_MKM.csv", containing the information about the parameters chosen for the simulation and the values of doses delivered and survival observed (a new line for each energy or dose evaluated).
            -# "cellValues" This kind of output is supported only by the MonteCarlo calculusType. It is a way to store the values of dose and survival obtained for each single cell irradiated during the monte carlo simulation. Then a directory will be created, named "PROJECTNAME_survival_data". In the directory the user will find a description file named "000_MonteCarlo_parameters.csv", listing the parameters used in the simulation, and a directory with the same name containing the corresponding data. In particular in the subdirectory some file will be created (a file for each level of dose imposed), each one containing two column with the dose delivered and the survival observed for each cell irradiated. When a new simulation is lauched with the same project name, the program will do a check over all the description files in the directory, if the parameters of the simulation are the same of another one already done, then it will enter the related subdirectory and append data there, if not a new description file (with progressive number) and corresponding subdirectory will be created.
            \warning Storing these informations could occupy a lot of memory in the computer, depending on the parameter set for the simulation. Use with caution!
            \note This parameter has to be specified. No default values are set.
        - \c -precision Supported only by the MonteCarlo calculusType, it's a \c double identifying the precision to be reached in the calculation. Two possibilities are provided, as the user can indicate:
            -# A positive integer: it will be taken as the number of cell to irradiate for each level of nominal dose imposed.
            -# A double between 0 and 1: it indicates the relative error on the simulated survival to be reached before stopping the calculation.
        - \c -parallelismType A positive integer indicating the level of parallelism to be used in the simulation. The user has the possibility to specify the number of threads to dedicate at the calculation, in particular:
            -# 1: Only 1 thread = parallelism disabled.
            -# N>1: The number of threads to be used.
            -# 0: Then the program will define a number of thread corresponding to the number of core of the computer executing the program. This is set as default value.
        - \c -model It's a string representing the model to use in the simulation. Five models are supported:
            -# "LEMI" Is the first formulation of the Local Effect Model, the published reference is:\n
                M. Scholz and G. Kraft, "Track structure and the calculation of biological effects of heavy charged particles", \a Advances \a in \a Space \a Research \b 18, 5-14 (1996)
            -# "LEMII" A reformulation of the Local Effect Model, as described in:\n
                T. Elsässer and M. Scholz, "Cluster effects within the local effect model", \a Radiation \a Research \b 167, 319-329 (2007)
            -# "LEMIII" A third formulation of the Local Effect Model, well described in the paper:\n
                T. Elsässer, M. Krämer and M. Scholz, "Accuracy of the local effect model for the prediction of biologic effects of carbon ion beams \em in \em vitro and \em in \em vivo", \a International \a Journal \a of \a Radiation \a Oncology-\a Biology-\a Physics \b 71, 866-872 (2008)
            -# "MKM" The Microdosimetric Kinetic Model, in the formulation of Hawkins, with the approach of Kase who suggest to use the Kiefer-Chatterjee amorphous track model. Some published references:\n
                R.B. Hawkins, "A Statistical Theory of Cell Killing by Radiation of Varying Linear Energy Transfer", \a Radiation \a Research \b 140, 366-374 (1994) -- And subsequent references\n
                Y. Kase, T. Kanai, N. Matsufuji, Y. Furusawa, T. Elsasser, and M. Scholz, "Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation", \a Physics \a in \a Medicine \a and \a Biology \b 53, 37-59 (2008)
            -# "tMKM" A monte carlo reformulation of the MKM, extended to the temporal dimension to evaluate the effect of the time structure of the irradiation on the LQ parameters, as described in:\n
                L. Manganaro, G. Russo, R. Cirio, F. Dalmasso, S. Giordanengo, V. Monaco, R. Sacchi, A. Vignati, A. Attili, "A novel formulation of the Microdosimetric Kinetic Model to account for dose-delivery time structure effects in ion beam therapy with application in treatment planning simulations", \a Medical \a Physics, (Submitted)\n
                The default value for this option is "MKM"
        - \c -calculusType A string identifying the type of calculus to be done. Some possibilities are available:
                -# "rapidScholz" It's an implementation of the method described in:\n
                M. Krämer and M. Scholz, "Rapid calculation of biological effects in ion radiotherapy", \a Physics \a in \a medicine \a and \a biology \b 51, 1959-1970 (2006)\n
                It's compatible only with LEMI-LEMII-LEMIII models.
                -# "rapidRusso" A new rapid method for LEM, proved to be more accurate, described in:\n
                G. Russo, "Develpment of a radiobiological database for carbon ion Treatment Planning Systems - Modelling and simulating the irradiation process", \a PhD \a Thesis, Università degli studi di Torino (2011)\n
                Also this method is compatible only with LEMI-LEMII-LEMIII models.
                -# "rapidMKM" An implementation of the original MKM calculation, described in:\n
                R.B. Hawkins, "A Statistical Theory of Cell Killing by Radiation of Varying Linear Energy Transfer", \a Radiation \a Research \b 140, 366-374 (1994) -- And subsequent references.\n
                It's compatible only with the MKM model and it's the default value for this option.
                -# "MonteCarlo" Compatible with all models implemented, performs a monte carlo simulation of the irradiation process to get the LQ parameters.
        - \c -cellType A string identifying the name of the cell lline used in the calculation. The default value is "Cell1"
                \note The cell line in reality is completely determined by the model parameters chosen. This is only a tag to indicate the cell but it isn't used in the simulation.
        - \c -MKM_alpha0 A double representing IDEALLY the linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$. The simulated \f$\alpha\f$ parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.312 \f$Gy^{-1}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -MKM_beta0 A double representing IDEALLY the linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$. The simulated \f$\beta\f$ parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.073 \f$Gy^{-2}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -MKM_rNucleus The radius of the cell expressed in \f$\mu m\f$. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 4.611 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -MKM_rDomain The radius of domains wich constitute the MKM nucleus expressed in \f$\mu m\f$. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.365 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -MKM_timeConst The time constant associated to the repair kinetics of the nucleus.
        - \c -LEM_alpha0 A double representing IDEALLY the linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$. The simulated \f$\alpha\f$ parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.312 \f$Gy^{-1}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -LEM_beta0 A double representing IDEALLY the linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$. The simulated \f$\beta\f$ parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.073 \f$Gy^{-2}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -LEM_rNucleus The radius of the cell expressed in \f$\mu m\f$. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 4.611 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.
        - \c -LEM_Dt The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 30 Gy.
        - \c -ion A string identifying the chemical symbol of the element, without mass number specifications. Ions from proton to neon are supported.
        - \c -energies A sequence of kinetic energies of the primary ions to be evaluated, expressed in MeV. This and the "-lets" option are mutually exclusive, but one of the two has to be specified.
        - \c -lets A sequence of LET of the primary ions to be evaluated, expressed in MeV. This and the "-energies" option are mutually exclusive, but one of the two has to be specified.
        - \c -doses The sequence of doses to be delivered. The default is 1 to 6 Gy (step 1), in order to construct a survival curve.
        - \c -nFraction Supported only by the tMKM model, it require an integer representing the number of fraction in which to divide each nominal dose to be delivered. The default is 1 (single fraction).
        - \c -timeSpacing Supported only by the tMKM model, it indicate the time spacing between consecutive fractions, expressed in hours. The default is 0 (no time spacing).
        - \c -fracDeliveryTime Supported only by the tMKM model, it indicate the fraction delivery time, expressed in hours. The default is 0 (istantaneous delivering).
        - \c -spectrum_file Useful for mixed fields evaluation. This option provide the possibility to irradiate the cellular population with different ions and different energies in a mixed field, described in an external file (see the TEMPLATE for the spectrum_file). The user has to specify the name of the file containing the spectrum complete with its relative or absolute path.
            \warning The "mixed fields option" hasn't been deeply tested!
        - \c -trackMode In the case of mixed fields this option require a string indicating the way in which to interpretate the spectrum specification. Two possibilities are provided:
            -# "histogram": then the spectrum will be considered as an histogram where each particle is a bin with its "weight" (see the TEMPLATE for the spectrum_file).
            -# "random": then the program will random extract with uniform probability (iteration by iteration) the particle to use.
 
    Typing \c --help a hint will be display, suggesting how to use the program.
 */
int main(int argc, char* argv[])
{
    // CONVERSION CONSTANT
    const double AMU2MEV = 931.494027;
    
    // SETTING DEFAULT VALUES:
    string projectName = "NewProject";
    string cellType = "Cell1";
    string model = "MKM";
    string trackType = "KieferChatterjee";
    string parametrizationType = "LQ_noDt";
    string calculusType = "rapidMKM";
    
    double precision = 0.05;
    int parallelismType = 0;
    
    vector<double> doses;
    
    vector<string> parameter_name;
    
    double MKM_alpha0 = 0.312;
    double MKM_beta0 = 0.073;
    double MKM_rNucleus = 4.611;
    double MKM_rDomain = 0.365;
    
    double tMKM_ac = 2.187;
    
    double LEM_alphaX = 0.312;
    double LEM_betaX = 0.073;
    double LEM_rNucleus = 4.611;
    double LEM_Dt = 30;
    
    string ionType = "C";
    int particleA = 12;
    int particleZ = 6;
    string trackMode = "histogram";
    string energyType = "energy";
    vector<double> energies;
    
    int nFraction = 1;
    double timeSpacing = 0.0;
    double fracDeliveryTime = 0.0;
    
    bool saveAlphaBeta = false;
    bool saveMeans = false;
    bool saveCell = false;
    
    bool mono = true;
    string spectrum_file="";
    
    // PARSING INPUT ARGUMENTS:
    parse(argc, argv, cellType, model, trackType, parametrizationType,
          calculusType, precision, parallelismType, doses, parameter_name,
          MKM_alpha0, MKM_beta0, MKM_rNucleus, MKM_rDomain, tMKM_ac,
          LEM_alphaX, LEM_betaX, LEM_rNucleus, LEM_Dt,
          ionType, particleA, particleZ, trackMode, energyType, energies,
          nFraction, timeSpacing, fracDeliveryTime,
          saveAlphaBeta, saveMeans, saveCell, projectName,
          mono, spectrum_file);
    
    // FILE TO SAVE ALPHA-BETA VALUES
    ostringstream filename_LQ;
    ofstream LQ_File;
    if (saveAlphaBeta) {
        filename_LQ << projectName << "_LQparameters_" << model << ".csv";
        
        bool already_exist = false;
        if(ifstream(filename_LQ.str().c_str()))
            already_exist = true;
        
        LQ_File.open(filename_LQ.str().c_str(), ios::app);
        if(!LQ_File)
        {
            cerr << "File " << filename_LQ.str() << " could not be opened." << endl;
            exit(1);
        }
        if(!already_exist)
        {
            LQ_File << "model,calculusType,cell,";
            for (size_t i=0; i<parameter_name.size(); i++)
                LQ_File << parameter_name[i] << ",";
            LQ_File << "radiation_type,particle,meanEnergy,meanEnergy_sigma,meanLET,meanLET_sigma,LETd,LETd_sigma,"
                    << "nFractions,timeSpacing,fracDeliveryTime,alpha,beta,alpha_sigma,beta_sigma" << endl;
        }
    }
    else
        filename_LQ << "noSave";
    
    // FILE TO SAVE MEAN DOSE AND SURVIVAL VALUES
    string title_means;
    ofstream means_File;
    if (saveMeans) {
        title_means = projectName + "_survival_" + model + ".csv";
        ostringstream filename_means;
        filename_means << title_means;
        
        bool already_exist = false;
        if (ifstream(filename_means.str().c_str()))
            already_exist = true;
        
        means_File.open(filename_means.str().c_str(), ios::app);
        if(!means_File)
        {
            cerr << "File " << filename_means.str() << " could not be opened." << endl;
            exit(1);
        }
        if (!already_exist)
        {
            means_File << "model,calculusType,cell,";
            for (size_t i=0; i<parameter_name.size(); i++)
                means_File << parameter_name[i] << ",";
            means_File << "radiation_type,particle,meanEnergy,meanEnergy_sigma,meanLET,meanLET_sigma,LETd,LETd_sigma,"
                       << "nFractions,timeSpacing,fracDeliveryTime,nominal_dose,"
                       << "mean_dose,mean_survival,mean_doseUncertainty,mean_survivalUncertainty" << endl;
        }
    }
    else
        title_means = "noSave";
    
    // FILE TO SAVE SINGLE CELL DOSE AND SURVIVAL VALUES
    string outputDir;
    if (saveCell) {
        outputDir = projectName + "_survival_data/";
        mkdir(outputDir.c_str());
    }
    
    // SETTING THE CELL LINE (USING THE INPUT PARAMETERS)
    CellLine *cellLine = new CellLine(cellType);
    if (model == "LEMI" || model == "LEMII" || model == "LEMIII")
    {
        cellLine->setNucleusRadius(LEM_rNucleus);
        cellLine->setDomainRadius(LEM_rNucleus); // not used...
        if (parametrizationType == "LQ") cellLine->addParametrization_LQ(LEM_alphaX, LEM_betaX, LEM_Dt);
        if (parametrizationType == "LQ2") cellLine->addParametrization_LQ2(LEM_alphaX, LEM_betaX, LEM_Dt, 0);
        if (parametrizationType == "LQ3") cellLine->addParametrization_LQ3(LEM_alphaX, LEM_betaX, LEM_Dt, 0);
    }
    else if (model == "MKM")
    {
        cellLine->setNucleusRadius(MKM_rNucleus);
        cellLine->setDomainRadius(MKM_rDomain);
        cellLine->addParametrization_LQ_noDt(MKM_alpha0, MKM_beta0);
    }
    else if (model == "tMKM")
    {
        cellLine->setNucleusRadius(MKM_rNucleus);
        cellLine->setDomainRadius(MKM_rDomain);
        cellLine->addParametrization_LQ_noDt_T(MKM_alpha0, MKM_beta0, tMKM_ac);
    }
    cellLine->setParametrization(parametrizationType);
                    
    // CREATING THE NUCLEUS
    Nucleus_Pixel nucleus(*cellLine);
    cout << "Nucleus of type PIXEL created." << endl;
    Nucleus_MKM nucleus_MKM(*cellLine);
    cout << "Nucleus of type MKM created." << endl;
    Nucleus_tMKM nucleus_tMKM(*cellLine);
    cout << "Nucleus of type tMKM created." << endl;
    
    // CREATING THE PARTICLE
    Particle particle;
    particle.type = ionType;
    particle.charge = particleZ;
    particle.A = particleA;
    particle.restEnergy = particle.A * AMU2MEV;  // MeV
    particle.x = 0.0; // mm
    particle.y = 0.0; // mm
    particle.z = 0.0; // mm
    particle.weight = 1.0;
    
    // MODEL PARAMETERS
    vector<double> parameters;
    if (model == "MKM" || model == "tMKM") {
        parameters.push_back(MKM_alpha0);
        parameters.push_back(MKM_beta0);
        parameters.push_back(MKM_rNucleus);
        parameters.push_back(MKM_rDomain);
        if (model == "tMKM")
            parameters.push_back(tMKM_ac);
    }
    else {
        parameters.push_back(LEM_alphaX);
        parameters.push_back(LEM_betaX);
        parameters.push_back(LEM_rNucleus);
        parameters.push_back(LEM_Dt);
    }
    
    // ALPHA, BETA VECTORS
    int NUMBER_OF_ITERATIONS;
    if (mono)
        NUMBER_OF_ITERATIONS = energies.size();
    else
        NUMBER_OF_ITERATIONS = 1;
    vector<double> alpha(NUMBER_OF_ITERATIONS, 0.0);
    vector<double> beta(NUMBER_OF_ITERATIONS, 0.0);
    vector<double> alpha_sigma(NUMBER_OF_ITERATIONS, 0.0);
    vector<double> beta_sigma(NUMBER_OF_ITERATIONS, 0.0);
                    
    // LOOP OVER LETS/ENERGIES
    for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) {
        
        Particles particles;
        
        if (mono) {
            if (energyType == "energy")
            {
                particle.e_c = energies[i]; // MeV
                particle.let = betheBloch_Srim(particle.type, particle.e_c * AMU2MEV / particle.restEnergy); // MeV/um
            }
            else if (energyType == "LET")
            {
                particle.let = energies[i] / 1000.; // MeV/um (si suppone venga data in keV/um)
                particle.e_c = betheBloch_inv_Srim(particle.type, particle.let) * particle.restEnergy / AMU2MEV; // MeV
            }
            
            clog << endl << "Setting particle kinetic energy: " << particle.e_c << " MeV  —  ("
                         << particle.e_c/(double)particle.A << " MeV/u)" << endl
                         << "                        and LET: " << particle.let * 1000. << " keV/um" << endl
                         << endl;
            
            particles << particle;
        }
        else {
            Particles particles_tmp(spectrum_file);
            particles = particles_tmp.getIons();
            particles.setSpectrumFile(spectrum_file);
            particles.reconstructIonLETandEnergy();
        }
        
        Tracks tracks(particles, trackType);

        // CREATING AND EXECUTING CALCULUS
        if (calculusType == "rapidScholz")
        {
            Calculus calculus(tracks, *cellLine, nucleus, filename_LQ.str(), model, parallelismType);
            calculus.rapidScholz_alphaIon_betaIon(alpha[i], beta[i]);
        }
        else if (calculusType == "rapidRusso")
        {
            Calculus calculus(tracks, *cellLine, nucleus, filename_LQ.str(), model, parallelismType);
            calculus.rapidRusso_alphaIon_betaIon(alpha[i], beta[i]);
        }
        else if (calculusType == "rapidMKM")
        {
            Calculus calculus(tracks, *cellLine, nucleus_MKM, filename_LQ.str(), model, parallelismType);
            calculus.rapidMKM_Kase_alphaIon_betaIon(alpha[i], beta[i]);
        }
        else if (calculusType == "MonteCarlo")
        {
            ostringstream prefix;
            if (saveCell)
            {
                int iter = 0;
                while (true) {
                    ostringstream find_file_path;
                    find_file_path << "./" << outputDir << (iter<10 ? "00" : (iter<100 ? "0" : "")) << iter << "_MonteCarlo_parameters.csv";
                    if (ifstream(find_file_path.str().c_str()))
                    {
                        ifstream read_file(find_file_path.str(), ifstream::in);
                        if (!read_file) {
                            cerr << "File " << find_file_path.str() << " could not be opened." << endl;
                            exit(1);
                        }
                        string read_line, read_col2;
                        int nLine = 1;
                        int incr = 0;
                        bool isEqual = true;
                        while (!read_file.eof()) {
                            getline(read_file, read_line, '\n');
                            read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                            if (nLine==2 && read_col2!=model)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==3 && read_col2!=parametrizationType)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==4 && read_col2!=calculusType)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==5 && read_col2!=cellType)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==6)
                            {
                                for (size_t a=0; a<parameters.size(); ++a) {
                                    ostringstream cfr;
                                    cfr << setprecision(3) << fixed << parameters[a];
                                    if(read_col2!=cfr.str()) {
                                        isEqual = false;
                                        break;
                                    }
                                    nLine++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                }
                                if (!isEqual)
                                    break;
                            }
                            if (nLine==10) {
                                if (mono) {
                                    if (read_col2!="Monoenergetic") {
                                        isEqual = false;
                                        break;
                                    }
                                    incr++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    if (read_col2!=ionType) {
                                        isEqual = false;
                                        break;
                                    }
                                    incr+=3;
                                    getline(read_file, read_line, '\n');
                                    getline(read_file, read_line, '\n');
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    nLine += incr;
                                }
                                else {
                                    if (read_col2!="Spectrum") {
                                        isEqual = false;
                                        break;
                                    }
                                    incr++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    if (read_col2!=spectrum_file) {
                                        isEqual = false;
                                        break;
                                    }
                                    incr++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    ostringstream cfr1;
                                    cfr1 << scientific << tracks.getDoseAveragedLet();
                                    if (read_col2!=cfr1.str()) {
                                        isEqual = false;
                                        break;
                                    }
                                    incr++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    ostringstream cfr2;
                                    cfr2 << scientific << tracks.getSigmaDoseAveragedLet();
                                    if (read_col2!=cfr2.str()) {
                                        isEqual = false;
                                        break;
                                    }
                                    incr++;
                                    getline(read_file, read_line, '\n');
                                    read_col2 = read_line.substr(read_line.find(",")+1,read_line.find("\n"));
                                    nLine += incr;
                                }
                            }
                            if (nLine==(10+incr) && read_col2!=trackType)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==(11+incr) && read_col2!=trackMode)
                            {
                                isEqual = false;
                                break;
                            }
                            if (nLine==(12+incr))
                            {
                                ostringstream cfr;
                                cfr << nFraction;
                                if(read_col2!=cfr.str()) {
                                    isEqual = false;
                                    break;
                                }
                            }
                            if (nLine==(13+incr))
                            {
                                ostringstream cfr;
                                cfr << setprecision(3) << fixed << timeSpacing;
                                if(read_col2!=cfr.str()) {
                                    isEqual = false;
                                    break;
                                }
                            }
                            if (nLine==(14+incr))
                            {
                                ostringstream cfr;
                                cfr << setprecision(3) << fixed << fracDeliveryTime;
                                if(read_col2!=cfr.str()) {
                                    isEqual = false;
                                    break;
                                }
                            }
                            nLine++;
                        }// while eof()
                        if (isEqual)
                        {
                            prefix << outputDir << (iter<10 ? "00" : (iter<100 ? "0" : "")) << iter
                                   << "_MonteCarlo_survival_data/" << "dose_survival_";
                            if (mono)
                                prefix << setprecision(2) << fixed << energies[i] << (energyType=="LET" ? "(keV-um)_" : "(MeV-u)_");
                            break;
                        }
                        else
                            iter++;
                    } // find
                    else {
                        ofstream Info_File(find_file_path.str().c_str());
                        if (!Info_File) {
                            cerr << "File " << find_file_path.str() << " could not be opened." << endl;
                            exit(1);
                        }
                        Info_File << "PARAMETER_NAME,VALUE" << endl
                                << "Model," << model << endl
                                << "Parametrization," << parametrizationType << endl
                                << "Calculus," << calculusType << endl
                                << "Cell," << cellType << endl << setprecision(3) << fixed;
                        for (size_t a=0; a<parameters.size(); a++)
                            Info_File << parameter_name[a] << ","  << parameters[a] << endl;
                        Info_File << "Radiation_type,";
                        if (mono) {
                            Info_File << "Monoenergetic" << endl
                                      << "Particle," << ionType << endl
                                      << "ParticleA," << particleA << endl
                                      << "ParticleZ," << particleZ << endl;
                        }
                        else {
                            Info_File << "Spectrum" << endl
                                      << "File," << spectrum_file << endl
                                      << "LETd," << tracks.getDoseAveragedLet() << endl
                                      << "LETd_sigma," << tracks.getSigmaDoseAveragedLet() << endl;
                        }
                        
                        Info_File
                                << "Track_Type," << trackType << endl
                                << "Track_Mode," << trackMode << endl
                                << "Number_of_Fractions," << nFraction << endl
                                << "Time_Spacing," << timeSpacing << endl
                                << "Fraction_Delivery_Time," << fracDeliveryTime;
                        Info_File.close();

                        prefix << outputDir << (iter<10 ? "00" : (iter<100 ? "0" : "")) << iter << "_MonteCarlo_survival_data/";
                        mkdir(prefix.str().c_str());
                        prefix << "dose_survival_";
                        if (mono)
                            prefix << setprecision(2) << fixed << energies[i] << (energyType=="LET" ? "(keV-um)_" : "(MeV)_");
                        break;
                    }
                }// while iter
            }
            else
                prefix << "noSave";
            
            if (model == "LEMI" || model == "LEMII" || model == "LEMIII")
            {
                Calculus calculus_LEM(tracks, *cellLine, nucleus, prefix.str(), model, parallelismType);
                calculus_LEM.slow_alphaIon_betaIon(trackMode, parameters, doses, precision,
                                                   alpha[i], alpha_sigma[i], beta[i], beta_sigma[i],
                                                   nFraction, timeSpacing, fracDeliveryTime,
                                                   saveAlphaBeta, saveMeans, saveCell, title_means);
            }
            else if (model == "MKM")
            {
                Calculus calculus_MKM(tracks, *cellLine, nucleus_MKM, prefix.str(), model, parallelismType);
                calculus_MKM.slow_alphaIon_betaIon(trackMode, parameters, doses, precision,
                                                   alpha[i], alpha_sigma[i], beta[i], beta_sigma[i],
                                                   nFraction, timeSpacing, fracDeliveryTime,
                                                   saveAlphaBeta, saveMeans, saveCell, title_means);
            }
            else if (model == "tMKM")
            {
                Calculus calculus_tMKM(tracks, *cellLine, nucleus_tMKM, prefix.str(), model, parallelismType);
                calculus_tMKM.slow_alphaIon_betaIon(trackMode, parameters, doses, precision,
                                                    alpha[i], alpha_sigma[i], beta[i], beta_sigma[i],
                                                    nFraction, timeSpacing, fracDeliveryTime,
                                                    saveAlphaBeta, saveMeans, saveCell, title_means);
            }
        }//if MonteCarlo
        
        if(calculusType!="MonteCarlo" && saveMeans) {
            double meanEn = tracks.getMeanEnergy();
            double meanEnSig = tracks.getSigmaMeanEnergy();
            double meanL = tracks.getMeanLet();
            double meanLSig = tracks.getSigmaMeanLet();
            double meanDoseAvL = tracks.getDoseAveragedLet();
            double meanDoseAvLSig = tracks.getSigmaDoseAveragedLet();
            
            for (size_t d=0; d<doses.size(); ++d) {
                double S = exp(-alpha[i]*doses[d]-beta[i]*doses[d]*doses[d]);
            
                means_File << model << "," << calculusType << "," << cellType << ",";
                for (size_t a=0; a<parameters.size(); a++)
                    means_File << parameters[a] << ",";
                if (mono)
                    means_File << "Monoenergetic," << ionType;
                else
                    means_File << ((meanEnSig==0 && meanLSig==0 && meanDoseAvLSig==0) ? "Monoenergetic," : "Spectrum,") << spectrum_file;
                means_File << "," << meanEn << "," << (meanEnSig<1e-12 ? 0 : meanEnSig) << ","
                           << meanL << "," << (meanLSig<1e-12 ? 0 : meanLSig) << ","
                           << meanDoseAvL << "," << (meanDoseAvLSig<1e-12 ? 0 : meanDoseAvLSig) << ","
                           << nFraction << "," << timeSpacing << "," << fracDeliveryTime << "," << doses[d] << ",";
                means_File << scientific << doses[d] << "," << S << "," << 0 << "," << 0 << endl;
                means_File.unsetf(std::ios::scientific);
            }
        }
        
        clog << "------------------------------------------------------" << endl;
        
        if (saveAlphaBeta) {
            double mEsig = tracks.getSigmaMeanEnergy();
            double mLsig = tracks.getSigmaMeanLet();
            double mDALsig = tracks.getSigmaDoseAveragedLet();
            LQ_File << model << "," << calculusType << "," << cellType << ",";
            for (size_t a=0; a<parameters.size(); a++)
                LQ_File << parameters[a] << ",";
            if (mono)
                LQ_File << "Monoenergetic," << ionType;
            else
                LQ_File << ((mEsig==0 && mLsig==0 && mDALsig==0) ? "Monoenergetic," : "Spectrum,") << spectrum_file;
            LQ_File << "," << tracks.getMeanEnergy() << "," << (mEsig<1e-12 ? 0 : mEsig) << ","
                    << tracks.getMeanLet() << "," << (mLsig<1e-12 ? 0: mLsig) << ","
                    << tracks.getDoseAveragedLet() << "," << (mDALsig<1e-12 ? 0 : mDALsig) << ","
                    << nFraction << "," << timeSpacing << "," << fracDeliveryTime << "," << alpha[i] << "," << beta[i] << ","
                    << alpha_sigma[i] << "," << beta_sigma[i] << endl;
        }
    }
    LQ_File.close();
         
    delete cellLine;
    
    return 0;
}
