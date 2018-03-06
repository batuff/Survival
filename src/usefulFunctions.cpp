#include "usefulFunctions.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_sf_gamma.h>

#include <map>
using std::map;
using std::pair;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::sqrt;

#include <fstream>
using std::ifstream;

#include <iostream>
using std::cerr;
using std::ios;
using std::endl;
using std::cout;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;
using std::getline;

#include <sstream>
using std::stringstream;

#include <sys/stat.h>

// for windows mkdir
#ifdef WIN32
#include <direct.h>
#endif

namespace Survival {
    
int _mkdir(const char* path)
{
    #ifdef WIN32
        return ::_mkdir(path);
    #else
    #if LINUX
        return ::mkdir(path, S_IRWXU | S_IRWXG | S_IRWXO);
    #else
    #if OSX
        return ::mkdir(path, 0755);
    #else
        cerr << "Error: Unknown Operating System" << endl;
        exit(1);
    #endif
    #endif
    #endif
}
    
//------------------------------------------------------------------------------

//! Class used to load precomputed values of LET for different ions of varying kinetic energy, evaluated by means of the Bethe-Bloch formula.
/*!
    \author Andrea Attili
    \author Lorenzo Manganaro
    \author Germano Russo
    \date 2008--2011
 
    Depending on the ion type chosen, it loads (in the constructor) the correspondent file from the "data" directory containing precomputed values of LET (linear energy transfer) each one associated to a different kinetic energy, storing them into two vectors (one for the kinetic energies #e_c and one for the LETs #let). It provides a method (GetLET()) to get the LET of the ions correspondent to a given kinetic energy.
 */
class BetheBlochTable
{
public:

    //! Constructor. Instantiates and sets the object.
    /*!
        Opens the data file containing the table correspondent to the chosen ion and store the values of kinetic energies and LETs in #e_c and #let respectively.
     
        \warning The execution of the program will be terminated if the data file doesn't exist or if the environment hasn't been set correctly.
     */
    BetheBlochTable(const string &ionType_)
    : ionType( ionType_ )
    {
        char *dataPath;
        dataPath = getenv("DATA");
        if( dataPath == NULL )
        {
            cerr << "DATA environment variable could not be determined." << endl;
            exit(1);
        }

        string betheBlochFile(dataPath);
        betheBlochFile += "/" + ionType + "_Srim.dat";
        ifstream inBetheBlochFile(betheBlochFile.c_str(), ios::in);
        if( !inBetheBlochFile )
        {
            cerr << "File " << betheBlochFile << " could not be opened." << endl;
            exit(1);
        }

        char tmp[100];
        for( int i = 1; i <= 4; i++ )
            inBetheBlochFile.getline(tmp, 100);

        double tmp_e_c, let_electrons, let_nuclear, value;
        while( inBetheBlochFile >> tmp_e_c >> let_electrons >> let_nuclear >> value >> value >> value )
        {
//            let = let_electrons + let_nuclear;  // MeV/um
            e_c.push_back( tmp_e_c );  // MeV/u
            let.push_back( let_electrons );  // MeV/um
        }

        inBetheBlochFile.close();
    };

    //! Returns the value of LET corresponding to the kinetic energy imposed.
    /*!
        Find the indexes in the #e_c vector of the two nearest neighbors of the energy fixed and interpolates them to get the correct value of LET.
     
        \param e_c_imposed The kinetic energy (in MeV/u) correspondent to the required LET.
     
        \return The value of LET corresponding to the kinetic energy imposed, expressed in MeV/um.
     
        \warning The execution of the program will be terminated if the kinetic energy imposed is out of the available range for the ion specified.
     
        \sa betheBloch_Srim()
     */
    double GetLET(double e_c_imposed) const
    {
        int i_1 = 0;
        int i_2 = e_c.size() - 1;

        //Check that the value is inside the range spanned by the array.
        if( e_c_imposed < e_c[0] || e_c.back() < e_c_imposed )
        {
            cerr << "betheBlochSrim(): kinetic energy (" << e_c_imposed << " MeV/u) out of available [" << e_c[0] << ", " << e_c.back() << "] range for " << ionType << " particles." << endl;
            exit(1);
        }

        //Find bounding indices.
        int index;
        while( i_2 - i_1 > 1 )
        {
            index = (i_1 + i_2) >> 1;  //equivalent to floor( (i_1 + i_2)/2 )
            if( e_c_imposed > e_c[index] )
                i_1 = index;
            else
                i_2 = index;
        }

        //Compute weights.
        double e_cWeight_1 = (e_c[i_2] != e_c[i_1]) ? (e_c[i_2] - e_c_imposed) / (e_c[i_2] - e_c[i_1]) : 1.;
        double e_cWeight_2 = 1 - e_cWeight_1;

        //Linear interpolation.
        return let[i_1] * e_cWeight_1 + let[i_2] * e_cWeight_2;
    };

private:
    
    //! A string defining the name of the ion (chemical symbol).
    string ionType;
    
    //! A vector to store the list of kinetic energies, expressed in MeV/u.
    vector<double> e_c;
    
    //! A vector to store the precomputed LETs, expressed in MeV/um.
    vector<double> let;
};

//------------------------------------------------------------------------------

double betheBloch_inv_Srim(string ionType, double let_imposed)  // MeV/u
{
    char *dataPath;
    dataPath = getenv("DATA");
    if( dataPath == NULL )
    {
        cerr << "DATA environment variable could not be determined." << endl;
        exit(1);
    }
    
    string betheBlochFile(dataPath);
    betheBlochFile += "/" + ionType + "_Srim.dat";
    ifstream inBetheBlochFile(betheBlochFile.c_str(), ios::in);
    if( !inBetheBlochFile )
    {
        cerr << "File " << betheBlochFile << " could not be opened." << endl;
        exit(1);
    }
    
    char tmp[100];
    
    for( int i = 1; i <= 5; i++ )
        inBetheBlochFile.getline(tmp, 100);
    
    double value;
    double let = 0., let_save = 0.;
    double e_c, e_c_save = 0.;
    double let_electrons, let_nuclear;
    
    bool isValid = false;
    while( inBetheBlochFile >> e_c >> let_electrons >> let_nuclear >> value >> value >> value )
    {
        //       let = let_electrons + let_nuclear; // MeV/um
        let = let_electrons; // MeV/um
        
        if( let <= let_imposed && let_save >= let_imposed ) // deve aver superato il massimo del let (si assumono energie cinetiche non troppo piccole...)
        {
            isValid = true;
            break;
        }
        
        e_c_save = e_c;
        let_save = let;
    }
    
    inBetheBlochFile.close();
    
    if( let_save != 0.0 && isValid )
        return e_c_save + (e_c - e_c_save) * (let_imposed - let_save) / (let - let_save);
    else
    {
        cerr << "betheBlochSrim(): let out of available data" << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

double betheBloch_Srim(string ionType, double e_c_imposed)  // MeV/u
{
    static map<string, BetheBlochTable> ionToTable;

    map<string, BetheBlochTable>::iterator iter = ionToTable.find(ionType);

    if( iter == ionToTable.end() )
    {
        pair<string, BetheBlochTable> newElement(ionType, BetheBlochTable(ionType));
        pair<map<string, BetheBlochTable>::iterator, bool> returnFromInsert;
        returnFromInsert = ionToTable.insert( newElement );

        iter = returnFromInsert.first;
    }

    BetheBlochTable &betheBlochTable = iter->second;

    return betheBlochTable.GetLET(e_c_imposed);
}

//------------------------------------------------------------------------------

void fit_LQ(vector<double> dose,
            vector<double> survival,
            vector<double> survivalUncertainty,
            double &alpha, double &alphaUncertainty,
            double &beta, double &betaUncertainty,
            double &chiSquared, double &incompleteGammaQ)
{
    size_t i;

    gsl_matrix *predictorVariableMatrix = gsl_matrix_alloc( dose.size(), 2 );
    gsl_vector *survivalVector = gsl_vector_alloc( survival.size() );
    gsl_vector *weightVector = gsl_vector_alloc( survivalUncertainty.size() );
    gsl_vector *coefficientVector = gsl_vector_alloc( 2 );
    gsl_matrix *covarianceMatrix = gsl_matrix_alloc( 2, 2 );

    for( i = 0; i < dose.size(); i++ )
    {
        gsl_matrix_set(predictorVariableMatrix, i, 0, dose[i]);
        gsl_matrix_set(predictorVariableMatrix, i, 1, pow(dose[i], 2));

        gsl_vector_set(survivalVector, i, -log(survival[i]));
        gsl_vector_set(weightVector, i, 1/pow(survivalUncertainty[i]/survival[i],2));
    }

    gsl_multifit_linear_workspace *workspace = gsl_multifit_linear_alloc(dose.size(), 2);
    gsl_multifit_wlinear(predictorVariableMatrix,
                         weightVector,
                         survivalVector,
                         coefficientVector,
                         covarianceMatrix,
                         &chiSquared,
                         workspace);

    alpha = gsl_vector_get(coefficientVector, 0);
    alphaUncertainty = sqrt( gsl_matrix_get(covarianceMatrix, 0, 0) );
    beta = gsl_vector_get(coefficientVector, 1);
    betaUncertainty = sqrt( gsl_matrix_get(covarianceMatrix, 1, 1) );

    incompleteGammaQ = gsl_sf_gamma_inc_Q( (dose.size()-2.0)/2.0, chiSquared/2.0 );

    gsl_matrix_free(predictorVariableMatrix);
    gsl_vector_free(survivalVector);
    gsl_vector_free(weightVector);
    gsl_vector_free(coefficientVector);
    gsl_matrix_free(covarianceMatrix);
    gsl_multifit_linear_free(workspace);
}
    
//------------------------------------------------------------------------------
    
bool folder_exists(string foldername)
{
    struct stat st;
    int ret = stat(foldername.c_str(), &st);
    return (ret == 0) && (st.st_mode & S_IFDIR) ? true: false;
}
    
//------------------------------------------------------------------------------
    
int mkdir(const char* path)
{
    string current_level = "";
    string level;
    stringstream ss(path);
        
    // split path using slash as a separator
    while (getline(ss, level, '/'))
    {
        current_level += level; // append folder to the current level
            
        // create current level
        if (!folder_exists(current_level) && _mkdir(current_level.c_str()) != 0)
            return -1;
            
        current_level += "/";
    }
        
    return 0;
}

//------------------------------------------------------------------------------

void parse(int argc,
           char* argv[],
           string& cellType,
           string& model,
           string& trackType,
           string& parametrizationType,
           string& calculusType,
           double& precision,
           int& parallelismType,
           vector<double>& doses,
           vector<string>& parameter_name,
           double& MKM_alpha0,
           double& MKM_beta0,
           double& MKM_rNucleus,
           double& MKM_rDomain,
           double& tMKM_ac,
           double& LEM_alpha0,
           double& LEM_beta0,
           double& LEM_rNucleus,
           double& LEM_Dt,
           string& ionType,
           int& particleA,
           int& particleZ,
           string& trackMode,
           string& energyType,
           vector<double>& energies,
           int& nFraction,
           double& timeSpacing,
           double& fracDeliveryTime,
           bool& saveAlphaBeta,
           bool& saveMeans,
           bool& saveCell,
           string& projectName,
           bool& mono,
           string& spectrum_file)
{
    if (argc == 1) {
        cerr << "No input arguments." << endl;
        Usage();
        exit(1);
    }
    string tmp;
    double dump;
    bool energy_specified = false;
    bool mono_pars_specified = false;
    bool spectrum_specified = false;
    bool temporal_pars_specified = false;
    for (int i=1; i<argc; i++){
        tmp = argv[i];
        if (tmp == "-h" || tmp == "-help" || tmp == "--help") {
            Usage();
            exit(1);
        }
        if (tmp == "-cellType"){
            if (i+1 != argc){
                cellType = argv[i+1];
                i++;
            }
            else {
                cerr << "No cellType specified." << endl;
                exit(1);
            }
        }
        else if (tmp == "-model"){
            if (i+1 != argc){
                model = argv[i+1];
                i++;
            }
            else {
                cerr << "No model specified." << endl;
                exit(1);
            }
        }
        else if (tmp == "-calculusType"){
            if (i+1 != argc){
                calculusType = argv[i+1];
                i++;
            }
            else {
                cerr << "No calculusType specified." << endl;
                exit(1);
            }
        }
        else if (tmp == "-precision"){
            if (i+1 != argc){
                precision = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No precision specified." << endl;
                exit(1);
            }
            if (precision <= 0) {
                cerr << "Negative or NULL precision not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-parallelismType"){
            double check;
            if (i+1 != argc){
                check = strtod(argv[i+1], 0);
                parallelismType = (int)check;
                i++;
            }
            else {
                cerr << "No parallelismType specified." << endl;
                exit(1);
            }
            if (parallelismType < 0) {
                cerr << "Negative parallelismType not acceptable." << endl;
                exit(1);
            }
            if (modf(check, &dump) != 0.0) {
                cerr << "Invalid parallelismType selected." << endl
                     << "  Supported options are:" << endl
                     << "    0: number of threads equal to the number of core" << endl
                     << "    1: parallelism disabled" << endl
                     << "    n: n threads" << endl;
                exit(1);
            }
        }
        else if (tmp == "-doses"){
            double next_dose;
            if (i+1 != argc){
                tmp = argv[i+1];
                if (tmp.find("-")!=string::npos) {
                    cerr << "No acceptable doses specified." << endl;
                    exit(1);
                }
            }
            else {
                cerr << "No doses specified." << endl;
                exit(1);
            }
            while (tmp.find("-")==string::npos) {
                next_dose = strtod(argv[i+1], 0);
                if (next_dose <= 0) {
                    cerr << "Invalid negative or NULL dose." << endl;
                    exit(1);
                }
                else {
                    doses.push_back(next_dose);
                    i++;
                    if (i+1 != argc)
                        tmp = argv[i+1];
                    else
                        break;
                }
            }
        }
        else if (tmp == "-MKM_alpha0"){
            if (i+1 != argc){
                MKM_alpha0 = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No MKM_alpha0 specified." << endl;
                exit(1);
            }
            if (MKM_alpha0 < 0) {
                cerr << "Negative MKM_alpha0 not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-MKM_beta0"){
            if (i+1 != argc){
                MKM_beta0 = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No MKM_beta0 specified." << endl;
                exit(1);
            }
            if (MKM_beta0 < 0) {
                cerr << "Negative MKM_beta0 not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-MKM_rNucleus"){
            if (i+1 != argc){
                MKM_rNucleus = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No MKM_rNucleus specified." << endl;
                exit(1);
            }
            if (MKM_rNucleus <= 0) {
                cerr << "Negative or NULL MKM_rNucleus not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-MKM_rDomain"){
            if (i+1 != argc){
                MKM_rDomain = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No MKM_rDomain specified." << endl;
                exit(1);
            }
            if (MKM_rDomain <= 0) {
                cerr << "Negative or NULL MKM_rDomain not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-tMKM_timeConst"){
            if (i+1 != argc){
                tMKM_ac = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No tMKM_timeConst specified." << endl;
                exit(1);
            }
            if (MKM_rDomain <= 0) {
                cerr << "Negative or NULL tMKM_timeConst not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-LEM_alpha0"){
            if (i+1 != argc){
                LEM_alpha0 = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No LEM_alpha0 specified." << endl;
                exit(1);
            }
            if (LEM_alpha0 < 0) {
                cerr << "Negative LEM_alpha0 not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-LEM_beta0"){
            if (i+1 != argc){
                LEM_beta0 = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No LEM_beta0 specified." << endl;
                exit(1);
            }
            if (LEM_beta0 < 0) {
                cerr << "Negative LEM_beta0 not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-LEM_rNucleus"){
            if (i+1 != argc){
                LEM_rNucleus = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No LEM_rNucleus specified." << endl;
                exit(1);
            }
            if (LEM_rNucleus <= 0) {
                cerr << "Negative or NULL LEM_rNucleus not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-LEM_Dt"){
            if (i+1 != argc){
                LEM_Dt = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No LEM_Dt specified." << endl;
                exit(1);
            }
            if (LEM_Dt < 0) {
                cerr << "Negative LEM_Dt not acceptable." << endl;
                exit(1);
            }
        }
        else if (tmp == "-ion"){
            if (spectrum_specified) {
                cerr << "Invalid \"-ion\" option: spectrum already specified." << endl;
                exit(1);
            }
            if (i+1 != argc){
                ionType = argv[i+1];
                i++;
            }
            else {
                cerr << "Ion not specified." << endl;
                exit(1);
            }
            if( ionType == "H" )
            {
                particleA = 1;
                particleZ = 1;
            }
            else if( ionType == "He" )
            {
                particleA = 4;
                particleZ = 2;
            }
            else if( ionType == "Li" )
            {
                particleA = 7;
                particleZ = 3;
            }
            else if( ionType == "Be" )
            {
                particleA = 9;
                particleZ = 4;
            }
            else if( ionType == "B" )
            {
                particleA = 10;
                particleZ = 5;
            }
            else if( ionType == "C" )
            {
                particleA = 12;
                particleZ = 6;
            }
            else if(ionType == "N"  )
            {
                particleA = 14;
                particleZ = 7;
            }
            else if( ionType == "O" )
            {
                particleA = 16;
                particleZ = 8;
            }
            else if(ionType == "F"  )
            {
                particleA = 19;
                particleZ = 9;
            }
            else if(ionType == "Ne"  )
            {
                particleA = 20;
                particleZ = 10;
            }
            else {
                cerr << "Ion: " << ionType << " not supported." << endl;
                exit(1);
            }
            mono_pars_specified = true;
        }
        else if (tmp == "-trackMode"){
            if (i+1 != argc){
                trackMode = argv[i+1];
                i++;
            }
            else {
                cerr << "No trackMode specified." << endl;
                exit(1);
            }
            if (trackMode != "histogram" && trackMode != "random"){
                cerr << "Track Mode: " << trackMode << " not supported. Choose between \"histogram\" and \"random\"." << endl;
                exit(1);
            }
        }
        else if (tmp == "-energies"){
            if (energy_specified) {
                cerr << "Invalid input \"-energies\": LETs already specified (\"-lets\" option)." << endl;
                exit(1);
            }
            if (spectrum_specified) {
                cerr << "Invalid \"-energies\" option: spectrum already specified." << endl;
                exit(1);
            }
            energyType = "energy";
            if (i+1 != argc){
                tmp = argv[i+1];
                if (tmp.find("-")!=string::npos) {
                    cerr << "No acceptable energies specified." << endl;
                    exit(1);
                }
            }
            else {
                cerr << "No energies specified." << endl;
                exit(1);
            }
            while (tmp.find("-")==string::npos) {
                energies.push_back(strtod(argv[i+1], 0));
                i++;
                if (i+1 != argc)
                    tmp = argv[i+1];
                else
                    break;
            }
            energy_specified = true;
            mono_pars_specified = true;
            
        }
        else if (tmp == "-lets"){
            if (energy_specified) {
                cerr << "Invalid \"-lets\" option: Energies already specified (\"-energies\" option)." << endl;
                exit(1);
            }
            if (spectrum_specified) {
                cerr << "Invalid \"-lets\" option: spectrum already specified." << endl;
                exit(1);
            }
            energyType = "LET";
            if (i+1 != argc){
                tmp = argv[i+1];
                if (tmp.find("-")!=string::npos) {
                    cerr << "No acceptable lets specified." << endl;
                    exit(1);
                }
            }
            else {
                cerr << "No lets specified." << endl;
                exit(1);
            }
            while (tmp.find("-")==string::npos) {
                energies.push_back(strtod(argv[i+1], 0));
                i++;
                if (i+1 != argc)
                    tmp = argv[i+1];
                else
                    break;
            }
            energy_specified = true;
            mono_pars_specified = true;
        }
        else if (tmp == "-nFraction"){
            double check;
            if (i+1 != argc){
                check = strtod(argv[i+1], 0);
                nFraction = check;
                i++;
            }
            else {
                cerr << "No nFraction specified." << endl;
                exit(1);
            }
            if (nFraction <= 0) {
                cerr << "Number of Fractions: Input not acceptable. Negative or NULL value." << endl;
                exit(1);
            }
            if (modf(check, &dump) != 0.0) {
                cerr << "nFraction has to be an integer." << endl;
                exit(1);
            }
            temporal_pars_specified = true;
        }
        else if (tmp == "-timeSpacing"){
            if (i+1 != argc){
                timeSpacing = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No timeSpacing specified." << endl;
                exit(1);
            }
            if (timeSpacing < 0) {
                cerr << "Time Spacing: Input not acceptable. Negative value." << endl;
                exit(1);
            }
            temporal_pars_specified = true;
        }
        else if (tmp == "-fracDeliveryTime"){
            if (i+1 != argc){
                fracDeliveryTime = strtod(argv[i+1], 0);
                i++;
            }
            else {
                cerr << "No fracDeliveryTime specified." << endl;
                exit(1);
            }
            if (fracDeliveryTime < 0) {
                cerr << "Fraction Delivery Time: Input not acceptable. Negative value." << endl;
                exit(1);
            }
            temporal_pars_specified = true;
        }
        else if (tmp == "-output"){
            vector<string> outOptions;
            if (i+1 != argc){
                tmp = argv[i+1];
                if (tmp.find("-")!=string::npos) {
                    cerr << "No acceptable output specified." << endl;
                    exit(1);
                }
            }
            else {
                cerr << "No output specified." << endl;
                exit(1);
            }
            while (tmp.find("-")==string::npos) {
                outOptions.push_back(argv[i+1]);
                i++;
                if (i+1 != argc)
                    tmp = argv[i+1];
                else
                    break;
            }
            for (unsigned int a=0; a<outOptions.size(); a++) {
                if (outOptions[a] == "LQ_pars")
                    saveAlphaBeta = true;
                else if (outOptions[a] == "meanValues")
                    saveMeans = true;
                else if (outOptions[a] == "cellValues")
                    saveCell = true;
                else if (outOptions[a] == "all") {
                    saveAlphaBeta = true;
                    saveMeans = true;
                    saveCell = true;
                }
                else {
                    cerr << "Output options " << outOptions[a] << " not supported." << endl
                         << "Choose between: \"LQ_pars\", \"meanValues\" and \"cellValues\"." << endl;
                    exit(1);
                }
            }
        }
        else if (tmp == "-projectName"){
            if (i+1 != argc){
                projectName = argv[i+1];
                i++;
            }
            else {
                cerr << "No projectName specified." << endl;
                exit(1);
            }
        }
        else if (tmp == "-spectrum_file"){
            if (mono_pars_specified) {
                cerr << "Invalid \"-spectrum_file\" option: Monoenergetic parameters already specified (\"-ionType\", \"-energies\" or \"-lets\")." << endl;
                exit(1);
            }
            if (i+1 != argc){
                spectrum_file = argv[i+1];
                i++;
            }
            else {
                cerr << "No spectrum specified." << endl;
                exit(1);
            }
            spectrum_specified = true;
            mono = false;
        }
        else {
            cerr << endl << "Invalid argument: " << tmp << endl;
            Usage();
            exit(1);
        }
    }
    if (!energy_specified && !spectrum_specified) {
        cerr << "No LETs/Energies specified." << endl;
        exit(1);
    }
    
    if (mono_pars_specified)
        mono = true;
    
    if (doses.size() == 0)
        for (int a=1; a<=6; a++)
            doses.push_back(a);
    
    if (fracDeliveryTime > timeSpacing && timeSpacing!=0) {
        cerr << "Time Spacing smaller than Fraction Delivery Time." << endl;
        exit(1);
    }
    
    if( model == "LEMI" )
    {
        trackType = "Scholz2000";
        parametrizationType = "LQ";
        
        if ( (calculusType != "rapidLEM_Russo2011") && (calculusType != "rapidLEM_Scholz2006") && (calculusType != "MonteCarlo") ) {
            cerr << "The indicated calculusType, " << calculusType << " is not compatible with the specified model, " << model << "." << endl
                 << "Choose between \"rapidLEM_Russo2011\", \"rapidLEM_Scholz2006\" and \"MonteCarlo\"." << endl;
            exit(1);
        }
        parameter_name.push_back("alpha_0");
        parameter_name.push_back("beta_0");
        parameter_name.push_back("r_nucleus");
        parameter_name.push_back("D_t");
    }
    else if ( model == "LEMII" )
    {
        trackType = "Elsasser2007";
        parametrizationType = "LQ2";
        
        if ( (calculusType != "rapidLEM_Russo2011") && (calculusType != "rapidLEM_Scholz2006") && (calculusType != "MonteCarlo")) {
            cerr << "The indicated calculusType, " << calculusType << ", is not compatible with the specified model, " << model << "." << endl
                 << "Choose between \"rapidLEM_Russo2011\", \"rapidLEM_Scholz2006\" and \"MonteCarlo\"." << endl;
            exit(1);
        }
        parameter_name.push_back("alpha_0");
        parameter_name.push_back("beta_0");
        parameter_name.push_back("r_nucleus");
        parameter_name.push_back("D_t");
    }
    else if ( model == "LEMIII" )
    {
        trackType = "Elsasser2008";
        parametrizationType = "LQ3";
        
        if ( (calculusType != "rapidLEM_Russo2011") && (calculusType != "rapidLEM_Scholz2006") && (calculusType != "MonteCarlo")) {
            cerr << "The indicated calculusType, " << calculusType << ", is not compatible with the specified model, " << model << "." << endl
                 << "Choose between \"rapidLEM_Russo2011\", \"rapidLEM_Scholz2006\" and \"MonteCarlo\"." << endl;
            exit(1);
        }
        parameter_name.push_back("alpha_0");
        parameter_name.push_back("beta_0");
        parameter_name.push_back("r_nucleus");
        parameter_name.push_back("D_t");
    }
    else if ( model == "MKM" )
    {
        trackType = "KieferChatterjee";
        parametrizationType = "LQ_noDt";
        
        if ( (calculusType != "rapidMKM_Kase2008")
               && (calculusType != "rapidMKM_Attili2013")
               && (calculusType != "rapidMKM_Kase2008_corrected_beta")
               && (calculusType != "rapidMKM_Attili2013_corrected_beta")
               && (calculusType != "MonteCarlo")) {
            cerr << "The indicated calculusType, " << calculusType << " is not compatible with the specified model, " << model << "." << endl
                 << "Choose between \"rapidMKM_Kase2008\",  \"rapidMKM_Attili2013\",  \"rapidMKM_Kase2008_corrected_beta\",  \"rapidMKM_Attili2013_corrected_beta\" and \"MonteCarlo\"." << endl;
            exit(1);
        }
        parameter_name.push_back("alpha_0");
        parameter_name.push_back("beta_0");
        parameter_name.push_back("r_nucleus");
        parameter_name.push_back("r_domain");
    }
    else if ( model == "tMKM_Manganaro2017" )
    {
        trackType = "KieferChatterjee";
        parametrizationType = "LQ_noDt_T";
        
        if ( calculusType != "MonteCarlo") {
            cerr << "The indicated calculusType, " << calculusType << " is not compatible with the specified model, " << model << "." << endl
                 << "Only \"MonteCarlo\" supported." << endl;
            exit(1);
        }
        parameter_name.push_back("alpha_0");
        parameter_name.push_back("beta_0");
        parameter_name.push_back("r_nucleus");
        parameter_name.push_back("r_domain");
        parameter_name.push_back("time_constant");
    }
    else
    {
        cerr << "The indicated model is not provided." << endl
             << "Choose between LEMI, LEMII, LEMIII, MKM and tMKM_Manganaro2017." << endl;
        exit(1);
    }
    
    if (temporal_pars_specified && model!="tMKM_Manganaro2017") {
        cerr << "The specified model (" << model << ") doesn't support temporal parameters." << endl
             << "Choose \"tMKM_Manganaro2017\" model or remove temporal specifications (\"-nFraction\", \"-timeSpacing\" or \"-fracDeliveryTime\")." << endl;
        exit(1);
    }
    
    if ( (calculusType == "rapidMKM_Kase2008" || calculusType == "rapidMKM_Attili2013" || calculusType == "rapidMKM_Kase2008_corrected_beta" || calculusType == "rapidMKM_Attili2013_corrected_beta"
            || calculusType == "rapidLEM_Scholz2006" || calculusType == "rapidLEM_Russo2011")
            && saveCell ) {
        cerr << "Output not supported from the indicated calculusType: " << calculusType << "." << endl
             << "Choose \"LQ_pars\" or change calculusType." << endl;
        exit(1);
    }
    if (saveAlphaBeta && doses.size()<=1 && calculusType=="MonteCarlo"){
        cerr << "At least 2 values of dose needed to evaluate the LQ parameters." << endl;
        exit(1);
    }
    
    if (!saveAlphaBeta && !saveMeans && !saveCell){
        cerr << "No output specified." << endl;
        Usage();
        exit(1);
    }
    
    cout << endl << "-----------------------------------------------------" << endl << endl
         << "SELECTED PARAMETERS:" << endl
         << "   Cell Type: " << cellType << endl
         << endl
         << "Model and Calculus setting:" << endl
         << "   Model: " << model << endl
         << "   Model parameters: " << endl
         << "   (MKM and tMKM)" << ((model == "MKM" || model == "tMKM_Manganaro2017") ? " " : " - Not used") << endl
         << "      Alpha0: " << MKM_alpha0 << endl
         << "      Beta0:  " << MKM_beta0 << endl
         << "      Nucleus Radius: " << MKM_rNucleus << endl
         << "      Domain Radius:  " << MKM_rDomain << endl
         << "   (tMKM)" << (model == "tMKM_Manganaro2017" ? " " : " - Not used") << endl
         << "      Repair Time Constant: " << tMKM_ac << endl
         << "   (LEMI, LEMII and LEMIII)" << ((model == "MKM" || model == "tMKM_Manganaro2017") ? " - Not used" : " ") << endl
         << "      alpha0: " << LEM_alpha0 << endl
         << "      beta0:  " << LEM_beta0 << endl
         << "      Nucleus Radius: " << LEM_rNucleus << endl
         << "      Threshold dose: " << LEM_Dt << endl
         << "   Calculus Type: " << calculusType << endl
         << "   Precision: " << precision << endl
         << "   Parallelism Type: " << parallelismType << endl
         << "   Doses: ";
    for (unsigned int i=0; i<doses.size(); i++)
        cout << doses[i] << " ";
    cout << endl
         << "   Parametrization Type: " << parametrizationType << endl
         << endl
         << "Radiation settings:" << endl
         << "   Radiation type: ";
    if (mono)
    {
        cout << "Monoenergetic" << endl
             << "   Ion: " << ionType << endl
             << "   Particle A: " << particleA << endl
             << "   Particle Z: " << particleZ << endl
             << "   Energy Type: " << energyType << " (" << energies.size() << " values specified)." << endl
             << (energyType == "LET" ? "   LETs:" : "   Energies:") << endl;
        for (unsigned int i=0; i<energies.size(); i++)
            cout << "     " << energies[i] << " " << (energyType == "LET" ? "keV/um" : "MeV") << endl;
        cout << endl;
    }
    else
        cout << "Spectrum" << endl
             << "   File: " << spectrum_file << endl;
    cout << "   Track Type: " << trackType << endl
         << "   Track Mode: " << trackMode << endl
         << endl
         << "Fractionation Scheme:" << endl
         << "   Number of Fractions: " << nFraction << endl
         << "   Time Spacing: " << timeSpacing << endl
         << "   Fraction Delivery Time: " << fracDeliveryTime << endl
         << endl
         << "Output:" << endl;
    if (saveAlphaBeta) cout << "   - LQ parameters." << endl;
    if (saveMeans) cout << "   - Mean dose and survival values." << endl;
    if (saveCell) cout << "   - Cell by cell dose and survival." << endl;
    cout << "   Project Name: " << projectName << endl
         << endl << "-----------------------------------------------------" << endl << endl;
    
    return;
}

//------------------------------------------------------------------------------

void Usage()
{
    cerr << endl
         << "USAGE: " << endl
         << endl
         << "Model and Calculus settings:" << endl
         << "     survival -model [string]" << endl
         << "          Supported: \"LEMI\" - \"LEMII\" - \"LEMIII\" - \"MKM\" - \"tMKM_Manganaro2017\"" << endl
         << "     survival -calculusType [string]" << endl
         << "          Supported:" << endl
         << "             \"rapidLEM_Scholz2006\" (compatible with LEMI, LEMII and LEMIII)" << endl
         << "             \"rapidLEM_Russo2011\"  (compatible with LEMI, LEMII and LEMIII)" << endl
         << "             \"rapidMKM_Kase2008\"    (compatible with MKM)" << endl
         << "             \"rapidMKM_Attili2013\"    (compatible with MKM)" << endl
         << "             \"rapidMKM_Kase2008_corrected_beta\"    (compatible with MKM)" << endl
         << "             \"rapidMKM_Attili2013_corrected_beta\"    (compatible with MKM)" << endl
         << "             \"MonteCarlo\"  (compatible with LEMI, LEMII, LEMIII, MKM and tMKM_Manganaro2017)" << endl
         << "     survival -precision [double]" << endl
         << "     survival -parallelismType [int]" << endl
         << "          Supported:" << endl
         << "             0: number of threads equal to the number of core" << endl
         << "             1: parallelism disabled" << endl
         << "             n: n threads" << endl
         << "     survival -doses [double] ... [double]" << endl
         << endl
         << "Model Parameters:" << endl
         << "   (MKM and tMKM)" << endl
         << "     survival -MKM_alpha0 [double]" << endl
         << "     survival -MKM_beta0  [double]" << endl
         << "     survival -MKM_rNucleus [double]" << endl
         << "     survival -MKM_rDomain  [double]" << endl
         << "   (tMKM)" << endl
         << "     survival -tMKM_timeConst [double]" << endl
         << "   (LEMI, LEMII and LEMIII)" << endl
         << "     survival -LEM_alpha0 [double]" << endl
         << "     survival -LEM_beta0  [double]" << endl
         << "     survival -LEM_rNucleus [double]" << endl
         << "     survival -LEM_Dt       [double]" << endl
         << endl
         << "Cell and Radiation settings:" << endl
         << "     survival -cellType [string]" << endl
         << "     survival -trackMode  string]" << endl
         << "          Supported: \"histogram\" - \"random\" " << endl
         << "   (Monoenergetic radiation)" << endl
         << "     survival -ion  [string]" << endl
         << "          Supported: from Z=1 to Z=10 (H, He, ..., Ne)" << endl
         << "     survival -energies [double] ... [double]" << endl
         << "          => Kinetic energies expressed in MeV. NOT compatible with \"-lets\" option." << endl
         << "     survival -lets [double] ... [double] (values in keV/um). NOT compatible with \"-energies\" option." << endl
         << "   (Spectrum)" << endl
         << "     survival -spectrum_file [string] (The name of the file containing the spectrum definition - With absolute or relative path)"
         << endl
         << "Fractionation Scheme:" << endl
         << "     survival -nFraction [int]" << endl
         << "     survival -timeSpacing [double]" << endl
         << "     survival -fracDeliveryTime [double]" << endl
         << endl
         << "Output:" << endl
         << "     survival -projectName [string]" << endl
         << "     survival -output [string] ... [string]" << endl
         << "          Supported: \"LQ_pars\", \"meanValues\" and \"cellValues\"" << endl
         << "                     Choose \"all\" to enable all output types." << endl
         << endl
         << "Help:" << endl
         << "     survival --help   ==>   Display this hint." << endl;
}
    
}
