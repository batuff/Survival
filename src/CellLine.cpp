#include "CellLine.h"
#include "usefulFunctions.h"

#include <fstream>
using std::ifstream;
using std::ofstream;

#define _USE_MATH_DEFINES
#include <cmath>
using std::pow;
using std::abs;

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
using std::flush;
using std::ios;
using std::clog;

#include <cstdlib>
using std::exit;

using std::string;
using std::vector;

using namespace Survival;

bool CellLine::needEtaGenerated = false;

CellLine::CellLine(const string cell_type,
                   const double nucleus_radius,
                   const double domain_radius)
: cellType( cell_type ),
  nucleusRadius( nucleus_radius ),
  domainRadius( domain_radius ),
  alpha_X( 0.0 ),
  beta_X( 0.0 ),
  isLQ_noDtLoaded( false ),
  ac( 0.0 ),
  isLQ_noDt_TLoaded( false ),
  alpha_X1( 0.0 ),
  beta_X1( 0.0),
  D_t( 0.0 ),
  s( 0.0 ),
  logS_t( 0.0 ),
  isLQloaded( false ),
  alpha_X2( 0.0 ),
  beta_X2( 0.0 ),
  D_t2( 0.0 ),
  s2( 0.0 ),
  logS_t2( 0.0 ),
  genomeLength( 0.0 ),
  alpha_SSB( 0.0 ),
  alpha_DSB( 0.0 ),
  base_Pairs( 0 ),
  isLQ2loaded( false ),
  alpha_X3( 0.0 ),
  beta_X3( 0.0 ),
  D_t3( 0.0 ),
  s3( 0.0 ),
  logS_t3( 0.0 ),
  isLQ3loaded( false )
{
    selectedParametrization = &CellLine::noParametrization;
    selectedParametrizationT = &CellLine::noParametrization;
}

//------------------------------------------------------------------------------

void CellLine::addParametrization_LQ_noDt(const double alphaX,
                                          const double betaX)
{
    alpha_X = alphaX;  //Gy^(-1)
    beta_X = betaX;  //Gy^(-2)
    
    isLQ_noDtLoaded = true;
    isLQ_noDt_TLoaded = true; //serve solo per creare il nucleo nel global scope...
}

//------------------------------------------------------------------------------

void CellLine::addParametrization_LQ_noDt_T(const double alphaX,
                                            const double betaX,
                                            const double ac_)
{
    alpha_X = alphaX;  //Gy^(-1)
    beta_X = betaX;  //Gy^(-2)
    ac = ac_;
    
    isLQ_noDtLoaded = true; //serve solo per creare il nucleo nel global scope...
    isLQ_noDt_TLoaded = true;
}

//------------------------------------------------------------------------------

void CellLine::addParametrization_LQ(const double alphaX1,
                                     const double betaX1,
                                     const double Dt)
{
    alpha_X1 = alphaX1;  //Gy^(-1)
    beta_X1 = betaX1;  //Gy^(-2)
    D_t = Dt;  //Gy
    
    logS_t = - alpha_X1 * D_t - beta_X1 * pow(D_t, 2);
    s = alpha_X1 + 2 * beta_X1 * D_t;  //Gy^(-1)
    
    isLQloaded = true;
    isLQ_noDtLoaded = true; //serve solo per creare il nucleo nel global scope...
    isLQ_noDt_TLoaded = true; //serve solo per creare il nucleo nel global scope...
}

//------------------------------------------------------------------------------

void CellLine::addParametrization_LQ2(const double alphaX2,
                                      const double betaX2,
                                      const double Dt2,
                                      const double genome_Length,
                                      const double alphaSSB,
                                      const double alphaDSB,
                                      long int basePairs)
{
    alpha_X2 = alphaX2;  //Gy^(-1)
    beta_X2 = betaX2;  //Gy^(-2)
    D_t2 = Dt2;  //Gy
    genomeLength = genome_Length;
    alpha_SSB = alphaSSB;
    alpha_DSB = alphaDSB;
    base_Pairs = basePairs;
    
    logS_t2 = exp(- alpha_X2 * D_t2 - beta_X2 * pow(D_t2, 2));
    s2 = alpha_X2 + 2 * beta_X2 * D_t2;  //Gy^(-1)
    
    isLQ2loaded = true;
    isLQ_noDtLoaded = true; //serve solo per creare il nucleo nel global scope...
    isLQ_noDt_TLoaded = true; //serve solo per creare il nucleo nel global scope...
}

//------------------------------------------------------------------------------

void CellLine::addParametrization_LQ3(const double alphaX3,
                                      const double betaX3,
                                      const double Dt3,
                                      const double genome_Length,
                                      const double alphaSSB,
                                      const double alphaDSB,
                                      long int basePairs)
{
    alpha_X3 = alphaX3;  //Gy^(-1)
    beta_X3 = betaX3;  //Gy^(-2)
    D_t3 = Dt3;  //Gy
    genomeLength = genome_Length;
    alpha_SSB = alphaSSB;
    alpha_DSB = alphaDSB;
    base_Pairs = basePairs;
    
    logS_t3 = exp(- alpha_X3 * D_t3 - beta_X3 * pow(D_t3, 2));
    s3 = alpha_X3 + 2 * beta_X3 * D_t3;  //Gy^(-1)
    
    isLQ3loaded = true;
    isLQ_noDtLoaded = true; //serve solo per creare il nucleo nel global scope...
    isLQ_noDt_TLoaded = true; //serve solo per creare il nucleo nel global scope...
}

//------------------------------------------------------------------------------

double CellLine::analyticDamageEnhancement(const double dose) const
{
    double N, i;
    double factN;
    double a1 = alpha_SSB * (base_Pairs/genomeLength);
    double a2 = alpha_DSB * (base_Pairs/genomeLength);
    double analyticEta;
    
    analyticEta = 1;
    for (N=1;N<=base_Pairs;N++)
    {
        factN = 1;
        for (i=2;i<=N;i++)
            factN = i * factN;
        analyticEta = analyticEta + (1/factN) * (exp(-dose * ( a1 + a2 ))) * (pow((a1 * dose),double(N)))*(1-pow(double(2),double(1-N)))/(a2 * dose);
    }
    return analyticEta;
}

//------------------------------------------------------------------------------

int CellLine::DNA[10000000];
bool CellLine::SSB1[10000000];
bool CellLine::SSB2[10000000];
bool CellLine::DSB[10000000];

//------------------------------------------------------------------------------

double CellLine::damageEnhancement(const double dose) const
{
    if( dose >= 1e2 )
    {
        int bpSelection = 10000000; // must be less than 10000000
        double MeanNSSB;
        double MeanNDSB;
        double which;
        
        double eta;
        
        long int position, position_min, position_max;
        int i, j;
        int Found, Found2;
        long int FoundPositions[2*base_Pairs];
        long int NSSB, NDSB, N2SB, nNSSB, nNDSB;
        
        NSSB=0;
        NDSB=0;
        N2SB=0;
        
        for (i=0;i<bpSelection;i++) // generation of DSBs
        {
            DNA[i]=0;
            SSB1[i]=false;
            SSB2[i]=false;
            DSB[i]=false;
        }
        MeanNSSB = alpha_SSB * dose * (double(bpSelection)/double(genomeLength));
        MeanNDSB = alpha_DSB * dose * (double(bpSelection)/double(genomeLength));
        nNSSB = (long int) MeanNSSB;
        nNDSB = (long int) MeanNDSB;
        
        NSSB += nNSSB;
        NDSB += nNDSB;
        
        for (i=0;i<nNDSB;i++)
        {
            position = (long int) floor((double(rand())/double(RAND_MAX)) * double(bpSelection));
            DSB[position]=true;
        }
        for (i=0;i<nNSSB;i++) // generation of SSBs on strands
        {
            position = (long int) floor((double(rand())/double(RAND_MAX)) * double(bpSelection));
            which = double(rand())/double(RAND_MAX);
            
            if (which<0.5) // case of break on first strand
            {
                SSB1[position]=true;
                
                if (DNA[position]==0)
                {
                    Found=0;
                    Found2=0;
                    FoundPositions[Found]=position;
                    position_min=position-base_Pairs/2;
                    position_max=position+base_Pairs/2;
                    
                    if (position-base_Pairs/2<0) {position_min=0;}
                    if (position_max>bpSelection) {position_max=bpSelection-1;}
                    
                    for (j=position_min;j<=position_max;j++)
                    {
                        if (SSB2[j]==true && DNA[j]==0)
                        {
                            Found++;
                            FoundPositions[Found]=j;
                        }
                        if (DSB[j]==true) {Found2++;}
                    }
                    
                    if (Found>0 && Found2==0)
                    {
                        position=(FoundPositions[1]+position)/2;
                        position_min=position-base_Pairs/2;
                        position_max=position+base_Pairs/2;
                        if (position-base_Pairs<0) {position_min=0;}
                        if (position_max>bpSelection) {position_max=bpSelection-1;}
                        for (j=position_min;j<=position_max;j++)
                            DNA[j]=Found;
                        N2SB++;
                    }
                }
            }
            else  // case of break on second strand
            {
                SSB2[position]=true;
                
                if (DNA[position]==0)
                {
                    Found=0;
                    Found2=0;
                    FoundPositions[Found]=position;
                    position_min=position-base_Pairs;
                    position_max=position+base_Pairs;
                    
                    if (position-base_Pairs<0) {position_min=0;}
                    if (position_max>bpSelection) {position_max=bpSelection-1;}
                    
                    for (j=position_min;j<=position_max;j++)
                    {
                        if (SSB1[j]==true && DNA[j]==0)
                        {
                            Found++;
                            FoundPositions[Found]=j;
                        }
                        if (DSB[j]==true) {Found2++;}
                    }
                    
                    if (Found>0 && Found2==0)
                    {
                        position=position+abs(static_cast<double>(FoundPositions[1]-position));
                        position_min=position-base_Pairs/2;
                        position_max=position+base_Pairs/2;
                        if (position-base_Pairs<0) {position_min=0;}
                        if (position_max>bpSelection) {position_max=bpSelection-1;}
                        for (j=position_min;j<=position_max;j++)
                            DNA[j]=Found;
                        N2SB++;
                    }
                }
            }
        }
        eta = (double(NDSB) + double(N2SB))/double(NDSB);
        
        return eta;
    }
    else
        return 1;
}

//------------------------------------------------------------------------------

double CellLine::getLogSurvival_X(const double dose) const
{
    return (*this.*selectedParametrization)(dose);
}

//------------------------------------------------------------------------------

double CellLine::getLogSurvival_X(const std::vector<double>doses,
                                  const std::vector<double>times) const
{
    if (selectedParametrizationT == &CellLine::parametrization_LQ_noDt_T)
        return (*this.*selectedParametrizationT)(doses, times);
    else {
        cerr << "Selected parametrization does not match \"getLogSurvival_X\" operations." << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

void CellLine::getParameters(double &returnAlpha_X,
                             double &returnBeta_X,
                             double &returnD_t) const
{
    if( selectedParametrization == &CellLine::parametrization_LQ_noDt ||
       selectedParametrizationT == &CellLine::parametrization_LQ_noDt_T )
    {
        returnAlpha_X = alpha_X;
        returnBeta_X = beta_X;
        returnD_t = -1.0;
    }
    else if( selectedParametrization == &CellLine::parametrization_LQ )
    {
        returnAlpha_X = alpha_X1;
        returnBeta_X = beta_X1;
        returnD_t = D_t;
    }
    else if( selectedParametrization == &CellLine::parametrization_LQ3 )
    {
        returnAlpha_X = alpha_X3;
        returnBeta_X = beta_X3;
        returnD_t = D_t3;
    }
    else
    {
        returnAlpha_X = alpha_X2;
        returnBeta_X = beta_X2;
        returnD_t = D_t2;
    }
}

//------------------------------------------------------------------------------

void CellLine::getParameters_LQ_noDt(double &returnAlpha_X,
                                     double &returnBeta_X) const
{
    if( isLQ_noDtLoaded )
    {
        returnAlpha_X = alpha_X;
        returnBeta_X = beta_X;
    }
    else
    {
        cerr << "CellLine: the LQ parametrization has not been loaded." << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

void CellLine::getParameters_LQ_noDt_T(double &returnAlpha_X,
                                       double &returnBeta_X,
                                       double &ac_) const
{
    if( isLQ_noDt_TLoaded )
    {
        returnAlpha_X = alpha_X;
        returnBeta_X = beta_X;
        ac_ = ac;
    }
    else
    {
        cerr << "CellLine: the LQ parametrization has not been loaded." << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

void CellLine::getParameters_LQ2(double &returnAlpha_X2,
                                 double &returnBeta_X2,
                                 double &returnD_t2,
                                 double &returnGenome_Length,
                                 double &returnAlpha_SSB,
                                 double &returnAlpha_DSB,
                                 long int &returnBase_Pairs) const
{
    if( isLQ2loaded )
    {
        returnAlpha_X2 = alpha_X2;
        returnBeta_X2 = beta_X2;
        returnD_t2 = D_t2;
        returnGenome_Length = genomeLength;
        returnAlpha_SSB = alpha_SSB;
        returnAlpha_DSB = alpha_DSB;
        returnBase_Pairs = base_Pairs;
    }
    else
    {
        cerr << "CellLine: the LQ2 parametrization has not been loaded." << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

void CellLine::getParameters_LQ3(double &returnAlpha_X3,
                                 double &returnBeta_X3,
                                 double &returnD_t3,
                                 double &returnGenome_Length,
                                 double &returnAlpha_SSB,
                                 double &returnAlpha_DSB,
                                 long int &returnBase_Pairs) const
{
    if( isLQ3loaded )
    {
        returnAlpha_X3 = alpha_X3;
        returnBeta_X3 = beta_X3;
        returnD_t3 = D_t3;
        returnGenome_Length = genomeLength;
        returnAlpha_SSB = alpha_SSB;
        returnAlpha_DSB = alpha_DSB;
        returnBase_Pairs = base_Pairs;
    }
    else
    {
        cerr << "CellLine: the LQ3 parametrization has not been loaded." << endl;
        exit(1);
    }
}

//------------------------------------------------------------------------------

double CellLine::interpolatedDamageEnhancement(const double dose) const
{
    double interpolationSlope;
    int points = 200;
    
    if( dose > 1e2 && dose <= 5e6)
    {
        int k;
        k = int(ceil( ( (log10(dose) - log10(1e2))/(log10(5e6) - log10(1e2)) ) * points ));
        
        interpolationSlope = ( log10(etaPre[k]) - log10(etaPre[k-1]) ) / ( log10(doseForEta[k]) - log10(doseForEta[k-1]) );
        
        return pow( double(10) , (interpolationSlope * (log10(dose) - log10(doseForEta[k-1])) + log10(etaPre[k-1])) );
    }
    else
        return 1;
}

//------------------------------------------------------------------------------

double CellLine::noParametrization(const double) const
{
    cerr << "CellLine: There is no parametrization loaded or selected." << endl;
    exit(1);
    
    return 0.0;
}

//------------------------------------------------------------------------------

double CellLine::noParametrization(const vector<double>,
                                   const vector<double>) const
{
    cerr << "CellLine: There is no parametrization loaded or selected." << endl;
    exit(1);
    
    return 0.0;
}

//------------------------------------------------------------------------------

double CellLine::parametrization_LQ_noDt(const double dose) const
{
    return - alpha_X * dose - beta_X * dose*dose;
}

//------------------------------------------------------------------------------

double CellLine::parametrization_LQ_noDt_T(const vector<double>doses,
                                           const vector<double>times) const
{
    double logSurv=0.0, sumD=0.0;
    if (times.size()!=0) {
        for (unsigned long int i=0; i<times.size()-1; ++i) {
            for (unsigned long int j=i+1; j<times.size(); ++j)
                logSurv += 2*beta_X * doses[i]*doses[j] * ( 1-exp(-ac*(times[j]-times[i])) );
            sumD+=doses[i];
        }
        sumD += doses[times.size()-1];
    }
    logSurv += - alpha_X * sumD - beta_X * sumD*sumD;
    
    return logSurv;
}

//------------------------------------------------------------------------------

double CellLine::parametrization_LQ(const double dose) const
{
    if (dose <= D_t)
        return - alpha_X1 * dose - beta_X1 * dose*dose;
    else
        return logS_t - s * (dose - D_t);
}

//------------------------------------------------------------------------------

double CellLine::parametrization_LQ2(const double dose) const
{
    double eta = (*this.*selectedDamageEnhancement)(dose);
    
    if (dose <= D_t2)
        return - alpha_X2 * dose - beta_X2 * dose*dose;
    else
        return logS_t2 - s2 * (eta * dose - D_t2);
}

//------------------------------------------------------------------------------

double CellLine::parametrization_LQ3(const double dose) const
{
    double eta = (*this.*selectedDamageEnhancement)(dose);
    
    if (dose <= D_t3)
        return - alpha_X3 * dose - beta_X3 * dose*dose;
    else
        return logS_t3 - s3 * (eta * dose - D_t3);
}

//------------------------------------------------------------------------------

double CellLine::readDamageEnhancement(const double dose) const
{
    double eta, getD = 1e2;
    //double putD = 1e2;
    
    char *dataPath;
    dataPath = getenv("DATA");
    if( dataPath == NULL )
    {
        cerr << "DATA environment variable could not be determined." << endl;
        exit(1);
    }
    
    string etaFile(dataPath);
    etaFile += "/eta_Elsasser2007.dat";
    
    ifstream inEtaFile(etaFile.c_str(), ios::in);
    if (!inEtaFile)
    {
        cout << "File " << etaFile << " could not be opened." << endl;
        exit(1);
    }
    
    do {
        inEtaFile >> getD >> eta;
    } while( getD < dose );
    inEtaFile.close();
    
    return eta;
}

//------------------------------------------------------------------------------

void CellLine::setParametrization(const string parametrization_type)
{
    if(parametrization_type == "LQ_noDt" && isLQ_noDtLoaded)
        selectedParametrization = &CellLine::parametrization_LQ_noDt;
    
    else if(parametrization_type == "LQ_noDt_T" && isLQ_noDt_TLoaded)
        selectedParametrizationT = &CellLine::parametrization_LQ_noDt_T;
    
    else if(parametrization_type == "LQ" && isLQloaded)
        selectedParametrization = &CellLine::parametrization_LQ;
    
    else if((parametrization_type == "LQ2" || parametrization_type == "LQ2_interpolated_readfile") && isLQ2loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ2;
        selectedDamageEnhancement = &CellLine::interpolatedDamageEnhancement;
        selectedEtaGeneration = &CellLine::readDamageEnhancement;
        needEtaGenerated = true;
    }
    
    else if(parametrization_type == "LQ3" && isLQ3loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ3;
        selectedDamageEnhancement = &CellLine::interpolatedDamageEnhancement;
        selectedEtaGeneration = &CellLine::readDamageEnhancement;
        needEtaGenerated = true;
    }
    
    else if(parametrization_type == "LQ2_interpolated_analytic" && isLQ2loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ2;
        selectedDamageEnhancement = &CellLine::interpolatedDamageEnhancement;
        selectedEtaGeneration = &CellLine::analyticDamageEnhancement;
        needEtaGenerated = true;
    }
    
    else if(parametrization_type == "LQ2_interpolated_MC" && isLQ2loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ2;
        selectedDamageEnhancement = &CellLine::interpolatedDamageEnhancement;
        selectedEtaGeneration = &CellLine::damageEnhancement;
        needEtaGenerated = true;
    }
    
    else if(parametrization_type == "LQ2_punctual_analytic" && isLQ2loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ2;
        selectedDamageEnhancement = &CellLine::analyticDamageEnhancement;
    }
    
    else if(parametrization_type == "LQ2_punctual_MC" && isLQ2loaded)
    {
        selectedParametrization = &CellLine::parametrization_LQ2;
        selectedDamageEnhancement = &CellLine::damageEnhancement;
    }
    
    else
    {
        cerr << "CellLine: The selected X-ray parametrization is not provided or has not been loaded." << endl;
        exit(1);
    }

    //clog << "Selected " << parametrization_type << " parametrization for X-rays." << endl;

    if(needEtaGenerated)
    {
        cout << "Generating damage enhancement... ";
        double d = 100;
        int etaCounter = 0;
        int points = 200;
        do {
            doseForEta[etaCounter] = d;
            etaPre[etaCounter] = (*this.*selectedEtaGeneration)(d);
            etaCounter++;
            d = pow( double(10), 2 + etaCounter * (4 + log10(5))/(points - 1) );
        } while( d <= 5e6 );
        needEtaGenerated = false;
        cout << "done." << endl
             << endl << flush;
    }
}
