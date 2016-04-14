Survival
========


This software was developed by the INFN (Istituto Nazionale di Fisica Nucleare) in collaboration with the University of Torino (UniTO, Physics Department) and provides different implementations of some radiobiological models to predict the cell survival after irradiation. The implemented models are (for the moment): LEMI, LEMII, LEMIII, MKM and MCt-MKM.

The code is written in C++ and makes use of the GSL (GNU Scientific Libraries) and OpenMP (Open Multi-Processing) external libraries.

### Full description
A full description of the code can be found in the User's Manual: "Documentation/Survival_REFERENCE_MANUAL.pdf".

### Reference published paper:
*Writing in progress!!*

### Papers describing the implemented models
A detailed description of the implemented models can be found in the following papers:
 - LEM I: M. Scholz and G. Kraft, "Track structure and the calculation of biological effects of heavy charged particles", *Advances in Space Research* **18**, 5-14 (1996).
 - LEM I rapid calculation: M. Krämer and M. Scholz, "Rapid calculation of biological effects in ion radiotherapy", *Physics in medicine and biology* **51**, 1959-1970 (2006).
 - LEM II: T. Elsässer and M. Scholz, "Cluster effects within the local effect model", *Radiation Research* **167**, 319-329 (2007).
 - LEM III: T. Elsässer, M. Krämer and M. Scholz, "Accuracy of the local effect model for the prediction of biologic effects of carbon ion beams *in vitro* and *in vivo*", *International Journal of Radiation Oncology-Biology-Physics* **71**, 866-872 (2008).
 - MKM:
   1. Original MKM formulation: R.B. Hawkins, "A Statistical Theory of Cell Killing by Radiation of Varying Linear Energy Transfer", *Radiation Research* **140**, 366-374 (1994). [Some corrections and improvements were made over the subsequent years].
   2. Kiefer-Chatterjee amorphous track structure introduction: Y. Kase, T. Kanai, N. Matsufuji, Y. Furusawa, T. Elsasser, and M. Scholz, "Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation", *Physics in Medicine and Biology* **53**, 37-59 (2008).
 - MCt-MKM: L. Manganaro, G. Russo, R. Cirio, F. Dalmasso, S. Giordanengo, V. Monaco, R. Sacchi, A. Vignati and A. Attili, "A novel formulation of the Microdosimetric Kinetic Model to account for dose-delivery time structure effects in ion beam therapy with application in treatment planning simulations", *Medical Physics*, **Submitted**.

### Usage (Unix and Unix-like systems)
To execute the program the user has to call from the command line:  
`$ source setenv.sh`  
`$ ./survival -SIMULATION_OPTION CHOSEN_VALUES ...`  
Hence, the user has the possibility to set a number of physical (and not only physical) parameters by using the syntax:  
`-PARAMETER_NAME PARAMETER_VALUE`.
 
Here is the complete list of parameters and their meaning:
 - `-projectName` It's a string representing the prefix to give at any file and directories that will be created in the simulation. The default value is "NewProject".

 - `-output` The user has the possibility to choose between three kinds of output (and all possible combination between them):
   1. "LQ_pars" Then a file will be created, named "PROJECTNAME_LQparameters_MKM.csv", containing the information about the parameters chosen for the simulation and the values of the simulated LQ \f$\alpha\f$ and \f$\beta\f$ parameters (a new line for each energy evaluated).
   2. "meanValues" Then a file will be created, named "PROJECTNAME_survival_MKM.csv", containing the information about the parameters chosen for the simulation and the values of doses delivered and survival observed (a new line for each energy or dose evaluated).
   3. "cellValues" This kind of output is supported only by the MonteCarlo calculusType. It is a way to store the values of dose and survival obtained for each single cell irradiated during the monte carlo simulation. Then a directory will be created, named "PROJECTNAME_survival_data". In the directory the user will find a description file named "000_MonteCarlo_parameters.csv", listing the parameters used in the simulation, and a directory with the same name containing the corresponding data. In particular in the subdirectory some file will be created (a file for each level of dose imposed), each one containing two column with the dose delivered and the survival observed for each cell irradiated. When a new simulation is lauched with the same project name, the program will do a check over all the description files in the directory, if the parameters of the simulation are the same of another one already done, then it will enter the related subdirectory and append data there, if not a new description file (with progressive number) and corresponding subdirectory will be created.
   
   **Warning**: these informations could occupy a lot of memory in the computer, depending on the parameter set for the simulation. Use with caution!

   **Note**: This parameter has to be specified. No default values are set.

 - `-precision` Supported only by the MonteCarlo calculusType, it's a `double` identifying the precision to be reached in the calculation. Two possibilities are provided, as the user can indicate:
   1. A positive integer: it will be taken as the number of cell to irradiate for each level of nominal dose imposed.
   2. A double between 0 and 1: it indicates the relative error on the simulated survival to be reached before stopping the calculation.

 - `-parallelismType` A positive integer indicating the level of parallelism to be used in the simulation. The user has the possibility to specify the number of threads to dedicate at the calculation, in particular:
   1. 1: Only 1 thread = parallelism disabled.
   2. N>1: The number of threads to be used.
   3. 0: Then the program will define a number of thread corresponding to the number of core of the computer executing the program. This is set as default value.
 - `-model` It's a string representing the model to use in the simulation. Five models are supported (see the references in the section named *"Papers describing the implemented models"*):
   1. "LEMI"
   2. "LEMII"
   3. "LEMIII"
   4. "MKM"
   5. "tMKM"

 The default value for this option is "MKM"
 - `-calculusType` A string identifying the type of calculus to be done. Some possibilities are available:
   1. "rapidScholz" It's an implementation of the method described in: M. Krämer and M. Scholz, "Rapid calculation of biological effects in ion radiotherapy", *Physics in medicine and biology* **51**, 1959-1970 (2006). It's compatible only with LEMI, LEMII and LEMIII models.
   2. "rapidRusso" A new rapid method for LEM, proved to be more accurate, described in: G. Russo, "Develpment of a radiobiological database for carbon ion Treatment Planning Systems - Modelling and simulating the irradiation process", *Ph.D. Thesis*, Università degli studi di Torino (2011). It is compatible only with LEMI-LEMII-LEMIII models.
   3. "rapidMKM" An implementation of the original MKM calculation. It's compatible only with the MKM model and it's the default value for this option.
   4. "MonteCarlo" Compatible with all models implemented, performs a monte carlo simulation of the irradiation process to get the LQ parameters.

 - `-cellType` A string identifying the name of the cell lline used in the calculation. The default value is "Cell1".

   **Note**: The cell line in reality is completely determined by the model parameters chosen. This is only a tag to indicate the cell but it isn't used in the simulation.

 - `-MKM_alpha0` A double representing IDEALLY the linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$. The simulated \f$\alpha\f$ parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.312 \f$Gy^{-1}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_beta0` A double representing IDEALLY the linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$. The simulated \f$\beta\f$ parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.073 \f$Gy^{-2}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_rNucleus` The radius of the cell expressed in \f$\mu m\f$. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 4.611 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_rDomain` The radius of domains wich constitute the MKM nucleus expressed in \f$\mu m\f$. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.365 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_timeConst`: The time constant associated to the repair kinetics of the cell.

 - `-LEM_alpha0` A double representing IDEALLY the linear quadratic \f$\alpha\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-1}\f$. The simulated \f$\alpha\f$ parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.312 \f$Gy^{-1}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_beta0` A double representing IDEALLY the linear quadratic \f$\beta\f$-parameter characteristic for X-rays, expressed in \f$Gy^{-2}\f$. The simulated \f$\beta\f$ parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.073 \f$Gy^{-2}\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_rNucleus` The radius of the cell expressed in \f$\mu m\f$. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 4.611 \f$\mu m\f$, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_Dt` The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 30 Gy.

 - `-ion` A string identifying the chemical symbol of the element, without mass number specifications. Ions from proton to neon are supported.

 - `-energies` A sequence of kinetic energies of the primary ions to be evaluated, expressed in MeV. This and the `-lets` option are mutually exclusive, but one of the two has to be specified.

 - `-lets` A sequence of LET of the primary ions to be evaluated, expressed in MeV. This and the `-energies` option are mutually exclusive, but one of the two has to be specified.

 - `-doses` The sequence of doses to be delivered. The default is 1 to 6 Gy (step 1), in order to construct a survival curve.

 - `-nFraction` Supported only by the tMKM model, it require an integer representing the number of fraction in which to divide each nominal dose to be delivered. The default is 1 (single fraction).

 - `-timeSpacing` Supported only by the tMKM model, it indicate the time spacing between consecutive fractions, expressed in hours. The default is 0 (no time spacing).

 - `-fracDeliveryTime` Supported only by the tMKM model, it indicate the fraction delivery time, expressed in hours. The default is 0 (istantaneous delivering).

 - `-spectrum_file` Useful for mixed fields evaluation. This option provide the possibility to irradiate the cellular population with different ions and different energies in a mixed field, described in an external file (see the TEMPLATE for the spectrum_file). The user has to specify the name of the file containing the spectrum complete with its relative or absolute path.
 
 **Warning**: The *"mixed fields"* option hasn't been deeply tested!

 - `-trackMode` In the case of mixed fields this option require a string indicating the way in which to interpretate the spectrum specification. Two possibilities are provided:
   1. "histogram": then the spectrum will be considered as an histogram where each particle is a bin with its "weight" (see the TEMPLATE for the spectrum_file).
   2. "random": then the program will random extract with uniform probability (iteration by iteration) the particle to use.
 
Typing `--help` a hint will be display, suggesting how to use the program.
