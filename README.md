Survival
========


This software was developed by the INFN (Istituto Nazionale di Fisica Nucleare) in collaboration with the University of Torino (UniTO, Physics Department) and provides different implementations of some radiobiological models to predict the cell survival after irradiation with ion beams. The implemented models are (for the moment): LEMI, LEMII, LEMIII, MKM and MCt-MKM.

The code is written in C++ and makes use of the GSL (GNU Scientific Libraries) and OpenMP (Open Multi-Processing) external libraries.

### Documentation
 - A full description of the code can be found in the [User's Manual](https://github.com/batuff/Survival/tree/master/Documentation/Survival_USERS_MANUAL.pdf).
 - Once you have downloaded the *Survival* folder, an interactive html manual is available opening the `../Documentation/html/index.html` file.

### Licensing
The *Survival* code is distributed under terms of the [GNU General Public Licence](https://github.com/batuff/Survival/edit/master/LICENSE)

### Reference paper:
L. Manganaro, G. Russo, F. Bourhaleb, F. Fausti, S. Giordanengo, V. Monaco, R. Sacchi, A. Vignati, R. Cirio and A. Attili. "Survival": a simulation toolkit introducing a modular approach for radiobiological evaluations in ion beam therapy. *Physics in Medicine and Biology* **63**(8), 08NT01 (2018). https://doi.org/10.1088/1361-6560/aab697

### A selection of papers describing the implemented models
A detailed description of the implemented models can be found in the following papers:
 
 - Local Effect Model (LEM):
 
      - LEM I: M. Scholz and G. Kraft, "Track structure and the calculation of biological effects of heavy charged particles", *Advances in Space Research* **18**, 5-14 (1996).
   
      - LEM II: T. Elsässer and M. Scholz, "Cluster effects within the local effect model", *Radiation Research* **167**, 319-329 (2007).
 
      - LEM III: T. Elsässer, M. Krämer and M. Scholz, "Accuracy of the local effect model for the prediction of biologic effects of carbon ion beams *in vitro* and *in vivo*", *International Journal of Radiation Oncology-Biology-Physics* **71**, 866-872 (2008).
      
      - LEM I,II,III rapid calculation (GSI approach): M. Krämer and M. Scholz, “Rapid calculation of biological effects in ion radiotherapy”, *Physics in Medicine and Biology* **51**, 1959-1970 (2006).

      - LEM I,II,III rapid calculation (INFN approach): G. Russo, "Develpment of a radiobiological database for carbon ion Treatment Planning Systems - Modelling and simulating the irradiation process", *Ph.D. Thesis*, Università degli studi di Torino (2011).
      
 
 - Microdosimetric Kinetic Model (MKM):
 
      - Original MKM formulation: R.B. Hawkins, "A Statistical Theory of Cell Killing by Radiation of Varying Linear Energy Transfer", *Radiation Research* **140**, 366-374 (1994).

      - MKM + amorphous track structure: Y. Kase, T. Kanai, N. Matsufuji, Y. Furusawa, T. Elsasser, and M. Scholz, “Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation”, *Physics in Medicine and Biology* **53**, 37-59 (2008).

      - MCt-MKM: Manganaro, L., Russo, G., Cirio, R., Dalmasso, F., Giordanengo, S., Monaco, V., … Attili, A. (2017). A Monte Carlo approach to the microdosimetric kinetic model to account for dose rate time structure effects in ion beam therapy with application in treatment planning simulations. *Medical Physics*, **44**(4), 1577–1589.


### Usage (Unix and Unix-like systems)

#### Setting the environment

In order to run the program, the user has to set a couple of environment variables. The BASH script `setenv.sh` automates this procedure. In the script the user has to define the variable `$install_folder` to the install folder where the program is located.

#### Running the program

The user has to call from the command line:

      source setenv.sh

to set the system variables and then

      survival -SIMULATION_OPTION CHOSEN_VALUES ...

to execute the program.

The user has the possibility to set a number of physical (and not only physical) parameters by using the syntax: `-PARAMETER_NAME PARAMETER_VALUE`.
 
Typing `survival --help` a simple help text will be displayed, suggesting how to use the program. In the following sections a complete list of parameters and their meaning is reported.


#### Output options

 - `-projectName` It's a string representing the prefix to give at any file and directories that will be created in the simulation. The default value is "NewProject".

 - `-output` The user has the possibility to choose between three kinds of output (and all possible combination between them):
 
      1. "LQ_pars" Then a file will be created, named "PROJECTNAME_LQparameters_MKM.csv", containing the information about the parameters chosen for the simulation and the values of the simulated LQ &#945; and &#946; parameters (a new line for each energy evaluated).
   
      2. "meanValues" Then a file will be created, named "PROJECTNAME_survival_MKM.csv", containing the information about the parameters chosen for the simulation and the values of doses delivered and survival observed (a new line for each energy or dose evaluated).
   
      3. "cellValues" This kind of output is supported only by the MonteCarlo calculusType. It is a way to store the values of dose and survival obtained for each single cell irradiated during the monte carlo simulation. Then a directory will be created, named "PROJECTNAME_survival_data". In the directory the user will find a description file named "000_MonteCarlo_parameters.csv", listing the parameters used in the simulation, and a directory with the same name containing the corresponding data. In particular in the subdirectory some file will be created (a file for each level of dose imposed), each one containing two column with the dose delivered and the survival observed for each cell irradiated. When a new simulation is lauched with the same project name, the program will do a check over all the description files in the directory, if the parameters of the simulation are the same of another one already done, then it will enter the related subdirectory and append data there, if not a new description file (with progressive number) and corresponding subdirectory will be created.
      
      **Warning**: these informations could occupy a lot of memory in the computer, depending on the parameter set for the simulation. Use with caution!
      
      **Note**: This parameter has to be specified. No default values are set.
      
   - `-cellType` A string identifying the name of the cell lline used in the calculation. The default value is "Cell1".

   **Note**: the cell line is only determined by the chosen model parameters. `-cellType` is only a tag to indicate the cell useful for bookeeping, but it isn't used in the simulation.




#### Model Selection

   - `-model` It's a string representing the model to use in the simulation. Five models are supported (see the references in the section named *"Papers describing the implemented models"*):
   
      1. "LEMI"
      
      2. "LEMII"
      
      3. "LEMIII"
      
      4. "MKM"
      
      5. "tMKM_Manganaro2017"

   The default value for this option is "MKM"

##### LEM parameters

 - `-LEM_alpha0` A double representing associated to the linear quadratic &#945; parameter characteristic for X-rays, expressed in Gy<sup>-1</sup>. The simulated &#945; parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.312 Gy<sup>-1</sup>, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_beta0` A double representing associated to the linear quadratic &#946; parameter characteristic for X-rays, expressed in Gy<sup>-2</sup>. The simulated &#946; parameter will tend to this value for low LET. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 0.073 Gy<sup>-2</sup>, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_rNucleus` The radius of the cell expressed in &#956;m. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 4.611 &#956;m, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-LEM_Dt` The transition dose beyond which the standard linear quadratic parametrization is no more valid, expressed in Gy. This option is compatible only with the LEM (I, II and III) model, it won't be use if a different model is chosen. The default value is 30 Gy.
   
##### MKM parameters

   - `-MKM_alpha0` A double representing associated to the linear quadratic &#945; parameter characteristic for X-rays, expressed in Gy<sup>-1</sup>. The simulated &#945; parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.312 Gy<sup>-1</sup>, a tipical value representing the Human Salivary Gland (HSG) cell line.

   - `-MKM_beta0` A double representing associated to the linear quadratic &#946; parameter characteristic for X-rays, expressed in Gy<sup>-2</sup>. The simulated &#946; parameter will tend to this value for low LET. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.073 Gy<sup>-2</sup>, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_rNucleus` The radius of the cell expressed in &#956;m. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 4.611 &#956;m, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_rDomain` The radius of domains wich constitute the MKM nucleus expressed in &#956;m. This option is compatible only with the MKM (and tMKM) model, it won't be use if a different model is chosen. The default value is 0.365 &#956;m, a tipical value representing the Human Salivary Gland (HSG) cell line.

 - `-MKM_timeConst`: The time constant associated to the repair kinetics of the cell. This parameters is only used in the "tMKM_Manganaro2017" model.
      

#### Evaluation method options
 
   - `-calculusType` A string identifying the type of numerical approach ("calculus") to be used in the simulations. Different possibilities are available:
 
      1. "rapidLEM_Scholz2006" It's an implementation of the method described in: M. Krämer and M. Scholz, "Rapid calculation of biological effects in ion radiotherapy", *Physics in medicine and biology* **51**, 1959-1970 (2006). It's compatible only with `LEMI`, `LEMII` and `LEMIII` models.
   
      2. "rapidLEM_Russo2011" A rapid evaluation method for LEM, with better agreement to the Monte Carlo approach, described in: G. Russo, "Develpment of a radiobiological database for carbon ion Treatment Planning Systems - Modelling and simulating the irradiation process", *Ph.D. Thesis*, Università degli studi di Torino (2011). It's compatible only with `LEMI`, `LEMII` and `LEMIII` models.
   
      3. "rapidMKM_Kase2008" A fast implementation of the MKM calculation as described in: Kase, Y., Kanai, T., Matsufuji, N., Furusawa, Y., Elsässer, T., & Scholz, M. (2008). Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation. *Physics in Medicine and Biology*, **53**(1), 37–59. It's compatible only with the `MKM`.
   
      4. "rapidMKM_Kase2008_corrected_beta" An extension of the "rapidMKM_Kase2008" method in which a non-Poissonian correction factor is added also for the &#946; parameter.  It's compatible only with the `MKM`.
      
      5. "rapidMKM_Attili2013" A fast original implementation of the MKM model, combining the methods described in: Hawkins_2003 Hawkins, R. B. (2003). A microdosimetric-kinetic model for the effect of non-Poisson distribution of lethal lesions on the variation of RBE with LET. *Radiation Research*, 160(1), 61–69, and Kase, Y., Kanai, T., Matsufuji, N., Furusawa, Y., Elsässer, T., & Scholz, M. (2008). Biophysical calculation of cell survival probabilities using amorphous track structure models for heavy-ion irradiation. *Physics in Medicine and Biology*, **53**(1), 37–59.  It's compatible only with the `MKM`.
   
      6. "rapidMKM_Attili2013_corrected_beta" An extension of the "rapidMKM_Attili2013" method in which a non-Poissonian correction factor is added also for the &#946; parameter. See Survival::Calculus::rapidMKM_Attili2017_corrected_beta() for details.  It's compatible only with the `MKM`.
   
      7. "MonteCarlo" Compatible with all the implemented models, it performs a Monte Carlo simulation of the irradiation process.
      
   **Note**: the model `tMKM_Manganaro2017` accepts only the Monte Carlo ("MonteCarlo") type of calculus.


##### Monte Carlo Options

   - `-precision` Supported only by the MonteCarlo calculusType, it's a `double` identifying the precision to be reached in the calculation. Two possibilities are provided, as the user can indicate:
 
      1. A positive integer: it will be taken as the number of cell to irradiate for each level of nominal dose imposed.
      
      2. A double between 0 and 1: it indicates the relative error on the simulated survival to be reached before stopping the calculation.

   - `-parallelismType` A positive integer indicating the level of parallelism to be used in the simulation. The user has the possibility to specify the number of threads to dedicate at the calculation, in particular:
   
      1. 1: Only 1 thread = parallelism disabled.
      
      2. N>1: The number of threads to be used.
      
      3. 0: Then the program will define a number of thread corresponding to the number of core of the computer executing the program. This is set as default value.


   
#### Source definition

 - `-ion` A string identifying the ion irradiated to the cells, by means of the chemical symbol of the element, without mass number specifications. Ions from proton to neon are supported.

 - `-energies` A sequence of kinetic energies of the ion to be evaluated, expressed in MeV. The energy values should be separated by one o more spaces. This and the `-lets` option are mutually exclusive, but one of the two has to be specified.

 - `-lets` A sequence of LET of the ion to be evaluated, expressed in keV/&#956;m. The LET values should be separated by one o more spaces. This and the `-energies` option are mutually exclusive, but one of the two has to be specified.

 - `-doses` The sequence of doses to be delivered. The LET values should be separated by one o more spaces. The default is 1 to 6 Gy ("1 2 3 4 5 6"), in order to construct a survival curve.

 - `-nFraction` Supported only by the "tMKM_Manganaro2017" model, it require an integer representing the number of fraction in which to divide each nominal dose to be delivered. The default is 1 (single fraction).

 - `-timeSpacing` Supported only by the "tMKM_Manganaro2017" model, it indicate the time spacing between consecutive fractions, expressed in hours. The default is 0 (no time spacing).

 - `-fracDeliveryTime` Supported only by the "tMKM_Manganaro2017" model, it indicate the fraction delivery time, expressed in hours. The default is 0 (istantaneous delivering).

 - `-spectrum_file` Useful for mixed fields evaluation. This option provide the possibility to irradiate the cellular population with different ions and different energies in a mixed field, described in an external file (see the TEMPLATE for the spectrum_file). The user has to specify the name of the file containing the spectrum complete with its relative or absolute path.
 
 **Warning**: The `-spectrum_file` option hasn't been deeply tested!

   - `-trackMode` In the case of mixed fields this option require a string indicating the way in which to interpretate the spectrum specification. Two possibilities are provided:
 
      1. "histogram": then the spectrum will be considered as an histogram where each particle is a bin with its "weight" (see the TEMPLATE for the spectrum_file).
   
      2. "random": then the program will random extract with uniform probability (iteration by iteration) the particle to use.
