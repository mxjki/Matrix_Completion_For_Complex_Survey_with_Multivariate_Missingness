#Date: 2019-12-29

Description: this contains the implementations of the proposed method and other compared methods in the simulation study.

########################################
Demo_first_simulation.R: setup for the first simulation study*.
########################################

########################################
SMC_functions.R: Auxiliary functions*.
########################################
cv_validate_sim:    This function calculates the validation error to tuning 				 parameters.

SVTE_alpha:        This function performs the singular value soft-thresholding                              			and scaling procedures.

SMCfit_cv:          This function performs the matrix completion part in the 				proposed method with the tunning parameter choosed by 				cross-validation.

SMCfit:            This function performs the matrix completion part in 				 proposed method with fix tunning parameters.

sample.gen:        This function is used to generate a sample.  

sample.miss:       	This function is used to generate the missingness for the 				 sample.

var.double.robust: 	This function is used to apply the plug-in variance 					 estimator.   

est.double.r: 	 This function is used to apply the doubly robust 					 estimator.  

one.iter:		 This is one iteration of the proposed method. 


*: More details, including comments, the arguments and outputs, are provided in the corresponding file. 