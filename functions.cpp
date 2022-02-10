/*
***********************************************************************************************************************************************************

												Levenberg-Marquardt Algorithm
													Fit Routine Functions

Author: Omkar Junnarkar, Semester-3 MSc. Computational Material Science
Matriculation Nr.: 66157	Email: omkar.junnarkar@student.tu-freiberg.de
IDE : Microsoft Visual Studio 2019

iostream: For Input-Output of the C++ Program
Eigen/Dense: For Using Dense Matrix Interface of Linear Algebra Library Eigen
iomanip: To manipulate the number of decimal places in output
fstream: To create the stream of file and write
functions.h: Contains the Fitting Routine Header

*/

#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<math.h>
#include<fstream>
#include"functions.h"

/*
To reduce the effort of specifying libraries/class for each command/data type
*/
using namespace std;
using namespace Eigen;

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	Arguments : Parameters of Equation, Input(x) for Equation 
	para( , ) : List of Parameters
	x( , ) : List of Input values for the equation

	Output : Matrix 'y', data points of the funcction y=f(parameters,x) 
*/
MatrixXd function_y(MatrixXd para, MatrixXd x) {
	int number_of_data_points = x.rows();
	MatrixXd y(number_of_data_points, 1);
	for (int i = 0; i < number_of_data_points; i++) {

		//y(i, 0) = para(0, 0) * pow(x(i, 0), 2) + para(1, 0) * exp(pow(z(i, 0), 2) + 1) + para(2, 0);
		//y(i, 0) = para(0, 0) * pow(x(i, 0), 4) + para(1, 0) * exp(para(2, 0) * x(i, 0));

		y(i, 0) = para(0, 0) * pow(x(i, 0), 3) + para(1, 0) * pow(x(i, 0), 2) + para(2, 0) * x(i, 0) + para(3, 0) * exp(para(4, 0) * x(i, 0) + para(5, 0));
	}
	return y;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
/*
	Jacobian Matrix : the matrix containing values of First Order Partial Derivatives w.r.t all parameters at All Data Points
	Size : (Number of Data Points, Parameters)
	The Partial Derivatives of the functions with respect to each Parameter can be estimated by using Finite Difference Method as follows :

	Df/Da = [f(a + DELa) - f(a)]/DELa , Df/Db = [f(b + DELb) - f(b)]/DELb , Df/Dc = [f(c + DELc) - f(c)]/DELc  
	where a,b,c are paramerters and DELa,DELb,DELc are the deflections (given by initial deflection)

	Thus by computing function values of these derivates at all data points 'x', All Elements of Jacobian can be obtained.
	[ Refer Report/Manual for Details ]
	y_deflected : Function value obtained by deflecting parameter
	
	Arguments: Estimated Parameters, Ínitial Deflection, 'y' measured, Input data 'x'
	Output : Jacobian Matrix
*/
MatrixXd getJacobianMatrix(MatrixXd para_est, MatrixXd deflection, MatrixXd ym, MatrixXd input) {

	
	MatrixXd Jacobian_Matrix(ym.rows(), para_est.rows());
	MatrixXd y = function_y(para_est, input);
	MatrixXd y_deflected(ym.rows(), 1);

	for (int i = 0; i < para_est.rows(); i++) {					/* Iterating through all parameters */

		para_est(i, 0) = para_est(i, 0) + deflection(i, 0);		/* Changing the parameters one by one */

		y_deflected = function_y(para_est, input);				/* Computing the deflected function arrray */
		for (int j = 0; j < input.rows(); j++) {

			// [f(v, p + dp) - f(v, p) ] / [dp] 

			Jacobian_Matrix(j, i) = (y_deflected(j, 0) - y(j, 0)) / deflection(i, 0);
		}
		para_est(i, 0) = para_est(i, 0) - deflection(i, 0);		/* Bringing back the parametes to original value */
	}

	return Jacobian_Matrix;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

/*
	Levenberg-Marquardt Algorithm (Damped Least Squares Method) :
	
	Arguments : Estimated Parameters, Initial Deflection of Parameters, Measured 'y' values, Data Input 'x' 
	Output : True Parameters

	Computation Strategy & Variables:

	-> Hessian Matrix (H) and difference between measured y values and estimated y values is computed (d) along with error
	-> Hessian for LMA (H_lm) is computed by Hessian Matrix + Lambda*IdentitiyMatrix
	-> Change in Parameters (dp) is computed by H_lm.inverse*J.transpose*d -> Parameters are updated and new error is computed
	-> If error is reduced, then Lambda(Damping Factor) is reduced by a factor of 10, because solution moves closer to the real solution 
	-> Else Lambda (Trust Regin Radius) is increased and parameters are not updated and Iterations continue untill condition is fulfilled

	Lambda has influence on the behaviour of routine. Initial value is set to 10. In case of complex curvy functions, Lambda can be Increased.
	[ Refer Report/Manual for Details ]

*/
MatrixXd LevenbergMarquardtFit(MatrixXd para_guess, MatrixXd deflection, MatrixXd ym, MatrixXd input) {

	cout << "-> Entered LevenbergMarquardtFit\n";

	ofstream errorfile("ErrorNorm.csv");
	ofstream hessianfile("Hessians.csv");
	ofstream lambdafile("DampingFactorPropagation.csv");

	/* Number of Parameters */
	int npara = para_guess.rows(), ndata = ym.rows();

	MatrixXd IdentityMat = MatrixXd::Identity(npara, npara);			
	MatrixXd H(npara, npara);	
	MatrixXd d(ndata, 1);
	MatrixXd J(ndata, npara);
	double error,error_lm;

	MatrixXd y_init = function_y(para_guess, input);
	
	double lambda = 10;					/* Damping Factor */
	int updateJ = 1;					/* Flag */
	MatrixXd para_est = para_guess;		/* Estimated Parameters*/
	int maxiter = 100, counter=0;		/* Maximum Iterations Allowed, Counter for Iterations*/
	

	while (counter < maxiter) {
		
		cout << "--> Iteration : " << counter << endl;

		if (updateJ == 1) {
			
			J = getJacobianMatrix(para_est, deflection, ym,input);
			
			MatrixXd y_est = function_y(para_est, input);
			d = ym - y_est;
			
			H = J.transpose() * J;
			cout << "H: \n" << H << endl;
			
			if (counter == 0) {							/* Only to be done for First Iteration */
				MatrixXd temp1 = d.transpose() * d;
				error = temp1(0,0);
				
			}
		}

		MatrixXd H_lm = H + lambda * IdentityMat;			/* Updating Hessian */
		MatrixXd dp = H_lm.inverse() * J.transpose() * d;	/*Computing change in parameters*/
		//cout << "dp: \n" << dp;
		MatrixXd para_lm = para_est + dp;					/* temperory update of parameters*/
		MatrixXd y_est_lm = function_y(para_lm, input);		/* temperory 'y' values from updated parameters*/
		MatrixXd d_lm = ym - y_est_lm;						/* difference betwen measured and computed values*/
		MatrixXd temp2 = d_lm.transpose() * d_lm;			/* Present Error */
		error_lm = temp2(0,0);

		if (error_lm < error) {					

			/* Accepting step, Parameters and reducing trust region radius */
			lambda = lambda / 10;
			para_est = para_lm;
			error = error_lm;
			updateJ = 1;
		}
		else {

			/* Rejecting step and temperory update of parameters, Increasing trust region radius */
			updateJ = 0;
			lambda = lambda * 10;
		}

		/* Stopping Criterion */
		if (dp.norm() < 1e-4) {
			counter = maxiter;
		}
		else counter++;
		
		/* Writing Files */
		hessianfile << H << endl << endl;
		errorfile << error_lm << endl;
		lambdafile << lambda << endl;

	}

	hessianfile.close();
	errorfile.close();
	lambdafile.close();

	/* Returning Final Set of Obtained Parameters */

	return para_est;
};

/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
