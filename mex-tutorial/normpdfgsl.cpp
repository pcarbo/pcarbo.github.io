#include "mex.h"
#include "gsl/gsl_randist.h"

using namespace std;

extern void _main();

const int numInputArgs  = 3;
const int numOutputArgs = 1;

// Function declarations.
// -----------------------------------------------------------------
double  getMatlabScalar    (const mxArray* ptr);
double& createMatlabScalar (mxArray*& ptr);

// Function definitions.
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[]) {
  
  // Check to see if we have the correct number of input and output
  // arguments.
  if (nrhs != numInputArgs)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != numOutputArgs)
    mexErrMsgTxt("Incorrect number of output arguments");

  // Get the inputs.
  double x  = getMatlabScalar(prhs[0]);
  double mu = getMatlabScalar(prhs[1]);
  double v  = getMatlabScalar(prhs[2]);

  // Create the output. It is also a double-precision scalar.
  double& p = createMatlabScalar(plhs[0]);

  // Compute the value of the univariate Normal at x.
  p = gsl_ran_gaussian_pdf(x-mu,sqrt(v));
}

double getMatlabScalar (const mxArray* ptr) {

  // Make sure the input argument is a scalar in double-precision.
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("The input argument must be a double-precision scalar");

  return *mxGetPr(ptr);
}

double& createMatlabScalar (mxArray*& ptr) { 
  ptr = mxCreateDoubleMatrix(1,1,mxREAL);
  return *mxGetPr(ptr);
}
