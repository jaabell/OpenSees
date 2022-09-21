#include <iostream>
#include "LinearHardeningScalar_EV.h"

#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

int main(int argc, char const *argv[])
{
	
	LinearHardeningScalar_EV k(1.0, 1.0);

	cout << "k = " << k << endl;

	VoightTensor6 depsilon 
		= VoightTensor6::Zero();//{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	VoightTensor6 m 
		= {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	VoightTensor6 sigma 
		= VoightTensor6::Zero();//{0.0,, 0.0,, 0.0,, 0.0,, 0.0,, 0.0};


	double h = k.getDerivative(depsilon, m, sigma);



	cout << "depsilon = \n" << depsilon << endl;
	cout << "m = \n" << m << endl;
	cout << "sigma = \n" << sigma << endl;

	cout << "h = " << h << endl;



	return 0;
}