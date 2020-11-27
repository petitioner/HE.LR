/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef IDASH2017_MYMETHODS_H_
#define IDASH2017_MYMETHODS_H_
#include "MyTools.h"
#include "Scheme.h"

class MyMethods {
public:
	
	static double* testCryptoFullBatchNAGwithG(double** traindata, double* trainlabel, long factorDim, long trainSampleDim, long numIter, double** testdata, double* testlabel, long testSampleDim, string resultpath);
	static double* testCryptoMiniBatchNAGwithG(double** traindata, double* trainlabel, long factorDim, long trainSampleDim, long numIter, double** testdata, double* testlabel, long testSampleDim, string resultpath);	

};

#endif /* MYMETHODS_H_ */
