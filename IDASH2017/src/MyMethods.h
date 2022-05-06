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
