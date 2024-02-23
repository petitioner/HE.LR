#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>

#include "MyMethods.h"
#include "MyTools.h"

using namespace std;
using namespace NTL;


int main(int argc, char **argv) {

	SetNumThreads(8);

	long numIter = 35;

	//string trainfile = "../data/Credit_train.csv";
	//string testfile = "../data/Credit_test.csv";
	string trainfile = "../data/MNIST_train.txt";
	string testfile = "../data/MNIST_test.txt";

	long trainSampleDim = 0, testSampleDim = 0, trainfactorDim = 0,	testfactorDim = 0;
	double **traindataset, **testdataset;
	double *traindatalabel, *testdatalabel;

	double **zData = MyTools::dataFromFile(trainfile, trainfactorDim, trainSampleDim, traindataset, traindatalabel);
	double **zDate = MyTools::dataFromFile(testfile, testfactorDim, 	testSampleDim, testdataset, testdatalabel);
	if (trainfactorDim != testfactorDim) {
		cout << " WARNING : trainfactorDim != testfactorDim" << endl;
		exit(-1);
	}

	MyTools::normalizeZData(traindataset, trainfactorDim, trainSampleDim);
	MyTools::normalizeZData(testdataset, testfactorDim, testSampleDim);


	string pathNesterovAGwithXTXasG =     "../result/20201123_FGCS";

		   //pathNesterovAGwithXTXasG.append("_MNIST");
		   pathNesterovAGwithXTXasG.append("_Credit");
		   
	//pathNesterovAGwithXTXasG.append("_MiniBatch_");
	pathNesterovAGwithXTXasG.append("_FullBatch_");

	// Step 1. clear the former result data stored in the four*3(AUC,MLE,TIME) different files.
	// SHOULD BE IN CONSSITENT WITH THE FILE PATH IN THE EACH ALGORITHM ! eg. "TrainAUC.csv"...
	std::ofstream ofs;

	ofs.open(pathNesterovAGwithXTXasG + "TrainAUC.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "TrainMLE.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "TestAUC.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "TestMLE.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "TIME.csv",			std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "TIMELabel.csv",	std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "CurrMEM.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();
	ofs.open(pathNesterovAGwithXTXasG + "PeakMEM.csv",		std::ofstream::out | std::ofstream::trunc);	ofs.close();

	// Step 2. Test Training and Testing ONLY one Time!
	//         each time we choose 1579 samples from the training data set at random
	//         use the whole test data set to validate the performance of module
	long testSampleNum = testSampleDim;
	long trainSampleNum = trainSampleDim;

	double **traindata, **testdata;
	double *trainlabel, *testlabel;

	traindata = new double*[trainSampleNum];
	trainlabel = new double[trainSampleNum];

	testdata = testdataset;
	testlabel = testdatalabel;



		// randomly choose 1579 samples from the train data set
		MyTools::shuffleDataSync(traindataset, trainfactorDim, trainSampleDim, traindatalabel);
		for (long i = 0; i < trainSampleDim; ++i) {
			traindata[i] = traindataset[i];
			trainlabel[i] = traindatalabel[i];
		}

		cout << "!!! Crypto: Full-Batch NAG With 0.25XTX as Quadratic Gradient !!!" << endl;
		MyMethods::testCryptoFullBatchNAGwithG(traindata, trainlabel, trainfactorDim,
				trainSampleDim, numIter, testdata, testlabel, testSampleDim,	pathNesterovAGwithXTXasG);
		cout << "@@@ Crypto: Full-Batch NAG With 0.25XTX as Quadratic Gradient @@@" << endl;

		// cout << "### Crypto: Mini-Batch NAG With 0.25XTX as Quadratic Gradient ###" << endl;
		// MyMethods::testCryptoMiniBatchNAGwithG(traindata, trainlabel, trainfactorDim,
		// 	trainSampleDim, numIter, testdata, testlabel, testSampleDim,	pathNesterovAGwithXTXasG);
		// cout << "$$$ Crypto: Mini-Batch NAG With 0.25XTX as Quadratic Gradient $$$" << endl;


	cout << endl << "END OF THE PROGRAMM" << endl;
	return 0;
}
