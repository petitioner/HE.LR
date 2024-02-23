#include "MyMethods.h"

#include "Ciphertext.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "TestScheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "MyTools.h"
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <iomanip>

#include <algorithm>    // std::shuffle
#include <vector>        // std::vector
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#include <unistd.h>



double* MyMethods::testCryptoMiniBatchNAGwithG(double** traindata, double* trainlabel, long factorDim, long trainSampleDim, long numIter, double** testdata, double* testlabel, long testSampleDim, string resultpath)
{

	long wBits = 30;                              
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kdeg = 7;
	long kBits = (long)ceil(log2(kdeg));                    

	long logN = MyTools::suggestLogN(80, logQ);  
	long slots = 1 << (logN - 1);  
	long sBits = (long)ceil(log2(slots));   
	long fdimBits = (long)ceil(log2(factorDim));          
	long sdimBits = (long)ceil(log2(trainSampleDim));       
	
	long batch = 1 << 5;         // Basically, batch is the Number of several factor dimensions.
	long bBits = (long)ceil(log2(batch));          // 2^batchBits = min( 2^logN / 2^sdimBits / 2, 2^fdimBits ) = min( N/2 /n, factorDim ) ;
	
	// the size of batch should be small than fatctorDim
	if (batch >= factorDim) {cout << "batch >= factorDim!" << endl; exit(0);}

	//long factorNum = 1 << (long)ceil(log2(factorDim));
	long cnum = (long)ceil(double(factorDim)/batch);  // To Divide the whole Train Data into Several Batches (cnum Ciphertexts).
	if( cnum > (1 << fdimBits) ) { cout << "cnum should be no more than factorDim!" << endl; exit(0);}

	long minbatchsize = slots / batch; // (each min-batch should be small to save the useless space of the last batch)
	long minBatchDimBits = (long)ceil(log2(minbatchsize));
	  
	long rnum = (long)ceil((double)trainSampleDim / minbatchsize);  // rnum : the number of min-batch, namely number of row
	
	cout << "logQ = " << logQ << ", logN = " << logN << endl;
	cout << "slots = " << slots << ", batch = " << batch << endl;
	cout << "rnum = " << rnum << ", cnum = " << cnum << endl;
	cout << "min-batch size = " << minbatchsize << ", minBatchDimBits = " << minBatchDimBits << endl;
	



MyTools::printData(traindata, factorDim, trainSampleDim);
double** Binv = MyTools::zInvBFromFile(traindata, factorDim, trainSampleDim);
MyTools::printData(Binv, factorDim, trainSampleDim);
cout << endl;
exit(0);



	string path = resultpath;
	ofstream openFileTrainAUC(path+"TrainAUC.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainACC(path+"TrainACC.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainMLE(path+"TrainMLE.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTestAUC(path+"TestAUC.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTestACC(path+"TestACC.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTestMLE(path+"TestMLE.csv",   std::ofstream::out | std::ofstream::app);
	ofstream openFileTIME(path+"TIME.csv", std::ofstream::out | std::ofstream::app);
	ofstream openFileTIMELabel(path+"TIMELabel.csv", std::ofstream::out | std::ofstream::app);
	ofstream openFileCurrMEM(path+"CurrMEM.csv", std::ofstream::out | std::ofstream::app);
	ofstream openFilePeakMEM(path+"PeakMEM.csv", std::ofstream::out | std::ofstream::app);

	if(!openFileTrainAUC.is_open()) cout << "Error: cannot read Train AUC file" << endl;
	if(!openFileTrainACC.is_open()) cout << "Error: cannot read Train ACC file" << endl;
	if(!openFileTrainMLE.is_open()) cout << "Error: cannot read Train MLE file" << endl;
	if(!openFileTestAUC.is_open())  cout << "Error: cannot read Test AUC file" << endl;
	if(!openFileTestACC.is_open())  cout << "Error: cannot read Test ACC file" << endl;
	if(!openFileTestMLE.is_open())  cout << "Error: cannot read Test MLE file" << endl;
	if(!openFileTIME.is_open())     cout << "Error: cannot read TIME file" << endl;
	if(!openFileTIMELabel.is_open())cout << "Error: cannot read TIME Label file" << endl;
	if(!openFileTIMELabel.is_open())cout << "Error: cannot read TIME Label file" << endl;
	if(!openFileCurrMEM.is_open())  cout << "Error: cannot read Current MEMORY file" << endl;
	if(!openFilePeakMEM.is_open())  cout << "Error: cannot read Peak MEMORY file" << endl;

	openFileTrainAUC<<"TrainAUC";	openFileTrainAUC.flush();
	openFileTrainACC<<"TrainACC";	openFileTrainACC.flush();
	openFileTrainMLE<<"TrainMLE";	openFileTrainMLE.flush();
	openFileTestAUC<<"TestAUC";	openFileTestAUC.flush();
	openFileTestACC<<"TestACC";	openFileTestACC.flush();
	openFileTestMLE<<"TestMLE";	openFileTestMLE.flush();
	openFileTIME<<"TIME";	    openFileTIME.flush();
	openFileTIMELabel<<"TIMELabel";	openFileTIMELabel.flush();
	openFileCurrMEM<<"MEMORY(GB)";	openFileCurrMEM.flush();
	openFilePeakMEM<<"MEMORY(GB)";	openFilePeakMEM.flush();


	TimeUtils timeutils;
	timeutils.start("Scheme generating");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Scheme generating";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();

	cout << "NOW THE PEAK RSS IS: " << ( MyTools::getPeakRSS() >> 20 ) << endl;


    timeutils.start("Bootstrap Key generating");
    long logT = 3; long logI = 4;
    long bootlogq = 30  +10;
    long lognslots = (long)ceil(log2(batch));  //batch = factorDim / cnum;
    scheme.addBootKey(secretKey, lognslots, bootlogq + 4);
    timeutils.stop("Bootstrap Key generated");

    openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
    openFileTIMELabel<<","<<"Bootstrap Key generating";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();

	
	cout << "NOW THE CURR RSS IS: " << ( MyTools::getCurrentRSS() >> 20 ) << endl;
	cout << "NOW THE PEAK RSS IS: " << ( MyTools::getPeakRSS() >> 20 ) << endl;
	
	// Basically, rpoly is used to calculate the sum{row}
	timeutils.start("Polynomial generating...");
	ZZ* dummy = new ZZ[N];;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	complex<double>* pvals = new complex<double> [slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j] = 1.0;
	}
	scheme.ring.encode(dummy, pvals, slots, pBits);
	delete[] pvals;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	timeutils.stop("Polynomial generation");
	
	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Polynomial generating";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();

	double* cwData = new double[factorDim]; 
	double* cvData = new double[factorDim];

	Ciphertext* encTrainData = new Ciphertext[rnum * cnum];
	Ciphertext* encTrainLabel= new Ciphertext[rnum * cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];
	//without [exit(0);], error: *** glibc detected *** ./iDASH2017: double free or corruption (!prev): 0x0000000001bbf490 ***

	/* - - - - - - - - - - - - - - - - - - - - - - - - Client and Server - - - - - - - - - - - - - - - - - - - - - - - - */
	/* cipherGD.encZData(encZData, zDataTrain, slots, factorDim, trainSampleDim, batch, cnum, wBits, logQ);  */
	timeutils.start("Encrypting trainlabel...");
	// encrypt the trainlabel
	for (long r = 0; r < rnum - 1; ++r) {
		for (long i = 0; i < cnum - 1; ++i) {

			complex<double>* pzLabel = new complex<double>[slots];
			for (long j = 0; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzLabel[batch * j + l].real(trainlabel[r*minbatchsize + j]);
					pzLabel[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainLabel[r*cnum+i], pzLabel, slots, wBits, logQ);
		}
			// i == cnum - 1       - the last cnum in each row
			complex<double>* pzLabel = new complex<double>[slots];
			for (long j = 0; j < minbatchsize; ++j) {
				long rest = factorDim - batch * (cnum - 1);

				for (long l = 0; l < rest; ++l) {
					pzLabel[batch * j + l].real(trainlabel[r*minbatchsize + j]);
					pzLabel[batch * j + l].imag(0);
				}
				for (long l = rest; l < batch; ++l) {
					pzLabel[batch * j + l].real(0);
					pzLabel[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainLabel[r*cnum+ cnum-1], pzLabel, slots, wBits, logQ);

    }
    // The last min-batch may consists of several ( trainSampleDim - minbatchsize * (rnum-1) ) rows of zeors.

        // r == rnum - 1       - the last rnum (the last min-batch)
    	long restrownum = trainSampleDim - minbatchsize * (rnum-1);
		for (long i = 0; i < cnum - 1; ++i) {			

			complex<double>* pzLabel = new complex<double>[slots];
			for (long j = 0; j < restrownum; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzLabel[batch * j + l].real(trainlabel[(rnum-1)*minbatchsize + j]);
					pzLabel[batch * j + l].imag(0);
				}
			}
			for (long j = restrownum; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzLabel[batch * j + l].real(0);
					pzLabel[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainLabel[(rnum - 1)*cnum+i], pzLabel, slots, wBits, logQ);
		}
			// i == cnum - 1       - the last cnum in each row
			complex<double>* pzLabel = new complex<double>[slots];
			for (long j = 0; j < restrownum; ++j) {
				long rest = factorDim - batch * (cnum - 1);
				for (long l = 0; l < rest; ++l) {
					pzLabel[batch * j + l].real(trainlabel[(rnum-1)*minbatchsize + j]);
					pzLabel[batch * j + l].imag(0);
				}
				for (long l = rest; l < batch; ++l) {
					pzLabel[batch * j + l].real(0);
					pzLabel[batch * j + l].imag(0);
				}
			}
			for (long j = restrownum; j < minbatchsize; ++j) {
				//long rest = factorDim - batch * (cnum - 1);
				for (long l = 0; l < batch; ++l) {
					pzLabel[batch * j + l].real(0);
					pzLabel[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainLabel[(rnum - 1)*cnum+ cnum-1], pzLabel, slots, wBits, logQ);
			delete[] pzLabel;
	timeutils.stop("trainlabel encryption");

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Encrypting trainlabel";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();


	timeutils.start("Encrypting traindata...");
	// encrypt the traindata
	for (long r = 0; r < rnum - 1; ++r) {
		for (long i = 0; i < cnum - 1; ++i) {

			complex<double>* pzData = new complex<double>[slots];
			for (long j = 0; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(traindata[r*minbatchsize + j][batch * i + l]);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainData[r*cnum+i], pzData, slots, wBits, logQ);
		}
			// i == cnum - 1       - the last cnum in each row
			complex<double>* pzData = new complex<double>[slots];
			for (long j = 0; j < minbatchsize; ++j) {
				long rest = factorDim - batch * (cnum - 1);

				for (long l = 0; l < rest; ++l) {
					pzData[batch * j + l].real(traindata[r*minbatchsize + j][batch * (cnum - 1) + l]);
					pzData[batch * j + l].imag(0);
				}
				for (long l = rest; l < batch; ++l) {
					pzData[batch * j + l].real(0);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainData[r*cnum+ cnum-1], pzData, slots, wBits, logQ);

    }
    // The last min-batch may consists of several ( trainSampleDim - minbatchsize * (rnum-1) ) rows of zeors.

        // r == rnum - 1       - the last rnum (the last min-batch)
    	restrownum = trainSampleDim - minbatchsize * (rnum-1);
		for (long i = 0; i < cnum - 1; ++i) {			

			complex<double>* pzData = new complex<double>[slots];
			for (long j = 0; j < restrownum; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(traindata[(rnum-1)*minbatchsize + j][batch * i + l]);
					pzData[batch * j + l].imag(0);
				}
			}
			for (long j = restrownum; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(0);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainData[(rnum - 1)*cnum+i], pzData, slots, wBits, logQ);
		}
			// i == cnum - 1       - the last cnum in each row
			complex<double>* pzData = new complex<double>[slots];
			for (long j = 0; j < restrownum; ++j) {
				long rest = factorDim - batch * (cnum - 1);
				for (long l = 0; l < rest; ++l) {
					pzData[batch * j + l].real(traindata[(rnum-1)*minbatchsize + j][batch * (cnum - 1) + l]);
					pzData[batch * j + l].imag(0);
				}
				for (long l = rest; l < batch; ++l) {
					pzData[batch * j + l].real(0);
					pzData[batch * j + l].imag(0);
				}
			}
			for (long j = restrownum; j < minbatchsize; ++j) {
				//long rest = factorDim - batch * (cnum - 1);
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(0);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encTrainData[(rnum - 1)*cnum+ cnum-1], pzData, slots, wBits, logQ);
			delete[] pzData;
	timeutils.stop("traindata encryption");

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Encrypting traindata";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();


	// zDataTrain and zDataTest are used for testing in the plaintext environment 
	// zData = (Y,Y@X)
	double** zDataTrain = new double*[trainSampleDim];
	for(int i=0;i<trainSampleDim;++i)
	{
		double* zi = new double[factorDim];
		zi[0] = trainlabel[i];
		for(int j=1;j<factorDim;++j)
			zi[j] = zi[0]*traindata[i][j];
		zDataTrain[i] = zi;
	}
	// zDataTest is only used for Cross-Validation test, not necessary for training LG model.
	// zData = (Y,Y@X)
	double** zDataTest = new double*[testSampleDim];
	for(int i=0;i<testSampleDim;++i)
	{
		double* zi = new double[factorDim];
		zi[0] = testlabel[i];
		for(int j=1;j<factorDim;++j)
			zi[j] = zi[0]*testdata[i][j];
		zDataTest[i] = zi;
	}

	// To compute the start value x0 for Newton's Method to invert the diagonal elements of fixed-Hessian matrix B.
	timeutils.start("Encrypting x0:=sumxij...");
	// To store the inv(B) in the quadratic gradient with the help of x0 
	// To Get the sum(Xij) to construct x0
	double* x0 = new double[rnum*cnum];
	for (long r = 0; r < rnum; ++r) {
		double sumxij = 0.0;
		for (long i = r*minbatchsize; 
			      i < ((r+1)*minbatchsize > trainSampleDim?trainSampleDim:(r+1)*minbatchsize); ++i) {
			for (long j = 0; j < factorDim; ++j)   sumxij += traindata[i][j];
		}
			
		sumxij = .25*sumxij; // the 1-st B[i][i], namely B[0][0]
		double x0r = 2.0 / sumxij  *    .6; // x0 < 2/a, set x0 := 1.8/a
		for (long c = 0; c < cnum; ++c)  
			x0[r*cnum + c] = x0r;

    }

	Ciphertext* encBinv = new Ciphertext[rnum*cnum]; // to store the final result inv(B)
	NTL_EXEC_RANGE(rnum*cnum, first, last);
	for (long i = first; i < last; ++i) {

		scheme.encryptSingle(encBinv[i], x0[i], wBits+wBits, logQ);
		encBinv[i].n = slots;
	}
	NTL_EXEC_RANGE_END;
	timeutils.stop("Encrypting x0:=sumxij...");         

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Encrypting x0:=sumxij";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();

	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */
	timeutils.start("Encrypting wData and vData...");
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		// scheme.encryptZeros(encWData[i], slots, wBits, encZData[0].logq); // To Make encVData[0].logq==encZData[0].logq
		scheme.encryptSingle(encWData[i], 0.0, wBits, logQ);
		encWData[i].n = slots;

		encVData[i].copy(encWData[i]);
	}
	NTL_EXEC_RANGE_END;
	timeutils.stop("wData and vData encryption");         

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Encrypting weight vector vData";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();
	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                                                    //
	//                        Client sent (encTrainData, encTrainLabel, enc(x0), encWData, and encVData) to Server                        //
    //                                                                                                                                    //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////// From now on, the server starts its work on what client sent to it. //////////////////////////////////
	
	/* --------------------- To Calculate the B from suming each row of  X.T*X (G) --------------------- */
	cout<<"--------------------- To Calculate the inv(B) in quadratic gradient --------------------- "<<endl<<endl<<endl;
	timeutils.start("Calculating the inverses of B[i][i]");
	/* Step 1. Sum Each Row To Its First Element */
	// encXTX[i] is to hold the sum of each row in every min-batch
 	Ciphertext* encXTX = new Ciphertext[rnum];
	for(long i=0; i<rnum; ++i)	encXTX[i].copy(encTrainData[i*cnum]);

	/* Sum the Batchs in the same row of the same min-batches To Get (rnum) Batches */
	for (long i = 0; i < rnum; ++i) {
		for (long j = 1; j < cnum; ++j) {
			scheme.addAndEqual(encXTX[i], encTrainData[i*cnum + j]);
		}
	}

	/* For Each Batch (ciphertext), Sum Itself Inside */
	NTL_EXEC_RANGE(rnum, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext rot;
		for (long l = 0; l < bBits; ++l) {
			scheme.leftRotateFast(rot, encXTX[i], (1 << l));
			scheme.addAndEqual(encXTX[i], rot);
		}
		rot.kill();
	}
	NTL_EXEC_RANGE_END


	/* Filtering out the first columns of each encXTX[i] and spanning it to the full rows */
	NTL_EXEC_RANGE(rnum, first, last);
	for (long i = first; i < last; ++i) {
		//scheme.multByPolyNTTAndEqual(encXTX[i], rpoly, pBits, pBits);
		scheme.multByPolyAndEqual(encXTX[i], dummy, pBits); //> logp: 2 * wBits + pBits
		Ciphertext tmp;
		for (long l = 0; l < bBits; ++l) {
			scheme.rightRotateFast(tmp, encXTX[i], (1 << l));
			scheme.addAndEqual(encXTX[i], tmp);
		}
		tmp.kill();
	}
	NTL_EXEC_RANGE_END
	/* THIS WILL INCREASE logp BY 'pBits', BUT WILL KEEP logq STILL */

	//Now, each row of encXTX[i*cnum] consists of the same value (sum(XTX[i][*])

	cout<<" --------------------- Print The Sum of Each Row --------------------- "<<endl;
	complex<double>* dcip = scheme.decrypt(secretKey, encXTX[0]);
	for (long j = 0; j < 7; ++j) {
		for (long l = 0; l < batch; ++l) {
			cout << setiosflags(ios::fixed) << setprecision(10) << dcip[batch * j + l].real() << "  ";
		}
		cout << endl;
	}
	cout<<" --------------------- Print The Sum of Each Row --------------------- "<<endl<<endl;

	/* Step 2. Each (i-th) Column (Inner*Product) the result of Step 1. */
	/*         This step has used/realized the special order of Bonte's method            */
	Ciphertext* encB = new Ciphertext[rnum*cnum];
	for (long i = 0; i < rnum*cnum; ++i) encB[i].copy(encTrainData[i]); 
	
	for (long i = 0; i < rnum; ++i) {

		NTL_EXEC_RANGE(cnum, first, last);
		for(long j = first; j < last; ++j){
			scheme.multAndEqual(encB[cnum*i + j], encXTX[i]);

			scheme.reScaleByAndEqual(encB[cnum*i + j], encXTX[i].logp); // will decrease the logq
			// put off this operation to keep its precision
		}
		NTL_EXEC_RANGE_END

	}

	/* To sum inside the  Batches from the each min-batch, so get every B[i][i] in each row */
	NTL_EXEC_RANGE(rnum*cnum, first, last);
	for (long j = first; j < last; ++j) {
		Ciphertext rot;
		long batch = 1 << bBits;
		for (long l = 0; l < minBatchDimBits; ++l) {
			scheme.leftRotateFast(rot, encB[j], (batch << l));
			scheme.addAndEqual(encB[j], rot);
		}
		rot.kill();
	}
	NTL_EXEC_RANGE_END
	/* NOW, EACH ROW OF encTrainData[j] for j=[0,cnum) HAS EVERY B[i][i] */
	/* Now, each column of encXXX[i] consists of the same value B[i][i]  */

	for(long i=0; i<rnum; ++i) encXTX[i].kill();
	delete[] encXTX;

	cout<<" --------------- Print The Sum of Each Column(B[i][i]) --------------- "<<endl;
	complex<double>* dctp = scheme.decrypt(secretKey, encB[0]);
	for (long j = 0; j < 7; ++j) {
		for (long l = 0; l < batch; ++l) {
			cout << setiosflags(ios::fixed) << setprecision(10) << dctp[batch * 0 + l].real() << "  ";
		}
		cout << endl;
	}
	cout<<" --------------- Print The Sum of Each Column(B[i][i]) --------------- "<<endl<<endl;


	/* Next : DONT forget to * .25 and add the epsion to make it positive */
	double epsion = 1e-8;   epsion *= .25;
	NTL_EXEC_RANGE(rnum*cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multByConstAndEqual(encB[i], .25, pBits);

		scheme.addConstAndEqual(encB[i], epsion, encB[i].logp);

		scheme.reScaleByAndEqual(encB[i], pBits);
	}
	NTL_EXEC_RANGE_END

	cout<<" ------------------- Print .25*B[i][i]+ .25*(1e-8) ------------------- "<<endl;
	complex<double>* dctq = scheme.decrypt(secretKey, encB[0]);
	for (long j = 0; j < 7; ++j) {
		for (long l = 0; l < batch; ++l) {
			cout << setiosflags(ios::fixed) << setprecision(10) << dctq[batch * j + l].real() << "  ";
		}
		cout << endl << setiosflags(ios::fixed) << setprecision(7);
	}
	cout<<" ------------------- Print .25*B[i][i]+ .25*(1e-8) ------------------- "<<endl<<endl;


	/* Step 3. Use Newton Method To Calculate the inv(B) */
	/*         x[k+1] = 2*x[k] - a*x[k]*x[k]             */
	/*                = x[k] * (2 - a*x[k])              */

	// To Get the sum(Xij) to construct x0

	cout<<" ---------- Use Newton Method To Calculate the inv(B) ---------- "<<endl;
	long NewtonIter = 9;
	// mod down the initial logq of encBinv[i] for Newton Iteration
	NTL_EXEC_RANGE(rnum*cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(encBinv[i], encB[i].logq);
	}
	NTL_EXEC_RANGE_END;
	
	for (long it = 0; it < NewtonIter; ++it)
	{
		// (i+1)-th iteration Newton Method For Zero Point (encBinv[i])
		Ciphertext* encTemp = new Ciphertext[rnum*cnum];

		NTL_EXEC_RANGE(rnum*cnum, first, last);
		for (long i = first; i < last; ++i) {
			encTemp[i].copy(encBinv[i]);

			// square... may avoid the error : x*x<0, do not use mult...
			scheme.squareAndEqual(encTemp[i]);                                     // encTemp = x1 * x1

			scheme.multAndEqual(encTemp[i], encB[i]);                              // encTemp = a * x1 * x1

			scheme.addAndEqual(encBinv[i], encBinv[i]);                            // encBinv[i] = 2 * x1
			                                                                       // avoid to use *, use + instead

            scheme.reScaleByAndEqual(encTemp[i],encTemp[i].logp-encBinv[i].logp);  // MAKE SURE : encBinv[i] and encTemp[i] share
			scheme.modDownToAndEqual(encBinv[i], encTemp[i].logq);                 // the same logp and logq, so they can add or *

			scheme.subAndEqual(encBinv[i], encTemp[i]);                            // encBinv = 2 * x1  - a * x1 * x1

			//scheme.reScaleByAndEqual(encBinv[i], 10);// make encBinv[i].logp equal to the encGrad[i].logp below, but encBinv@encGrad don't need that!
		}
		NTL_EXEC_RANGE_END

		for(long i=0; i<cnum; ++i) encTemp[i].kill();
		delete[] encTemp;

		cout << "after "<<(it+1)<<"-th iteration: encBinv[0].logq = " << encBinv[0].logq <<  ", encBinv[0].logp = " << encBinv[0].logp << endl;
		complex<double>* dcth = scheme.decrypt(secretKey, encBinv[0]);
		for (long j = 0; j < 7; ++j) {
			for (long l = 0; l < batch; ++l) {
				cout << setiosflags(ios::fixed) << setprecision(10)	<< dcth[batch * j + l].real() << "  ";
			}
			cout << endl;
		}
		cout << "after "<<(it+1)<<"-th iteration: encBinv[0].logq = " << encBinv[0].logq <<  ", encBinv[0].logp = " << encBinv[0].logp<<endl<<endl;

	}
	
	for(long i=0; i<rnum*cnum; ++i) encB[i].kill();
	delete[] encB;
	/* Each ([i][0~batch)-th) column of encBinv[i] consists of the same value (inv(B[i][i]) */
	timeutils.stop("Calculating the inverses of B[i][i]");

	openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
	openFileTIMELabel<<","<<"Calculating the inverses of B[i][i]";  openFileTIMELabel.flush();
	openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
	openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();

				

	double alpha0, alpha1, eta, gamma;
	double enccor, encauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long iter = 0; iter < numIter; ++iter) {
			
		vector<int> randr;
  		for (int ir = 0; ir < rnum; ++ir) randr.push_back(ir);
  		cout << "randr[]: "; for (int& x: randr) cout << ' ' << x; cout << '\n';
  		unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  		shuffle (randr.begin(), randr.end(), default_random_engine(seed));
  		cout << "randr[]: "; for (int& x: randr) cout << ' ' << x; cout << '\n';

		//cout << endl << "NesterovWithG : "+ to_string(iter+1)+" -th iteration" << endl;
		for (long r = 0; r < rnum; ++r) {
			timeutils.start("NesterovWithGminBatch : "+ to_string(iter+1)+" -th iteration");

			eta = (1 - alpha0) / alpha1;
			double gamma = 1.0 / minbatchsize / (1 + iter + r);
			if (trainSampleDim % minbatchsize != 0 && randr[r] == rnum-1)
				gamma = 1.0 / (trainSampleDim % minbatchsize) / (1 + iter + r);

			cout << endl << endl;
			cout << " ----------------- the " << (iter+1) << "-th epoch ----------------- " << endl;
			cout << " ----------------- the " << (r+1) << "-th iteration ----------------- " << endl;
			cout << "choose the " << randr[r] << "-th minibatch to update the weight vector: " << endl; 

		    /* cipherGD.encZData(encZData, zDataTrain, slots, factorDim, trainSampleDim, batch, cnum, wBits, logQ);  */
			Ciphertext* encZData = new Ciphertext[cnum];
			//timeutils.start("encZData[i] = encTrainLabel[i] @ encTrainData[i] for i in range(cnum)");
			// To Get the encZData
			NTL_EXEC_RANGE(cnum, first, last);
			for(long i = first; i < last; ++i){
				encZData[i].copy(encTrainLabel[cnum*randr[r] + i]);
				scheme.multAndEqual(encZData[i], encTrainData[cnum*randr[r] + i]);
				scheme.reScaleByAndEqual(encZData[i], encTrainData[cnum*randr[r] + i].logp);
			}
			NTL_EXEC_RANGE_END
			//timeutils.stop("encZData[i] for i in range(cnum) is done");
			/* cipherGD.encZData(encZData, zDataTrain, slots, factorDim, trainSampleDim, batch, cnum, wBits, logQ);  */

				
			/* cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits); */

				/* CipherGD::encInnerProduct(encIP, encZData, encWData, rpoly, cnum, bBits, wBits, pBits); */
					Ciphertext* encIPvec = new Ciphertext[cnum];

					/* For Each Batch, Sum Itself Inside */
					NTL_EXEC_RANGE(cnum, first, last);
					for (long i = first; i < last; ++i) {
						// MAKE SURE : encZData[i].logq >= encVData[i].logq (and of course : encZData[i].logp == encVData[i].logp)
						encIPvec[i].copy(encZData[i]);
						if (encIPvec[i].logq > encVData[i].logq){ //cout << "FF" << endl;
							scheme.modDownToAndEqual(encIPvec[i], encVData[i].logq);
						}
						if (encIPvec[i].logq < encVData[i].logq){ //cout << "EE" << endl;
							scheme.modDownToAndEqual(encVData[i], encIPvec[i].logq);
						}
						// V is the final weights to store the result weights.
						scheme.multAndEqual(encIPvec[i], encVData[i]);                // encIPvec = ENC(zData) .* ENC(V)

						/* For Each Batch (==ciphertext), Sum Itself Inside, Result in Each Row consisting of the same value */
						Ciphertext rot;                                               // encIPvec = ENC(zData) @  ENC(V)
						for (long l = 0; l < bBits; ++l) {
							scheme.leftRotateFast(rot, encIPvec[i], (1 << l));
							scheme.addAndEqual(encIPvec[i], rot);
						}
						rot.kill();
					}
					NTL_EXEC_RANGE_END


					/* Sum All Batchs To Get One Batch */
					Ciphertext encIP; encIP.copy(encIPvec[0]);             // to store the sum of all batches
					for (long i = 1; i < cnum; ++i) {
						scheme.addAndEqual(encIP, encIPvec[i]);
					}

					/* Sum This Batch Inside To Get The Inner Product */
					//scheme.multByPolyNTTAndEqual(encIP, rpoly, pBits, pBits);
					scheme.multByPolyAndEqual(encIP, dummy, pBits); //> logp: 2 * wBits + pBits
					Ciphertext tmp;
					for (long l = 0; l < bBits; ++l) {
						scheme.rightRotateFast(tmp, encIP, (1 << l));
						scheme.addAndEqual(encIP, tmp);
					}
					tmp.kill();
					/* THIS WILL INCREASE logp BY 'pBits', BUT WILL KEEP logq STILL */
					scheme.reScaleByAndEqual(encIP, pBits); // -5 is to make logp equal to encGrad[i].logp

					for(long i=0; i<cnum; ++i) encIPvec[i].kill();
					delete[] encIPvec;
				 //CipherGD::encInnerProduct(encIP, encZData, encVData, rpoly, cnum, bBits, wBits, pBits); 

				Ciphertext* encGrad = new Ciphertext[cnum];

				/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */

						Ciphertext encIP2(encIP);
						scheme.multAndEqual(encIP2, encIP);

						// IT IS VERY IMPORT TO KEEP THE logp BIG ENOUGH TO PRESENT THE NUMBER !!  MAKE SURE encIP2.logp>=35
						scheme.reScaleByAndEqual(encIP2, encIP.logp);                // For now, encIP.logp is big enough


						if( iter * rnum  + r < 12 ){
							//////////////////////////////////////// when iteration < 05 ////////////////////////////////////////
							cout << endl << "INSIDE iter < 12; poly7 = ";
							cout << setiosflags(ios::showpos) << degree7_12[0] << " ";
							cout << setiosflags(ios::showpos) << degree7_12[1] << "x " ;
							cout << setiosflags(ios::showpos) << degree7_12[2] << "x^3 ";
							cout << setiosflags(ios::showpos) << degree7_12[3] << "x^5 ";
							cout << setiosflags(ios::showpos) << degree7_12[4] << "x^7 " << endl << endl;
							cout << std::noshowpos;

							Ciphertext encIP4;
							scheme.square(encIP4, encIP2);
							scheme.reScaleByAndEqual(encIP4, encIP2.logp);

							Ciphertext encIP2c;
							scheme.multByConst(encIP2c, encIP2, degree7_12[3] / degree7_12[4], wBits);
							scheme.reScaleByAndEqual(encIP2c, wBits);

							if(encIP4.logp != encIP2c.logp) {cout<<"encIP4.logp!=encIP2c.logp"; exit(0); }
							if(encIP4.logq > encIP2c.logq) scheme.modDownToAndEqual(encIP4, encIP2c.logq);
							if(encIP4.logq < encIP2c.logq) scheme.modDownToAndEqual(encIP2c, encIP4.logq);
							scheme.addAndEqual(encIP4, encIP2c);

							//scheme.addConstAndEqual(encIP4, degree7[2] / degree7[4], wBits + 10);
							scheme.addConstAndEqual(encIP4, degree7_12[2] / degree7_12[4], encIP4.logp);

							NTL_EXEC_RANGE(cnum, first, last);
							for (long i = first; i < last; ++i) {
								Ciphertext tmp;
								scheme.multByConst(tmp, encZData[i], (1+gamma)  * degree7_12[1], wBits);

								scheme.modDownToAndEqual(tmp, encIP.logq);

								if(tmp.logq != encIP.logq) {cout << "$$#$$" << endl;exit(0);}

								scheme.multAndEqual(tmp, encIP);
								scheme.reScaleByAndEqual(tmp, encIP.logp);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_12[0], wBits);
								//scheme.reScaleByAndEqual(encGrad[i], pBits);
								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);

								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);

								scheme.addAndEqual(tmp, encGrad[i]);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_12[4], wBits + wBits);
								scheme.reScaleByAndEqual(encGrad[i], wBits);

								scheme.modDownToAndEqual(encGrad[i], encIP.logq);

								scheme.multAndEqual(encGrad[i], encIP);

								Ciphertext ctIP2(encIP2);
								if(encGrad[i].logq > ctIP2.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP2.logq);
								if(encGrad[i].logq < ctIP2.logq)
									scheme.modDownToAndEqual(ctIP2, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], ctIP2);
								scheme.reScaleByAndEqual(encGrad[i], ctIP2.logp);

								Ciphertext ctIP4(encIP4);
								if(encGrad[i].logq > ctIP4.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP4.logq);
								if(encGrad[i].logq < ctIP4.logq)
									scheme.modDownToAndEqual(ctIP4, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], encIP4);
								scheme.reScaleByAndEqual(encGrad[i], ctIP4.logp);

								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);
								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);
								scheme.addAndEqual(encGrad[i], tmp);

								tmp.kill();
								ctIP2.kill();
								ctIP4.kill();

							}
							NTL_EXEC_RANGE_END;

							encIP4.kill();
							encIP2c.kill();

						}else if( iter * rnum  + r < 24 ){
							//////////////////////////////////////// when iteration < 10 ////////////////////////////////////////
							cout << endl << "INSIDE iter < 24; poly7 = ";
							cout << setiosflags(ios::showpos) << degree7_24[0] << " ";
							cout << setiosflags(ios::showpos) << degree7_24[1] << "x " ;
							cout << setiosflags(ios::showpos) << degree7_24[2] << "x^3 ";
							cout << setiosflags(ios::showpos) << degree7_24[3] << "x^5 ";
							cout << setiosflags(ios::showpos) << degree7_24[4] << "x^7 " << endl << endl;
							cout << std::noshowpos;

							Ciphertext encIP4;
							scheme.square(encIP4, encIP2);
							scheme.reScaleByAndEqual(encIP4, encIP2.logp);

							Ciphertext encIP2c;
							scheme.multByConst(encIP2c, encIP2, degree7_24[3] / degree7_24[4], wBits);
							scheme.reScaleByAndEqual(encIP2c, wBits);

							if(encIP4.logp != encIP2c.logp) {cout<<"encIP4.logp!=encIP2c.logp"; exit(0); }
							if(encIP4.logq > encIP2c.logq) scheme.modDownToAndEqual(encIP4, encIP2c.logq);
							if(encIP4.logq < encIP2c.logq) scheme.modDownToAndEqual(encIP2c, encIP4.logq);
							scheme.addAndEqual(encIP4, encIP2c);

							//scheme.addConstAndEqual(encIP4, degree7[2] / degree7[4], wBits + 10);
							scheme.addConstAndEqual(encIP4, degree7_24[2] / degree7_24[4], encIP4.logp);

							NTL_EXEC_RANGE(cnum, first, last);
							for (long i = first; i < last; ++i) {
								Ciphertext tmp;
								scheme.multByConst(tmp, encZData[i], (1+gamma)  * degree7_24[1], wBits);

								scheme.modDownToAndEqual(tmp, encIP.logq);

								if(tmp.logq != encIP.logq) {cout << "$$#$$" << endl;exit(0);}

								scheme.multAndEqual(tmp, encIP);
								scheme.reScaleByAndEqual(tmp, encIP.logp);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_24[0], wBits);
								//scheme.reScaleByAndEqual(encGrad[i], pBits);
								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);

								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);

								scheme.addAndEqual(tmp, encGrad[i]);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_24[4], wBits + wBits);
								scheme.reScaleByAndEqual(encGrad[i], wBits);

								scheme.modDownToAndEqual(encGrad[i], encIP.logq);

								scheme.multAndEqual(encGrad[i], encIP);

								Ciphertext ctIP2(encIP2);
								if(encGrad[i].logq > ctIP2.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP2.logq);
								if(encGrad[i].logq < ctIP2.logq)
									scheme.modDownToAndEqual(ctIP2, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], ctIP2);
								scheme.reScaleByAndEqual(encGrad[i], ctIP2.logp);

								Ciphertext ctIP4(encIP4);
								if(encGrad[i].logq > ctIP4.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP4.logq);
								if(encGrad[i].logq < ctIP4.logq)
									scheme.modDownToAndEqual(ctIP4, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], encIP4);
								scheme.reScaleByAndEqual(encGrad[i], ctIP4.logp);

								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);
								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);
								scheme.addAndEqual(encGrad[i], tmp);

								tmp.kill();
								ctIP2.kill();
								ctIP4.kill();

							}
							NTL_EXEC_RANGE_END;

							encIP4.kill();
							encIP2c.kill();

						}else{
							//////////////////////////////////////// when iteration < 30 ////////////////////////////////////////
							cout << endl << "INSIDE iter < 36; poly7 = ";
							cout << setiosflags(ios::showpos) << degree7_36[0] << " ";
							cout << setiosflags(ios::showpos) << degree7_36[1] << "x " ;
							cout << setiosflags(ios::showpos) << degree7_36[2] << "x^3 ";
							cout << setiosflags(ios::showpos) << degree7_36[3] << "x^5 ";
							cout << setiosflags(ios::showpos) << degree7_36[4] << "x^7 " << endl << endl;
							cout << std::noshowpos;

							if(iter * rnum  + r > 30){
								cout << endl << "THE NUMBER OF MAX ITERATION SHOULD BE LESS THAN 30!" << endl;
								//exit(0);
							}
							if(iter * rnum  + r > 36){
								cout << endl << "The Number of Max Iteration should be less than 35!" << endl;
								cout << "otherwise, the poly7 should be replaced with a larger range poly!" << endl;
								exit(0);
							}

							Ciphertext encIP4;
							scheme.square(encIP4, encIP2);
							scheme.reScaleByAndEqual(encIP4, encIP2.logp);

							Ciphertext encIP2c;
							scheme.multByConst(encIP2c, encIP2, degree7_36[3] / degree7_36[4], wBits);
							scheme.reScaleByAndEqual(encIP2c, wBits);

							if(encIP4.logp != encIP2c.logp) {cout<<"encIP4.logp!=encIP2c.logp"; exit(0); }
							if(encIP4.logq > encIP2c.logq) scheme.modDownToAndEqual(encIP4, encIP2c.logq);
							if(encIP4.logq < encIP2c.logq) scheme.modDownToAndEqual(encIP2c, encIP4.logq);
							scheme.addAndEqual(encIP4, encIP2c);

							//scheme.addConstAndEqual(encIP4, degree7[2] / degree7[4], wBits + 10);
							scheme.addConstAndEqual(encIP4, degree7_36[2] / degree7_36[4], encIP4.logp);

							NTL_EXEC_RANGE(cnum, first, last);
							for (long i = first; i < last; ++i) {
								Ciphertext tmp;
								scheme.multByConst(tmp, encZData[i], (1+gamma)  * degree7_36[1], wBits);

								scheme.modDownToAndEqual(tmp, encIP.logq);

								if(tmp.logq != encIP.logq) {cout << "$$#$$" << endl;exit(0);}

								scheme.multAndEqual(tmp, encIP);
								scheme.reScaleByAndEqual(tmp, encIP.logp);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_36[0], wBits);
								//scheme.reScaleByAndEqual(encGrad[i], pBits);
								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);

								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);

								scheme.addAndEqual(tmp, encGrad[i]);

								//////////////////////////////////////////////////////////////////////////////
								scheme.multByConst(encGrad[i], encZData[i], (1+gamma)  * degree7_36[4], wBits + wBits);
								scheme.reScaleByAndEqual(encGrad[i], wBits);

								scheme.modDownToAndEqual(encGrad[i], encIP.logq);

								scheme.multAndEqual(encGrad[i], encIP);

								Ciphertext ctIP2(encIP2);
								if(encGrad[i].logq > ctIP2.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP2.logq);
								if(encGrad[i].logq < ctIP2.logq)
									scheme.modDownToAndEqual(ctIP2, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], ctIP2);
								scheme.reScaleByAndEqual(encGrad[i], ctIP2.logp);

								Ciphertext ctIP4(encIP4);
								if(encGrad[i].logq > ctIP4.logq)
									scheme.modDownToAndEqual(encGrad[i], ctIP4.logq);
								if(encGrad[i].logq < ctIP4.logq)
									scheme.modDownToAndEqual(ctIP4, encGrad[i].logq);
								scheme.multAndEqual(encGrad[i], encIP4);
								scheme.reScaleByAndEqual(encGrad[i], ctIP4.logp);

								if(tmp.logp > encGrad[i].logp) scheme.reScaleByAndEqual(tmp,tmp.logp-encGrad[i].logp);
								if(tmp.logp < encGrad[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-tmp.logp);
								if(tmp.logq > encGrad[i].logq) scheme.modDownToAndEqual(tmp, encGrad[i].logq);
								if(tmp.logq < encGrad[i].logq) scheme.modDownToAndEqual(encGrad[i], tmp.logq);
								scheme.addAndEqual(encGrad[i], tmp);

								tmp.kill();
								ctIP2.kill();
								ctIP4.kill();

							}
							NTL_EXEC_RANGE_END;

							encIP4.kill();
							encIP2c.kill();

						}
						encIP2.kill();
						encIP.kill();

					// Sum Each Column of encGrad[i] To Get the Final gradient : (1 - sigm(yWTx)) * Y.T @ X
					NTL_EXEC_RANGE(cnum, first, last);
					for (long i = first; i < last; ++i) {
						Ciphertext tmp;
						for (long l = bBits; l < sBits; ++l) {
							scheme.leftRotateFast(tmp, encGrad[i], (1 << l));
							scheme.addAndEqual(encGrad[i], tmp);
						}

						Ciphertext ctBinv(encBinv[randr[r]*cnum+i]);
						if (encGrad[i].logq > ctBinv.logq)
							scheme.modDownToAndEqual(encGrad[i], ctBinv.logq);
						if (encGrad[i].logq < ctBinv.logq)
							scheme.modDownToAndEqual(ctBinv, encGrad[i].logq);

						scheme.multAndEqual(encGrad[i], encBinv[randr[r]*cnum+i]);
						scheme.reScaleByAndEqual(encGrad[i], encBinv[randr[r]*cnum+i].logp);

						tmp.kill();
						ctBinv.kill();
					}
					NTL_EXEC_RANGE_END;
					/* Each ([i][0~batch)-th) column of encGrad[i] consists of the same value (gamma * encGrad[i][0~batch)) */
				/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */


				/* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */
					NTL_EXEC_RANGE(cnum, first, last);
					for (long i = first; i < last; ++i) {

						if(encGrad[i].logp > encVData[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-encVData[i].logp);
						if(encGrad[i].logp < encVData[i].logp) scheme.reScaleByAndEqual(encVData[i], encVData[i].logp-encGrad[i].logp);
						scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);

						Ciphertext ctmpw;
						scheme.add(ctmpw, encVData[i], encGrad[i]); 					// encGrad[i] has already self-multiplied with gamma
						                                                                // ctmpw = encVData[i] - encGrad[i]

						scheme.multByConst(encVData[i], ctmpw, 1. - eta, pBits);        // encVData[i] = ( 1. - eta ) * ctmpw
						//scheme.reScaleByAndEqual(encVData[i], pBits-5);

						scheme.multByConstAndEqual(encWData[i], eta, pBits);            // encWData[i] = eta * encWData[i]
						//scheme.reScaleByAndEqual(encWData[i], pBits-5);

						if (encWData[i].logq > encVData[i].logq) scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
						if (encWData[i].logq < encVData[i].logq) scheme.modDownToAndEqual(encVData[i], encWData[i].logq);
						if (encWData[i].logp != encVData[i].logp) { cout << "logp != logp" ;exit(0); }

						scheme.addAndEqual(encVData[i], encWData[i]);                   // encVData[i] = encVData[i] + encWData[i]
						                                                 // encVData[i] = ( 1. - eta ) * ctmpw + eta * encWData[i]

						scheme.reScaleByAndEqual(encVData[i], pBits);
						encWData[i].copy(ctmpw);

						ctmpw.kill();
					}
					NTL_EXEC_RANGE_END;
		        /* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */

				for(long i=0;i<cnum;++i) encGrad[i].kill();
				delete[] encGrad;
			/* cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits); */


			timeutils.stop("NesterovWithG : "+ to_string(iter+1)+" -th iteration");

			openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
			openFileTIMELabel<<","<<"NesterovWithGminBatch : "+ to_string(iter+1) + " -th epoch." + to_string(r+1) + "-th iteration";  openFileTIMELabel.flush();
			openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
			openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();


			for(long i=0; i<cnum; ++i) encZData[i].kill();
			delete[] encZData;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			cout << endl << "---------- TEST : THE " << iter*rnum + r +1 << "-th ITERATION : Weights, AUC, MLE ----------" << endl;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);     */
				for (long i = 0; i < (cnum - 1); ++i) {
					complex<double>* dcvv = scheme.decrypt(secretKey, encVData[i]);
					for (long j = 0; j < batch; ++j) {
						cvData[batch * i + j] = dcvv[j].real();
					}
					delete[] dcvv;
				}
				complex<double>* dcvv = scheme.decrypt(secretKey, encVData[cnum-1]);
				long rest = factorDim - batch * (cnum - 1);
				for (long j = 0; j < rest; ++j) {
					cvData[batch * (cnum - 1) + j] = dcvv[j].real();
				}
				delete[] dcvv;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits); */
			cout << "Current cWdata (encVData) : " << endl;
			for(long i = 0; i < factorDim; ++i) cout << setiosflags(ios::fixed) << setprecision(12) << cvData[i] << ",\t"; cout << endl;
			
			openFileTestAUC<<","<<MyTools::calculateAUC(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc);    openFileTestAUC.flush();
			openFileTrainAUC<<","<<MyTools::calculateAUC(zDataTrain, cvData, factorDim, trainSampleDim, enccor, encauc); openFileTrainAUC.flush();

			openFileTestACC<<","<<MyTools::calculateACC(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc);    openFileTestACC.flush();
			openFileTrainACC<<","<<MyTools::calculateACC(zDataTrain, cvData, factorDim, trainSampleDim, enccor, encauc); openFileTrainACC.flush();

			openFileTestMLE<<","<<MyTools::calculateMLE(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc);    openFileTestMLE.flush();
			openFileTrainMLE<<","<<MyTools::calculateMLE(zDataTrain, cvData, factorDim, trainSampleDim, enccor, encauc); openFileTrainMLE.flush();
			
			cout << "Test MLE : " << MyTools::calculateMLE(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc) << endl;
			cout << "Test AUC : " << MyTools::calculateAUC(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc) << endl;
			cout << "Test ACC : " << MyTools::calculateACC(zDataTest, cvData, factorDim, testSampleDim, enccor, encauc) << endl;
			cout << " ----------------- end of the " << (iter+1) << "-th epoch ----------------- " << endl;
			cout << " ----------------- end of the " << (r+1) << "-th iteration ----------------- " << endl << endl; 
			cout << "--------------------------------------------------------------------------------" << endl;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			if ( encVData[0].logq < 450 ) {
			//if ( encVData[0].logq <= 450 && iter < numIter-1 || encVData[0].logq < wBits && iter == numIter-1) {

				timeutils.start("Use Bootstrap To Recrypt Ciphertext");

				NTL_EXEC_RANGE(cnum, first, last);
				for(long i = first; i < last; ++i){
				 	scheme.modDownToAndEqual(encWData[i], bootlogq);
				 	encWData[i].n = batch;
				 	scheme.bootstrapAndEqual(encWData[i], bootlogq, logQ, logT, logI);
				 	encWData[i].n = slots;
				 }
				NTL_EXEC_RANGE_END
				NTL_EXEC_RANGE(cnum, first, last);
				for(long i = first; i < last; ++i){
				 	scheme.modDownToAndEqual(encVData[i], bootlogq);
				 	encVData[i].n = batch;
				 	scheme.bootstrapAndEqual(encVData[i], bootlogq, logQ, logT, logI);
				 	encVData[i].n = slots;
				 }
				NTL_EXEC_RANGE_END

				timeutils.stop("Use Bootstrap To Recrypt Ciphertext");

				openFileTIME<<","<<timeutils.timeElapsed;  openFileTIME.flush();
				openFileTIMELabel<<","<<"Bootstrapping";  openFileTIMELabel.flush();
				openFileCurrMEM<<","<< ( MyTools::getCurrentRSS() >> 20 );  openFileCurrMEM.flush();
				openFilePeakMEM<<","<< ( MyTools::getPeakRSS() >> 20 );  openFilePeakMEM.flush();
			
			}
			/////////////////////////////////////////////////////////////////////////////
			//        BOOTSTRAPPING                                                    //
			//             Over and Out                                                //
			/////////////////////////////////////////////////////////////////////////////
			cout << "--------------------------------------------------------------------------------" << endl << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << endl << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl << endl;
		}
	}

	openFileTIME<<endl;      openFileTIME.flush();
	openFileTIMELabel<<endl; openFileTIMELabel.flush();
	openFileTestAUC<<endl ; openFileTestAUC.flush();
	openFileTrainAUC<<endl ; openFileTrainAUC.flush();
	openFileTestACC<<endl ; openFileTestACC.flush();
	openFileTrainACC<<endl ; openFileTrainACC.flush();
	openFileTestMLE<<endl ; openFileTestMLE.flush();
	openFileTrainMLE<<endl ; openFileTrainMLE.flush();
	openFileCurrMEM<<endl;  openFileCurrMEM.flush();
	openFilePeakMEM<<endl;  openFilePeakMEM.flush();


	openFileTIME.close();
	openFileTIMELabel.close();
	openFileTestAUC.close();
	openFileTrainAUC.close();
	openFileTestACC.close();
	openFileTrainACC.close();
	openFileTestMLE.close();
	openFileTrainMLE.close();
	openFileCurrMEM.close();
	openFilePeakMEM.close();
}

