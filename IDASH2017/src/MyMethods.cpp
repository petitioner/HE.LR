#include "MyMethods.h"
#include "SerializationUtils.h"  // error: ‘SerializationUtils’ has not been declared
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

double* MyMethods::testCryptoFullBatchNAGwithG(double **traindata,
        double *trainlabel, long factorDim, long trainSampleDim, long Epoch_Number,
        double **testdata, double *testlabel, long testSampleDim,
        string resultpath) {

	long wBits = 30;
	long pBits = 20;

	//long logN = MyTools::suggestLogN(80, logQ);
	long logN = 16;
	long logQ = 275; // 991.300840336 > logQ  to ensure 128-bit security level. Security Parameter λ
	long slots = 1 << (logN - 1);
	long sBits = (long) ceil(log2(slots));
	long fdimBits = (long) ceil(log2(factorDim));
	long sdimBits = (long) ceil(log2(trainSampleDim));

	long batch = 1 << 5; // Basically, batch is the Number of several factor dimensions.
	long bBits = (long) ceil(log2(batch)); // 2^batchBits = min( 2^logN / 2^sdimBits / 2, 2^fdimBits ) = min( N/2 /n, factorDim ) ;

	// the size of batch should be small than fatctorDim
	//if (batch >= factorDim) {
	//	cout << "batch >= factorDim!" << endl;
	//	exit(0);
	//}

	//long factorNum = 1 << (long)ceil(log2(factorDim));
	long cnum = (long) ceil(double(factorDim) / batch); // To Divide the whole Train Data into Several Batches (cnum Ciphertexts).
	if (cnum > (1 << fdimBits)) {
		cout << "cnum should be no more than factorDim!" << endl;
		exit(0);
	}

// We assume that the dataset is too large to be encrypted into one single ciphertext.
	long minbatchsize = slots / batch; // (each min-batch should be small to save the useless space of the last batch)
	long minBatchDimBits = (long) ceil(log2(minbatchsize));

	long rnum = (long) ceil((double) trainSampleDim / minbatchsize); // rnum : the number of min-batch, namely number of row

	cout << "logQ = " << logQ << ", logN = " << logN << endl;
	cout << "slots = " << slots << ", batch = " << batch << endl;
	cout << "rnum = " << rnum << ", cnum = " << cnum << endl;
	cout << "min-batch size = " << minbatchsize << ", minBatchDimBits = "
	     << minBatchDimBits << endl;

	string path = resultpath;
	ofstream openFileTrainAUC(path + "TrainAUC.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainACC(path + "TrainACC.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainMLE(path + "TrainMLE.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTestAUC(path + "TestAUC.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTestACC(path + "TestACC.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTestMLE(path + "TestMLE.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTIME(path + "TIME.csv",
	                      std::ofstream::out | std::ofstream::app);
	ofstream openFileTIMELabel(path + "TIMELabel.csv",
	                           std::ofstream::out | std::ofstream::app);
	ofstream openFileCurrMEM(path + "CurrMEM.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFilePeakMEM(path + "PeakMEM.csv",
	                         std::ofstream::out | std::ofstream::app);

	if (!openFileTrainAUC.is_open())
		cout << "Error: cannot read Train AUC file" << endl;
	if (!openFileTrainACC.is_open())
		cout << "Error: cannot read Train ACC file" << endl;
	if (!openFileTrainMLE.is_open())
		cout << "Error: cannot read Train MLE file" << endl;
	if (!openFileTestAUC.is_open())
		cout << "Error: cannot read Test AUC file" << endl;
	if (!openFileTestACC.is_open())
		cout << "Error: cannot read Test ACC file" << endl;
	if (!openFileTestMLE.is_open())
		cout << "Error: cannot read Test MLE file" << endl;
	if (!openFileTIME.is_open())
		cout << "Error: cannot read TIME file" << endl;
	if (!openFileTIMELabel.is_open())
		cout << "Error: cannot read TIME Label file" << endl;
	if (!openFileTIMELabel.is_open())
		cout << "Error: cannot read TIME Label file" << endl;
	if (!openFileCurrMEM.is_open())
		cout << "Error: cannot read Current MEMORY file" << endl;
	if (!openFilePeakMEM.is_open())
		cout << "Error: cannot read Peak MEMORY file" << endl;

	openFileTrainAUC << "TrainAUC";
	openFileTrainAUC.flush();
	openFileTrainACC << "TrainACC";
	openFileTrainACC.flush();
	openFileTrainMLE << "TrainMLE";
	openFileTrainMLE.flush();
	openFileTestAUC << "TestAUC";
	openFileTestAUC.flush();
	openFileTestACC << "TestACC";
	openFileTestACC.flush();
	openFileTestMLE << "TestMLE";
	openFileTestMLE.flush();
	openFileTIME << "TIME";
	openFileTIME.flush();
	openFileTIMELabel << "TIMELabel";
	openFileTIMELabel.flush();
	openFileCurrMEM << "MEMORY(GB)";
	openFileCurrMEM.flush();
	openFilePeakMEM << "MEMORY(GB)";
	openFilePeakMEM.flush();

	TimeUtils timeutils;
	timeutils.start("Scheme generating");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Scheme generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	cout << "NOW THE PEAK RSS IS: " << (MyTools::getPeakRSS() >> 20) << endl;

	timeutils.start("Bootstrap Key generating");
	long logT = 3;
	long logI = 4;
	long bootlogq = 30 + 10;
	long lognslots = (long) ceil(log2(batch));  //batch = factorDim / cnum;
	scheme.addBootKey(secretKey, lognslots, bootlogq + 4);
	timeutils.stop("Bootstrap Key generated");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Bootstrap Key generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	cout << "NOW THE CURR RSS IS: " << (MyTools::getCurrentRSS() >> 20) << endl;
	cout << "NOW THE PEAK RSS IS: " << (MyTools::getPeakRSS() >> 20) << endl;

	// Basically, rpoly is used to calculate the sum{row}
	timeutils.start("Polynomial generating...");
	ZZ *dummy = new ZZ[N];
	;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	complex<double> *pvals = new complex<double> [slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j] = 1.0;
	}
	scheme.ring.encode(dummy, pvals, slots, pBits);
	delete[] pvals;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	timeutils.stop("Polynomial generation");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Polynomial generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	double *cwData = new double[factorDim]();
	double *cvData = new double[factorDim]();

	Ciphertext *encTrainData = new Ciphertext[rnum * cnum];
	Ciphertext *encTrainLabel = new Ciphertext[rnum * cnum];

	Ciphertext *encXyZdata = new Ciphertext[rnum * cnum]; // TrainLabel @ TrainData = yX
	Ciphertext *encBinv = new Ciphertext[cnum];

	Ciphertext *encWData = new Ciphertext[cnum];
	Ciphertext *encVData = new Ciphertext[cnum];

	/* - - - - - - - - - - - - - - - - - - - - - - - - Client and Server - - - - - - - - - - - - - - - - - - - - - - - - */

	timeutils.start("Encrypting encXyZdata...");
	// encrypt the traindata
	for (long r = 0; r < rnum - 1; ++r) {
		for (long i = 0; i < cnum - 1; ++i) {

			complex<double> *pzData = new complex<double> [slots]();
			for (long j = 0; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(
					    trainlabel[r * minbatchsize + j]
					    * traindata[r * minbatchsize + j][batch * i
					                                      + l]);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encXyZdata[r * cnum + i], pzData, slots, wBits,
			               logQ);
			SerializationUtils::writeCiphertext(encXyZdata[r * cnum + i], "encXyZdata[" + std::to_string(r * cnum + i) + "].txt");
		}
		// i == cnum - 1       - the last cnum in each row
		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < minbatchsize; ++j) {
			long rest = factorDim - batch * (cnum - 1);

			for (long l = 0; l < rest; ++l) {
				pzData[batch * j + l].real(
				    trainlabel[r * minbatchsize + j]
				    * traindata[r * minbatchsize + j][batch
				                                      * (cnum - 1) + l]);
				pzData[batch * j + l].imag(0);
			}
			for (long l = rest; l < batch; ++l) {
				pzData[batch * j + l].real(0);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encXyZdata[r * cnum + cnum - 1], pzData, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encXyZdata[r * cnum + cnum - 1], "encXyZdata[" + std::to_string(r * cnum + cnum - 1) + "].txt");

	}
	// The last min-batch may consists of several ( trainSampleDim - minbatchsize * (rnum-1) ) rows of zeors.

	// r == rnum - 1       - the last rnum (the last min-batch)
	auto restrownum = trainSampleDim - minbatchsize * (rnum - 1);
	for (long i = 0; i < cnum - 1; ++i) {

		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < restrownum; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(
				    trainlabel[(rnum - 1) * minbatchsize + j]
				    * traindata[(rnum - 1) * minbatchsize + j][batch
				            * i + l]);
				pzData[batch * j + l].imag(0);
			}
		}
		for (long j = restrownum; j < minbatchsize; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(0);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encXyZdata[(rnum - 1) * cnum + i], pzData, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encXyZdata[(rnum - 1) * cnum + i], "encXyZdata[" + std::to_string((rnum - 1) * cnum + i) + "].txt");
	}
	// i == cnum - 1       - the last cnum in each row
	complex<double> *pzDatra = new complex<double> [slots]();
	for (long j = 0; j < restrownum; ++j) {
		long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < rest; ++l) {
			pzDatra[batch * j + l].real(
			    trainlabel[(rnum - 1) * minbatchsize + j]
			    * traindata[(rnum - 1) * minbatchsize + j][batch
			            * (cnum - 1) + l]);
			pzDatra[batch * j + l].imag(0);
		}
		for (long l = rest; l < batch; ++l) {
			pzDatra[batch * j + l].real(0);
			pzDatra[batch * j + l].imag(0);
		}
	}
	for (long j = restrownum; j < minbatchsize; ++j) {
		//long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < batch; ++l) {
			pzDatra[batch * j + l].real(0);
			pzDatra[batch * j + l].imag(0);
		}
	}
	scheme.encrypt(encXyZdata[(rnum - 1) * cnum + cnum - 1], pzDatra, slots,
	               wBits, logQ);
	SerializationUtils::writeCiphertext(encXyZdata[(rnum - 1) * cnum + cnum - 1], "encXyZdata[" + std::to_string((rnum - 1) * cnum + cnum - 1) + "].txt");
	delete[] pzDatra;
	timeutils.stop("encXyZdata encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting encXyZdata";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	double **Binv = MyTools::zInvBFromFile(traindata, factorDim,
	                                       trainSampleDim);
	for (long r = 0; r < trainSampleDim; ++r)
		for (long c = 0; c < factorDim; ++c)
			Binv[r][c] = min(Binv[r][c], 8.); // overflow:: beyond the precision!!
	timeutils.start("Encrypting Binv...");
	// encrypt the traindata

	for (long i = 0; i < cnum - 1; ++i) {

		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < minbatchsize; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(
				    Binv[j][batch * i + l]);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encBinv[i], pzData, slots, wBits, logQ);
		SerializationUtils::writeCiphertext(encBinv[i], "encBinv[" + std::to_string(i) + "].txt");
	}
	// i == cnum - 1       - the last cnum in each row
	complex<double> *pzData3 = new complex<double> [slots]();
	for (long j = 0; j < minbatchsize; ++j) {
		long rest = factorDim - batch * (cnum - 1);

		for (long l = 0; l < rest; ++l) {
			pzData3[batch * j + l].real(
			    Binv[j][batch * (cnum - 1) + l]);
			pzData3[batch * j + l].imag(0);
		}
		for (long l = rest; l < batch; ++l) {
			pzData3[batch * j + l].real(0);
			pzData3[batch * j + l].imag(0);
		}
	}
	scheme.encrypt(encBinv[cnum - 1], pzData3, slots, wBits,
	               logQ);
	SerializationUtils::writeCiphertext(encBinv[cnum - 1], "encBinv[" + std::to_string(cnum - 1) + "].txt");

	timeutils.stop("Binv encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting Binv";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	// zDataTrain and zDataTest are used for testing in the plaintext environment
	// zData = (Y,Y@X)
	double **zDataTrain = new double*[trainSampleDim];
	for (int i = 0; i < trainSampleDim; ++i) {
		double *zi = new double[factorDim]();
		zi[0] = trainlabel[i];
		for (int j = 1; j < factorDim; ++j)
			zi[j] = zi[0] * traindata[i][j];
		zDataTrain[i] = zi;
	}
	// zDataTest is only used for Cross-Validation test, not necessary for training LG model.
	// zData = (Y,Y@X)
	double **zDataTest = new double*[testSampleDim];
	for (int i = 0; i < testSampleDim; ++i) {
		double *zi = new double[factorDim]();
		zi[0] = testlabel[i];
		for (int j = 1; j < factorDim; ++j)
			zi[j] = zi[0] * testdata[i][j];
		zDataTest[i] = zi;
	}

	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */
	timeutils.start("Encrypting wData and vData...");
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		// scheme.encryptZeros(encWData[i], slots, wBits, encZData[0].logq); // To Make encVData[0].logq==encZData[0].logq
		scheme.encryptSingle(encWData[i], 0.0, wBits, logQ);
		encWData[i].n = slots;

		encVData[i].copy(encWData[i]);
	}
	NTL_EXEC_RANGE_END
	;
	timeutils.stop("wData and vData encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting weight vector vData";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();
	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                                                    //
	//                        Client sent (encTrainData, encTrainLabel, enc(x0), encWData, and encVData) to Server                        //
	//                                                                                                                                    //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// From now on, the server starts its work on what client sent to it. //////////////////////////////////
	cout << endl << endl << endl << endl << "scheme.leftRotKeyMap.size() = " << scheme.leftRotKeyMap.size() << endl << endl << endl << endl;
	cout << endl << endl << endl << endl;
	double min_lr = 1.0;
	double max_lr = 2.0;
	double total_steps = Epoch_Number * rnum;
	double exp_gamma = 2.5;
	cout << "min_lr = " << min_lr << endl;
	cout << "max_lr = " << max_lr << endl;
	cout << "total_steps = " << total_steps << endl;
	cout << "exp_gamma = " << exp_gamma << endl;
	cout << "len(Xminbatches) = rnum = " << rnum << endl;
	cout << endl << endl << endl << endl;

	double totalUpdatingWeightsTime = 0.0;
	double totalUpdatingWeights = 0.0;
	double totalRefreshingWeightsTime = 0.0;
	double totalBootstrappingOperations = 0.0;




	double alpha0, alpha1, eta, gamma;
	double enccor, encauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	cout << endl << "unsigned seed = chrono::system_clock::now().time_since_epoch().count();";
	cout << endl << "unsigned seed = " << (unsigned)seed << endl;

	for (long iter = 0; iter < Epoch_Number; ++iter) {

		vector<int> randr;
		for (int ir = 0; ir < rnum; ++ir)
			randr.push_back(ir);
		cout << "randr[]: ";
		for (int &x : randr)
			cout << ' ' << x;
		cout << '\n';

		shuffle(randr.begin(), randr.end(), default_random_engine(seed));
		cout << "randr[]: ";
		for (int &x : randr)
			cout << ' ' << x;
		cout << '\n';

		//cout << endl << "NesterovWithG : "+ to_string(iter+1)+" -th iteration" << endl;
		for (long r = 0; r < rnum; ++r) {
			timeutils.start("NesterovWithGminBatch : " + to_string(iter + 1) + " -th iteration");
			auto start_timeUpdatingWeights = std::chrono::high_resolution_clock::now();

			cout << endl << endl << endl << endl << "encVData[0].logq = "
			     << encVData[0].logq << endl << endl << endl << endl;


			eta = (1 - alpha0) / alpha1;

			auto iterations = iter * rnum + r;
			// self.max_lr - (self.max_lr - self.min_lr) * (self.current_step / self.total_steps)**self.gamma
			auto gamma = max_lr - (max_lr - min_lr) * pow(iterations / total_steps, exp_gamma );


			Ciphertext *encZData = new Ciphertext[rnum * cnum];
			NTL_EXEC_RANGE(rnum * cnum, first, last);
			for (long i = first; i < last; ++i) {
				//encZData[i].copy(encXyZdata[randr[r] * cnum + i]);
				encZData[i].copy(encXyZdata[i]);
			}
			NTL_EXEC_RANGE_END
			Ciphertext *encIPvec = new Ciphertext[rnum * cnum];

			/* For Each Batch, Sum Itself Inside */
			NTL_EXEC_RANGE(rnum * cnum, first, last);
			for (long i = first; i < last; ++i) {
				// MAKE SURE : encZData[i].logq >= encVData[i].logq (and of course : encZData[i].logp == encVData[i].logp)
				encIPvec[i].copy(encZData[i]);
				if (encIPvec[i].logq > encVData[i % cnum].logq) { //cout << "FF" << endl;
					scheme.modDownToAndEqual(encIPvec[i], encVData[i % cnum].logq);
				}
				if (encIPvec[i].logq < encVData[i % cnum].logq) { //cout << "EE" << endl;
					scheme.modDownToAndEqual(encVData[i % cnum], encIPvec[i].logq);
				}
				// V is the final weights to store the result weights.
				scheme.multAndEqual(encIPvec[i], encVData[i % cnum]); // encIPvec = ENC(zData) .* ENC(V)

			}
			NTL_EXEC_RANGE_END

			/* Sum All Batchs To Get One Batch */
			Ciphertext *encIP = new Ciphertext[rnum];

			NTL_EXEC_RANGE(rnum, first, last);
			for (long i = first; i < last; ++i) {
				encIP[i].copy(encIPvec[i * cnum]); // to store the sum of all batches

				for (long c = 1; c < cnum; ++c) {
					scheme.addAndEqual(encIP[i], encIPvec[i * cnum + c]);
				}

				/* For Each Batch (==ciphertext), Sum Itself Inside, Result in Each Row consisting of the same value */
				Ciphertext rot;               // encIPvec = ENC(zData) @  ENC(V)
				for (long l = 0; l < bBits; ++l) {
					scheme.leftRotateFast(rot, encIP[i], (1 << l));
					scheme.addAndEqual(encIP[i], rot);
				}
				rot.kill();

				/* Sum This Batch Inside To Get The Inner Product */
				//scheme.multByPolyNTTAndEqual(encIP, rpoly, pBits, pBits);
				scheme.multByPolyAndEqual(encIP[i], dummy, pBits);//> logp: 2 * wBits + pBits
				Ciphertext tmp;
				for (long l = 0; l < bBits; ++l) {
					scheme.rightRotateFast(tmp, encIP[i], (1 << l));
					scheme.addAndEqual(encIP[i], tmp);
				}
				tmp.kill();
				/* THIS WILL INCREASE logp BY 'pBits', BUT WILL KEEP logq STILL */
				scheme.reScaleByAndEqual(encIP[i], pBits); // -5 is to make logp equal to encGrad[i].logp

			}
			NTL_EXEC_RANGE_END

			for (long i = 0; i < rnum * cnum; ++i)
				encIPvec[i].kill();
			delete[] encIPvec;
			//CipherGD::encInnerProduct(encIP, encZData, encVData, rpoly, cnum, bBits, wBits, pBits);


			Ciphertext *encGrad = new Ciphertext[cnum];
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {
				scheme.encryptSingle(encGrad[i], 0.0, wBits, logQ);
				encGrad[i].n = slots;
			}
			NTL_EXEC_RANGE_END


			//////////////////////////////////////// when iteration < 05 ////////////////////////////////////////
			cout << endl << "INSIDE iter < 5;  poly3 = ";
			cout << setiosflags(ios::showpos) << degree3[0] << " ";
			cout << setiosflags(ios::showpos) << degree3[1] << "x ";
			cout << setiosflags(ios::showpos) << degree3[2] << "x^3 "
			     << endl << endl;
			cout << std::noshowpos;
			cout << "gamma = " << gamma << endl << endl << endl;


			for (long r = 0; r < rnum; ++r)
				//NTL_EXEC_RANGE(rnum, first, last);
				//for (long r = first; r < last; ++r)
			{

				/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */

				Ciphertext encIP2(encIP[r]);
				scheme.multAndEqual(encIP2, encIP[r]);

				// IT IS VERY IMPORT TO KEEP THE logp BIG ENOUGH TO PRESENT THE NUMBER !!  MAKE SURE encIP2.logp>=35
				scheme.reScaleByAndEqual(encIP2, encIP[r].logp); // For now, encIP.logp is big enough

				scheme.addConstAndEqual(encIP2, degree3[1] / degree3[2],
				                        encIP2.logp);                // encIP2 = a/b + yWTx*yWTx

				NTL_EXEC_RANGE(cnum, first, last);
				//long first = 0, last = cnum;
				for (long i = first; i < last; ++i) {

					Ciphertext encGradtemp;

					scheme.multByConst(encGradtemp, encZData[r * cnum + i], (gamma) * degree3[2], wBits + pBits);

					scheme.reScaleByAndEqual(encGradtemp, pBits); // encGrad = Y@X *gamma * b

					Ciphertext ctIP(encIP[r]);
					if (encGradtemp.logq > ctIP.logq)
						scheme.modDownToAndEqual(encGradtemp, ctIP.logq); /* whose logq should be ... */
					if (encGradtemp.logq < ctIP.logq)
						scheme.modDownToAndEqual(ctIP, encGradtemp.logq);

					scheme.multAndEqual(encGradtemp, ctIP); // encGrad = gamma * Y@X * b * yWTx
					scheme.reScaleByAndEqual(encGradtemp, ctIP.logp);

					Ciphertext ctIP2(encIP2);
					if (encGradtemp.logq > ctIP2.logq)
						scheme.modDownToAndEqual(encGradtemp, ctIP2.logq);
					if (encGradtemp.logq < ctIP2.logq)
						scheme.modDownToAndEqual(ctIP2, encGradtemp.logq);
					scheme.multAndEqual(encGradtemp, ctIP2);// encGrad = gamma * Y@X * (a * yWTx + b * yWTx ^3)
					scheme.reScaleByAndEqual(encGradtemp, ctIP2.logp);

					Ciphertext tmp;
					scheme.multByConst(tmp, encZData[r * cnum + i], (gamma) * degree3[0], wBits); // tmp = Y@X * gamma * 0.5

					scheme.modDownToAndEqual(tmp, encGradtemp.logq);// encGrad[i].logq == tmp.logq

					scheme.addAndEqual(encGradtemp, tmp);// encGrad = gamma * Y@X * (0.5 + a * yWTx + b * yWTx ^3)



					if (encGrad[i].logp > encGradtemp.logp)
						scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp - encGradtemp.logp);
					if (encGrad[i].logp < encGradtemp.logp)
						scheme.reScaleByAndEqual(encGradtemp, encGradtemp.logp - encGrad[i].logp);
					if (encGrad[i].logq > encGradtemp.logq)
						scheme.modDownToAndEqual(encGrad[i], encGradtemp.logq);
					if (encGrad[i].logq < encGradtemp.logq)
						scheme.modDownToAndEqual(encGradtemp, encGrad[i].logq);
					if (encGrad[i].logp != encGradtemp.logp) {cout << "logp != logp"; exit(0);}

					scheme.addAndEqual(encGrad[i], encGradtemp);

					tmp.kill();
					ctIP2.kill();
					ctIP.kill();
					encGradtemp.kill();
				}
				NTL_EXEC_RANGE_END
				;
				/* END OF if(kdeg == 3) {  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

				encIP2.kill();
				//encIP.kill();
			}
			//NTL_EXEC_RANGE_END Bugs Found Here ！

			for (long kk = 0; kk < rnum; ++kk) encIP[kk].kill();

			// Sum Each Column of encGrad[i] To Get the Final gradient : (1 - sigm(yWTx)) * Y.T @ X
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {
				Ciphertext tmp;
				for (long l = bBits; l < sBits; ++l) {
					scheme.leftRotateFast(tmp, encGrad[i], (1 << l));
					scheme.addAndEqual(encGrad[i], tmp);
				}

				Ciphertext ctBinv(encBinv[i]); ////~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~
				if (encGrad[i].logq > ctBinv.logq)
					scheme.modDownToAndEqual(encGrad[i], ctBinv.logq);
				if (encGrad[i].logq < ctBinv.logq)
					scheme.modDownToAndEqual(ctBinv, encGrad[i].logq);

				scheme.multAndEqual(encGrad[i], ctBinv);
				scheme.reScaleByAndEqual(encGrad[i], ctBinv.logp);

				tmp.kill();
				ctBinv.kill();
			}
			NTL_EXEC_RANGE_END
			;
			/* Each ([i][0~batch)-th) column of encGrad[i] consists of the same value (gamma * encGrad[i][0~batch)) */
			/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */


			cout << endl << "Quadratic Gradient: " << endl ;
			for (long c = 0; c < cnum; ++c) {
				complex<double> *dcvdddddv = scheme.decrypt(secretKey, encGrad[c]);
				for (long j = 0; j < batch; ++j)
					cout << setiosflags(ios::fixed) << setprecision(6) << dcvdddddv[j].real() << "\t";

				delete[] dcvdddddv;
			}



			/* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {


				if (encGrad[i].logp > encVData[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp - encVData[i].logp);
				if (encGrad[i].logp < encVData[i].logp) scheme.reScaleByAndEqual(encVData[i], encVData[i].logp - encGrad[i].logp);
				scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);

				Ciphertext ctmpw;
				scheme.add(ctmpw, encVData[i], encGrad[i]); // encGrad[i] has already self-multiplied with gamma
				// ctmpw = encVData[i] - encGrad[i]

				scheme.multByConst(encVData[i], ctmpw, 1. - eta, pBits);// encVData[i] = ( 1. - eta ) * ctmpw
				//scheme.reScaleByAndEqual(encVData[i], pBits-5);

				scheme.multByConstAndEqual(encWData[i], eta, pBits);// encWData[i] = eta * encWData[i]
				//scheme.reScaleByAndEqual(encWData[i], pBits-5);

				if (encWData[i].logq > encVData[i].logq) scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
				if (encWData[i].logq < encVData[i].logq) scheme.modDownToAndEqual(encVData[i], encWData[i].logq);
				if (encWData[i].logp != encVData[i].logp) {cout << "logp != logp"; exit(0);}

				scheme.addAndEqual(encVData[i], encWData[i]); // encVData[i] = encVData[i] + encWData[i]
				// encVData[i] = ( 1. - eta ) * ctmpw + eta * encWData[i]

				scheme.reScaleByAndEqual(encVData[i], pBits);
				encWData[i].copy(ctmpw);

				ctmpw.kill();
			}
			NTL_EXEC_RANGE_END

			;
			/* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */

	auto end_timeUpdatingWeights = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end_timeUpdatingWeights - start_timeUpdatingWeights;
        totalUpdatingWeightsTime += duration.count();  // Returns the elapsed time in seconds
        totalUpdatingWeights += 1;


			for (long i = 0; i < cnum; ++i)
				encGrad[i].kill();
			delete[] encGrad;
			/* cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits); */

			timeutils.stop("NesterovWithG : " + to_string(iter + 1)  + " -th iteration");

			openFileTIME << "," << timeutils.timeElapsed;
			openFileTIME.flush();
			openFileTIMELabel << ","
			                  << "NesterovWithGminBatch : " + to_string(iter + 1)
			                  + " -th epoch." + to_string(r + 1)
			                  + "-th iteration";
			openFileTIMELabel.flush();
			openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
			openFileCurrMEM.flush();
			openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
			openFilePeakMEM.flush();

			for (long i = 0; i < rnum * cnum; ++i)
				encZData[i].kill();
			delete[] encZData;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			cout << endl << "---------- TEST : THE " << iter * rnum + r + 1
			     << "-th ITERATION : Weights, AUC, MLE ----------" << endl;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);     */
			for (long i = 0; i < (cnum - 1); ++i) {
				complex<double> *dcvv = scheme.decrypt(secretKey, encVData[i]);
				for (long j = 0; j < batch; ++j) {
					cvData[batch * i + j] = dcvv[j].real();
				}
				delete[] dcvv;
			}
			complex<double> *dcvv = scheme.decrypt(secretKey,
			                                       encVData[cnum - 1]);
			long rest = factorDim - batch * (cnum - 1);
			for (long j = 0; j < rest; ++j) {
				cvData[batch * (cnum - 1) + j] = dcvv[j].real();
			}
			delete[] dcvv;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits); */
			cout << "Current cWdata (encVData) : " << endl;
			for (long i = 0; i < factorDim; ++i)
				cout << setiosflags(ios::fixed) << setprecision(12) << cvData[i]
				     << ",\t";
			cout << endl;

			openFileTestAUC << ","
			                << MyTools::calculateAUC(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestAUC.flush();
			openFileTrainAUC << ","
			                 << MyTools::calculateAUC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainAUC.flush();

			openFileTestACC << ","
			                << MyTools::calculateACC(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestACC.flush();
			openFileTrainACC << ","
			                 << MyTools::calculateACC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainACC.flush();

			openFileTestMLE << ","
			                << MyTools::calculateMLE(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestMLE.flush();
			openFileTrainMLE << ","
			                 << MyTools::calculateMLE(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainMLE.flush();

			cout << endl << endl << endl << endl << endl;

			cout << "TrainMLE : "
			     << MyTools::calculateMLE(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << "TrainAUC : "
			     << MyTools::calculateAUC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << "TrainACC : "
			     << MyTools::calculateACC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << endl;
			cout << "Test MLE : "
			     << MyTools::calculateMLE(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << "Test AUC : "
			     << MyTools::calculateAUC(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << "Test ACC : "
			     << MyTools::calculateACC(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << " ----------------- end of the " << (iter + 1)
			     << "-th epoch ----------------- " << endl;
			cout << " ----------------- end of the " << (r + 1)
			     << "-th iteration ----------------- " << endl << endl;
			cout
			        << "--------------------------------------------------------------------------------"
			        << endl;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			cout << endl << endl << endl << endl << "encVData[0].logq = "
			     << encVData[0].logq << endl << endl << endl << endl;





			if (encVData[0].logq < 220 + encVData[0].logp) {
				//if ( encVData[0].logq <= 450 && iter < Epoch_Number-1 || encVData[0].logq < wBits && iter == Epoch_Number-1) {

				timeutils.start("Use Bootstrap To Recrypt Ciphertext");
				auto start_timeRefreshingWeights = std::chrono::high_resolution_clock::now();

				NTL_EXEC_RANGE(cnum, first, last);
				for (long i = first; i < last; ++i) {
					scheme.modDownToAndEqual(encWData[i], bootlogq);
					encWData[i].n = batch;
					scheme.bootstrapAndEqual(encWData[i], bootlogq, 990, logT, logI);
					encWData[i].n = slots;
				}
				NTL_EXEC_RANGE_END
				NTL_EXEC_RANGE(cnum, first, last);
				for (long i = first; i < last; ++i) {
					scheme.modDownToAndEqual(encVData[i], bootlogq);
					encVData[i].n = batch;
					scheme.bootstrapAndEqual(encVData[i], bootlogq, 990, logT, logI);
					encVData[i].n = slots;
				}
				NTL_EXEC_RANGE_END

				auto end_timeRefreshingWeights = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duretion = end_timeRefreshingWeights - start_timeRefreshingWeights;
        totalRefreshingWeightsTime += duretion.count();  // Returns the elapsed time in seconds
        totalBootstrappingOperations += 1;



				timeutils.stop("Use Bootstrap To Recrypt Ciphertext");

				openFileTIME << "," << timeutils.timeElapsed;
				openFileTIME.flush();
				openFileTIMELabel << "," << "Bootstrapping";
				openFileTIMELabel.flush();
				openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
				openFileCurrMEM.flush();
				openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
				openFilePeakMEM.flush();

			}
			/////////////////////////////////////////////////////////////////////////////
			//        BOOTSTRAPPING                                                    //
			//             Over and Out                                                //
			/////////////////////////////////////////////////////////////////////////////
			cout
			        << "--------------------------------------------------------------------------------"
			        << endl << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << endl << " !!! STOP " << iter + 1 << " ITERATION !!! "
			     << endl << endl;
			     			     	cout << endl << endl << endl << endl;
			     	cout << "totalUpdatingWeightsTime = " << totalUpdatingWeightsTime << endl;
			     	cout << "totalUpdatingWeights = " << totalUpdatingWeights << endl;
				cout << "totalRefreshingWeightsTime = " << totalRefreshingWeightsTime << endl;
				cout << "totalBootstrappingOperations = " << totalBootstrappingOperations << endl;
				cout << endl << endl << endl << endl;
		}
	}

	openFileTIME << endl;
	openFileTIME.flush();
	openFileTIMELabel << endl;
	openFileTIMELabel.flush();
	openFileTestAUC << endl;
	openFileTestAUC.flush();
	openFileTrainAUC << endl;
	openFileTrainAUC.flush();
	openFileTestACC << endl;
	openFileTestACC.flush();
	openFileTrainACC << endl;
	openFileTrainACC.flush();
	openFileTestMLE << endl;
	openFileTestMLE.flush();
	openFileTrainMLE << endl;
	openFileTrainMLE.flush();
	openFileCurrMEM << endl;
	openFileCurrMEM.flush();
	openFilePeakMEM << endl;
	openFilePeakMEM.flush();

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




double* MyMethods::testCryptoMiniBatchNAGwithG(double **traindata,
        double *trainlabel, long factorDim, long trainSampleDim, long Epoch_Number,
        double **testdata, double *testlabel, long testSampleDim,
        string resultpath) {

	long wBits = 30;
	long pBits = 20;

	//long logN = MyTools::suggestLogN(80, logQ);
	long logN = 16;
	long logQ = 275; // 991.300840336 > logQ  to ensure 128-bit security level. Security Parameter λ
	long slots = 1 << (logN - 1);
	long sBits = (long) ceil(log2(slots));
	long fdimBits = (long) ceil(log2(factorDim));
	long sdimBits = (long) ceil(log2(trainSampleDim));

	long batch = 1 << 5; // Basically, batch is the Number of several factor dimensions.
	long bBits = (long) ceil(log2(batch)); // 2^batchBits = min( 2^logN / 2^sdimBits / 2, 2^fdimBits ) = min( N/2 /n, factorDim ) ;

	// the size of batch should be small than fatctorDim
	//if (batch >= factorDim) {
	//	cout << "batch >= factorDim!" << endl;
	//	exit(0);
	//}

	//long factorNum = 1 << (long)ceil(log2(factorDim));
	long cnum = (long) ceil(double(factorDim) / batch); // To Divide the whole Train Data into Several Batches (cnum Ciphertexts).
	if (cnum > (1 << fdimBits)) {
		cout << "cnum should be no more than factorDim!" << endl;
		exit(0);
	}

// We assume that the dataset is too large to be encrypted into one single ciphertext.
	long minbatchsize = slots / batch; // (each min-batch should be small to save the useless space of the last batch)
	long minBatchDimBits = (long) ceil(log2(minbatchsize));

	long rnum = (long) ceil((double) trainSampleDim / minbatchsize); // rnum : the number of min-batch, namely number of row

	cout << "logQ = " << logQ << ", logN = " << logN << endl;
	cout << "slots = " << slots << ", batch = " << batch << endl;
	cout << "rnum = " << rnum << ", cnum = " << cnum << endl;
	cout << "min-batch size = " << minbatchsize << ", minBatchDimBits = "
	     << minBatchDimBits << endl;

	string path = resultpath;
	ofstream openFileTrainAUC(path + "TrainAUC.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainACC(path + "TrainACC.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTrainMLE(path + "TrainMLE.csv",
	                          std::ofstream::out | std::ofstream::app);
	ofstream openFileTestAUC(path + "TestAUC.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTestACC(path + "TestACC.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTestMLE(path + "TestMLE.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFileTIME(path + "TIME.csv",
	                      std::ofstream::out | std::ofstream::app);
	ofstream openFileTIMELabel(path + "TIMELabel.csv",
	                           std::ofstream::out | std::ofstream::app);
	ofstream openFileCurrMEM(path + "CurrMEM.csv",
	                         std::ofstream::out | std::ofstream::app);
	ofstream openFilePeakMEM(path + "PeakMEM.csv",
	                         std::ofstream::out | std::ofstream::app);

	if (!openFileTrainAUC.is_open())
		cout << "Error: cannot read Train AUC file" << endl;
	if (!openFileTrainACC.is_open())
		cout << "Error: cannot read Train ACC file" << endl;
	if (!openFileTrainMLE.is_open())
		cout << "Error: cannot read Train MLE file" << endl;
	if (!openFileTestAUC.is_open())
		cout << "Error: cannot read Test AUC file" << endl;
	if (!openFileTestACC.is_open())
		cout << "Error: cannot read Test ACC file" << endl;
	if (!openFileTestMLE.is_open())
		cout << "Error: cannot read Test MLE file" << endl;
	if (!openFileTIME.is_open())
		cout << "Error: cannot read TIME file" << endl;
	if (!openFileTIMELabel.is_open())
		cout << "Error: cannot read TIME Label file" << endl;
	if (!openFileTIMELabel.is_open())
		cout << "Error: cannot read TIME Label file" << endl;
	if (!openFileCurrMEM.is_open())
		cout << "Error: cannot read Current MEMORY file" << endl;
	if (!openFilePeakMEM.is_open())
		cout << "Error: cannot read Peak MEMORY file" << endl;

	openFileTrainAUC << "TrainAUC";
	openFileTrainAUC.flush();
	openFileTrainACC << "TrainACC";
	openFileTrainACC.flush();
	openFileTrainMLE << "TrainMLE";
	openFileTrainMLE.flush();
	openFileTestAUC << "TestAUC";
	openFileTestAUC.flush();
	openFileTestACC << "TestACC";
	openFileTestACC.flush();
	openFileTestMLE << "TestMLE";
	openFileTestMLE.flush();
	openFileTIME << "TIME";
	openFileTIME.flush();
	openFileTIMELabel << "TIMELabel";
	openFileTIMELabel.flush();
	openFileCurrMEM << "MEMORY(GB)";
	openFileCurrMEM.flush();
	openFilePeakMEM << "MEMORY(GB)";
	openFilePeakMEM.flush();

	TimeUtils timeutils;
	timeutils.start("Scheme generating");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Scheme generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	cout << "NOW THE PEAK RSS IS: " << (MyTools::getPeakRSS() >> 20) << endl;

	timeutils.start("Bootstrap Key generating");
	long logT = 3;
	long logI = 4;
	long bootlogq = 30 + 10;
	long lognslots = (long) ceil(log2(batch));  //batch = factorDim / cnum;
	scheme.addBootKey(secretKey, lognslots, bootlogq + 4);
	timeutils.stop("Bootstrap Key generated");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Bootstrap Key generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	cout << "NOW THE CURR RSS IS: " << (MyTools::getCurrentRSS() >> 20) << endl;
	cout << "NOW THE PEAK RSS IS: " << (MyTools::getPeakRSS() >> 20) << endl;

	// Basically, rpoly is used to calculate the sum{row}
	timeutils.start("Polynomial generating...");
	ZZ *dummy = new ZZ[N];
	;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	complex<double> *pvals = new complex<double> [slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j] = 1.0;
	}
	scheme.ring.encode(dummy, pvals, slots, pBits);
	delete[] pvals;
	/* cipherGD.generateAuxPoly(rpoly, slots, batch, pBits); */
	timeutils.stop("Polynomial generation");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Polynomial generating";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	double *cwData = new double[factorDim]();
	double *cvData = new double[factorDim]();

	Ciphertext *encTrainData = new Ciphertext[rnum * cnum];
	Ciphertext *encTrainLabel = new Ciphertext[rnum * cnum];

	Ciphertext *encXyZdata = new Ciphertext[rnum * cnum]; // TrainLabel @ TrainData = yX
	Ciphertext *encBinv = new Ciphertext[rnum * cnum];

	Ciphertext *encWData = new Ciphertext[cnum];
	Ciphertext *encVData = new Ciphertext[cnum];

	/* - - - - - - - - - - - - - - - - - - - - - - - - Client and Server - - - - - - - - - - - - - - - - - - - - - - - - */

	timeutils.start("Encrypting encXyZdata...");
	// encrypt the traindata
	for (long r = 0; r < rnum - 1; ++r) {
		for (long i = 0; i < cnum - 1; ++i) {

			complex<double> *pzData = new complex<double> [slots]();
			for (long j = 0; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(
					    trainlabel[r * minbatchsize + j]
					    * traindata[r * minbatchsize + j][batch * i
					                                      + l]);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encXyZdata[r * cnum + i], pzData, slots, wBits,
			               logQ);
			SerializationUtils::writeCiphertext(encXyZdata[r * cnum + i], "encXyZdata[" + std::to_string(r * cnum + i) + "].txt");
		}
		// i == cnum - 1       - the last cnum in each row
		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < minbatchsize; ++j) {
			long rest = factorDim - batch * (cnum - 1);

			for (long l = 0; l < rest; ++l) {
				pzData[batch * j + l].real(
				    trainlabel[r * minbatchsize + j]
				    * traindata[r * minbatchsize + j][batch
				                                      * (cnum - 1) + l]);
				pzData[batch * j + l].imag(0);
			}
			for (long l = rest; l < batch; ++l) {
				pzData[batch * j + l].real(0);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encXyZdata[r * cnum + cnum - 1], pzData, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encXyZdata[r * cnum + cnum - 1], "encXyZdata[" + std::to_string(r * cnum + cnum - 1) + "].txt");

	}
	// The last min-batch may consists of several ( trainSampleDim - minbatchsize * (rnum-1) ) rows of zeors.

	// r == rnum - 1       - the last rnum (the last min-batch)
	auto restrownum = trainSampleDim - minbatchsize * (rnum - 1);
	for (long i = 0; i < cnum - 1; ++i) {

		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < restrownum; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(
				    trainlabel[(rnum - 1) * minbatchsize + j]
				    * traindata[(rnum - 1) * minbatchsize + j][batch
				            * i + l]);
				pzData[batch * j + l].imag(0);
			}
		}
		for (long j = restrownum; j < minbatchsize; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(0);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encXyZdata[(rnum - 1) * cnum + i], pzData, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encXyZdata[(rnum - 1) * cnum + i], "encXyZdata[" + std::to_string((rnum - 1) * cnum + i) + "].txt");
	}
	// i == cnum - 1       - the last cnum in each row
	complex<double> *pzDatra = new complex<double> [slots]();
	for (long j = 0; j < restrownum; ++j) {
		long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < rest; ++l) {
			pzDatra[batch * j + l].real(
			    trainlabel[(rnum - 1) * minbatchsize + j]
			    * traindata[(rnum - 1) * minbatchsize + j][batch
			            * (cnum - 1) + l]);
			pzDatra[batch * j + l].imag(0);
		}
		for (long l = rest; l < batch; ++l) {
			pzDatra[batch * j + l].real(0);
			pzDatra[batch * j + l].imag(0);
		}
	}
	for (long j = restrownum; j < minbatchsize; ++j) {
		//long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < batch; ++l) {
			pzDatra[batch * j + l].real(0);
			pzDatra[batch * j + l].imag(0);
		}
	}
	scheme.encrypt(encXyZdata[(rnum - 1) * cnum + cnum - 1], pzDatra, slots,
	               wBits, logQ);
	SerializationUtils::writeCiphertext(encXyZdata[(rnum - 1) * cnum + cnum - 1], "encXyZdata[" + std::to_string((rnum - 1) * cnum + cnum - 1) + "].txt");
	delete[] pzDatra;
	timeutils.stop("encXyZdata encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting encXyZdata";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	double **Binv = MyTools::zInvBFromFile(traindata, factorDim,
	                                       trainSampleDim);
	for (long r = 0; r < rnum - 1; ++r) {
		double **zInvB = new double*[minbatchsize];
		for (long i = 0; i < minbatchsize; ++i) {
			double *Bj = new double[factorDim]();
			for (long j = 0; j < factorDim; ++j)
				Bj[j] = traindata[r * minbatchsize + i][j];
			zInvB[i] = Bj;
		}
		auto zTemp = MyTools::zInvBFromFile(zInvB, factorDim, minbatchsize);

		for (long i = 0; i < minbatchsize; ++i) {
			for (long j = 0; j < factorDim; ++j)
				Binv[r * minbatchsize + i][j] = min(zTemp[i][j], 8.); // overflow:: beyond the precision!!
		}
		delete[] zInvB;
		delete[] zTemp;
	}
	auto restrows = trainSampleDim - minbatchsize * (rnum - 1);
	double **zInvB = new double*[restrows];
	for (long i = 0; i < restrows; ++i) {
		double *Bj = new double[factorDim]();
		for (long j = 0; j < factorDim; ++j)
			Bj[j] = traindata[(rnum - 1) * minbatchsize + i][j];
		zInvB[i] = Bj;
	}
	auto zTemp = MyTools::zInvBFromFile(zInvB, factorDim, restrows);

	for (long i = 0; i < restrows; ++i) {
		for (long j = 0; j < factorDim; ++j)
			Binv[(rnum - 1) * minbatchsize + i][j] = min(zTemp[i][j], 8.); // overflow:: beyond the precision!!
	}
	delete[] zInvB;
	delete[] zTemp;
	timeutils.start("Encrypting Binv...");
	// encrypt the traindata
	for (long r = 0; r < rnum - 1; ++r) {
		for (long i = 0; i < cnum - 1; ++i) {

			complex<double> *pzData = new complex<double> [slots]();
			for (long j = 0; j < minbatchsize; ++j) {
				for (long l = 0; l < batch; ++l) {
					pzData[batch * j + l].real(
					    Binv[r * minbatchsize + j][batch * i + l]);
					pzData[batch * j + l].imag(0);
				}
			}
			scheme.encrypt(encBinv[r * cnum + i], pzData, slots, wBits, logQ);
			SerializationUtils::writeCiphertext(encBinv[r * cnum + i], "encBinv[" + std::to_string(r * cnum + i) + "].txt");
		}
		// i == cnum - 1       - the last cnum in each row
		complex<double> *pzData3 = new complex<double> [slots]();
		for (long j = 0; j < minbatchsize; ++j) {
			long rest = factorDim - batch * (cnum - 1);

			for (long l = 0; l < rest; ++l) {
				pzData3[batch * j + l].real(
				    Binv[r * minbatchsize + j][batch * (cnum - 1) + l]);
				pzData3[batch * j + l].imag(0);
			}
			for (long l = rest; l < batch; ++l) {
				pzData3[batch * j + l].real(0);
				pzData3[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encBinv[r * cnum + cnum - 1], pzData3, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encBinv[r * cnum + cnum - 1], "encBinv[" + std::to_string(r * cnum + cnum - 1) + "].txt");

	}
	// The last min-batch may consists of several ( trainSampleDim - minbatchsize * (rnum-1) ) rows of zeors.

	// r == rnum - 1       - the last rnum (the last min-batch)
	restrownum = trainSampleDim - minbatchsize * (rnum - 1);
	for (long i = 0; i < cnum - 1; ++i) {

		complex<double> *pzData = new complex<double> [slots]();
		for (long j = 0; j < restrownum; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(
				    Binv[(rnum - 1) * minbatchsize + j][batch * i + l]);
				pzData[batch * j + l].imag(0);
			}
		}
		for (long j = restrownum; j < minbatchsize; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(0);
				pzData[batch * j + l].imag(0);
			}
		}
		scheme.encrypt(encBinv[(rnum - 1) * cnum + i], pzData, slots, wBits,
		               logQ);
		SerializationUtils::writeCiphertext(encBinv[(rnum - 1) * cnum + i], "encBinv[" + std::to_string((rnum - 1) * cnum + i) + "].txt");
	}
	// i == cnum - 1       - the last cnum in each row
	complex<double> *pzData7 = new complex<double> [slots]();
	for (long j = 0; j < restrownum; ++j) {
		long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < rest; ++l) {
			pzData7[batch * j + l].real(
			    Binv[(rnum - 1) * minbatchsize + j][batch * (cnum - 1) + l]);
			pzData7[batch * j + l].imag(0);
		}
		for (long l = rest; l < batch; ++l) {
			pzData7[batch * j + l].real(0);
			pzData7[batch * j + l].imag(0);
		}
	}
	for (long j = restrownum; j < minbatchsize; ++j) {
		//long rest = factorDim - batch * (cnum - 1);
		for (long l = 0; l < batch; ++l) {
			pzData7[batch * j + l].real(0);
			pzData7[batch * j + l].imag(0);
		}
	}
	scheme.encrypt(encBinv[(rnum - 1) * cnum + cnum - 1], pzData7, slots, wBits,
	               logQ);
	SerializationUtils::writeCiphertext(encBinv[(rnum - 1) * cnum + cnum - 1], "encBinv[" + std::to_string((rnum - 1) * cnum + cnum - 1) + "].txt");
	delete[] pzData7;
	timeutils.stop("Binv encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting Binv";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();

	// zDataTrain and zDataTest are used for testing in the plaintext environment
	// zData = (Y,Y@X)
	double **zDataTrain = new double*[trainSampleDim];
	for (int i = 0; i < trainSampleDim; ++i) {
		double *zi = new double[factorDim]();
		zi[0] = trainlabel[i];
		for (int j = 1; j < factorDim; ++j)
			zi[j] = zi[0] * traindata[i][j];
		zDataTrain[i] = zi;
	}
	// zDataTest is only used for Cross-Validation test, not necessary for training LG model.
	// zData = (Y,Y@X)
	double **zDataTest = new double*[testSampleDim];
	for (int i = 0; i < testSampleDim; ++i) {
		double *zi = new double[factorDim]();
		zi[0] = testlabel[i];
		for (int j = 1; j < factorDim; ++j)
			zi[j] = zi[0] * testdata[i][j];
		zDataTest[i] = zi;
	}

	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */
	timeutils.start("Encrypting wData and vData...");
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		// scheme.encryptZeros(encWData[i], slots, wBits, encZData[0].logq); // To Make encVData[0].logq==encZData[0].logq
		scheme.encryptSingle(encWData[i], 0.0, wBits, logQ);
		encWData[i].n = slots;

		encVData[i].copy(encWData[i]);
	}
	NTL_EXEC_RANGE_END
	;
	timeutils.stop("wData and vData encryption");

	openFileTIME << "," << timeutils.timeElapsed;
	openFileTIME.flush();
	openFileTIMELabel << "," << "Encrypting weight vector vData";
	openFileTIMELabel.flush();
	openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
	openFileCurrMEM.flush();
	openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
	openFilePeakMEM.flush();
	/* cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);  */

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                                                                                                                    //
	//                        Client sent (encTrainData, encTrainLabel, enc(x0), encWData, and encVData) to Server                        //
	//                                                                                                                                    //
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// From now on, the server starts its work on what client sent to it. //////////////////////////////////
	cout << endl << endl << endl << endl << "scheme.leftRotKeyMap.size() = " << scheme.leftRotKeyMap.size() << endl << endl << endl << endl;
	cout << endl << endl << endl << endl;
	double min_lr = 1.0;
	double max_lr = 2.0;
	double total_steps = Epoch_Number * rnum;
	double exp_gamma = 2.5;
	cout << "min_lr = " << min_lr << endl;
	cout << "max_lr = " << max_lr << endl;
	cout << "total_steps = " << total_steps << endl;
	cout << "exp_gamma = " << exp_gamma << endl;
	cout << "len(Xminbatches) = rnum = " << rnum << endl;
	cout << endl << endl << endl << endl;

	double totalUpdatingWeightsTime = 0.0;
	double totalUpdatingWeights = 0.0;
	double totalRefreshingWeightsTime = 0.0;
	double totalBootstrappingOperations = 0.0;




	double alpha0, alpha1, eta, gamma;
	double enccor, encauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	cout << endl << "unsigned seed = chrono::system_clock::now().time_since_epoch().count();";
	cout << endl << "unsigned seed = " << (unsigned)seed << endl;
	
	for (long iter = 0; iter < Epoch_Number; ++iter) {

		vector<int> randr;
		for (int ir = 0; ir < rnum; ++ir)
			randr.push_back(ir);
		cout << "randr[]: ";
		for (int &x : randr)
			cout << ' ' << x;
		cout << '\n';
		
		shuffle(randr.begin(), randr.end(), default_random_engine(seed));
		cout << "randr[]: ";
		for (int &x : randr)
			cout << ' ' << x;
		cout << '\n';

		//cout << endl << "NesterovWithG : "+ to_string(iter+1)+" -th iteration" << endl;
		for (long r = 0; r < rnum; ++r) {
			timeutils.start("NesterovWithGminBatch : " + to_string(iter + 1) + " -th iteration");
			auto start_timeUpdatingWeights = std::chrono::high_resolution_clock::now();

			cout << endl << endl << endl << endl << "encVData[0].logq = "
			     << encVData[0].logq << endl << endl << endl << endl;

			     
			eta = (1 - alpha0) / alpha1;

			auto iterations = iter * rnum + r;
			// self.max_lr - (self.max_lr - self.min_lr) * (self.current_step / self.total_steps)**self.gamma
			auto gamma = max_lr - (max_lr - min_lr) * pow(iterations / total_steps, exp_gamma );



			/* cipherGD.encZData(encZData, zDataTrain, slots, factorDim, trainSampleDim, batch, cnum, wBits, logQ);  */
			Ciphertext *encZData = new Ciphertext[cnum];
			//timeutils.start("encZData[i] = encTrainLabel[i] @ encTrainData[i] for i in range(cnum)");
			// To Get the encZData
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {
				//encZData[i].copy(encXyZdata[randr[r] * cnum + i]);
				encZData[i].copy(encXyZdata[r * cnum + i]);
			}
			NTL_EXEC_RANGE_END
			//timeutils.stop("encZData[i] for i in range(cnum) is done");
			/* cipherGD.encZData(encZData, zDataTrain, slots, factorDim, trainSampleDim, batch, cnum, wBits, logQ);  */

			/* cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits); */

			/* CipherGD::encInnerProduct(encIP, encZData, encWData, rpoly, cnum, bBits, wBits, pBits); */
			Ciphertext *encIPvec = new Ciphertext[cnum];

			/* For Each Batch, Sum Itself Inside */
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {
				// MAKE SURE : encZData[i].logq >= encVData[i].logq (and of course : encZData[i].logp == encVData[i].logp)
				encIPvec[i].copy(encZData[i]);
				if (encIPvec[i].logq > encVData[i].logq) { //cout << "FF" << endl;
					scheme.modDownToAndEqual(encIPvec[i], encVData[i].logq);
				}
				if (encIPvec[i].logq < encVData[i].logq) { //cout << "EE" << endl;
					scheme.modDownToAndEqual(encVData[i], encIPvec[i].logq);
				}
				// V is the final weights to store the result weights.
				scheme.multAndEqual(encIPvec[i], encVData[i]);// encIPvec = ENC(zData) .* ENC(V)

			}
			NTL_EXEC_RANGE_END

			/* Sum All Batchs To Get One Batch */
			Ciphertext encIP;
			encIP.copy(encIPvec[0]);          // to store the sum of all batches
			for (long i = 1; i < cnum; ++i) {
				scheme.addAndEqual(encIP, encIPvec[i]);
			}

			/* For Each Batch (==ciphertext), Sum Itself Inside, Result in Each Row consisting of the same value */
			Ciphertext rot;                   // encIPvec = ENC(zData) @  ENC(V)
			for (long l = 0; l < bBits; ++l) {
				scheme.leftRotateFast(rot, encIP, (1 << l));
				scheme.addAndEqual(encIP, rot);
			}
			rot.kill();

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

			for (long i = 0; i < cnum; ++i)
				encIPvec[i].kill();
			delete[] encIPvec;
			//CipherGD::encInnerProduct(encIP, encZData, encVData, rpoly, cnum, bBits, wBits, pBits);

			Ciphertext *encGrad = new Ciphertext[cnum];

			/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */

			Ciphertext encIP2(encIP);
			scheme.multAndEqual(encIP2, encIP);

			// IT IS VERY IMPORT TO KEEP THE logp BIG ENOUGH TO PRESENT THE NUMBER !!  MAKE SURE encIP2.logp>=35
			scheme.reScaleByAndEqual(encIP2, encIP.logp); // For now, encIP.logp is big enough

			//////////////////////////////////////// when iteration < 05 ////////////////////////////////////////
			cout << endl << "INSIDE iter < 5;  poly3 = ";
			cout << setiosflags(ios::showpos) << degree3[0] << " ";
			cout << setiosflags(ios::showpos) << degree3[1] << "x ";
			cout << setiosflags(ios::showpos) << degree3[2] << "x^3 " << endl
			     << endl;
			cout << std::noshowpos;
			cout << "gamma = " << gamma << endl << endl << endl;

			scheme.addConstAndEqual(encIP2, degree3[1] / degree3[2],
			                        encIP2.logp);                // encIP2 = a/b + yWTx*yWTx

			NTL_EXEC_RANGE(cnum, first, last);
			//long first = 0, last = cnum;
			for (long i = first; i < last; ++i) {

				scheme.multByConst(encGrad[i], encZData[i], (gamma) * degree3[2], wBits + pBits);

				scheme.reScaleByAndEqual(encGrad[i], pBits); // encGrad = Y@X *gamma * b

				Ciphertext ctIP(encIP);
				if (encGrad[i].logq > ctIP.logq)
					scheme.modDownToAndEqual(encGrad[i], ctIP.logq); /* whose logq should be ... */
				if (encGrad[i].logq < ctIP.logq)
					scheme.modDownToAndEqual(ctIP, encGrad[i].logq);

				scheme.multAndEqual(encGrad[i], ctIP); // encGrad = gamma * Y@X * b * yWTx
				scheme.reScaleByAndEqual(encGrad[i], ctIP.logp);

				Ciphertext ctIP2(encIP2);
				if (encGrad[i].logq > ctIP2.logq)
					scheme.modDownToAndEqual(encGrad[i], ctIP2.logq);
				if (encGrad[i].logq < ctIP2.logq)
					scheme.modDownToAndEqual(ctIP2, encGrad[i].logq);
				scheme.multAndEqual(encGrad[i], ctIP2);// encGrad = gamma * Y@X * (a * yWTx + b * yWTx ^3)
				scheme.reScaleByAndEqual(encGrad[i], ctIP2.logp);

				Ciphertext tmp;
				scheme.multByConst(tmp, encZData[i], (gamma) * degree3[0], wBits);// tmp = Y@X * gamma * 0.5

				scheme.modDownToAndEqual(tmp, encGrad[i].logq);// encGrad[i].logq == tmp.logq

				scheme.addAndEqual(encGrad[i], tmp);// encGrad = gamma * Y@X * (0.5 + a * yWTx + b * yWTx ^3)

				tmp.kill();
				ctIP2.kill();
				ctIP.kill();
			}
			NTL_EXEC_RANGE_END
			;
			/* END OF if(kdeg == 3) {  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */

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

				Ciphertext ctBinv(encBinv[r * cnum + i]); ////~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~//~~~~~~
				if (encGrad[i].logq > ctBinv.logq)
					scheme.modDownToAndEqual(encGrad[i], ctBinv.logq);
				if (encGrad[i].logq < ctBinv.logq)
					scheme.modDownToAndEqual(ctBinv, encGrad[i].logq);

				scheme.multAndEqual(encGrad[i], ctBinv);
				scheme.reScaleByAndEqual(encGrad[i], ctBinv.logp);

				tmp.kill();
				ctBinv.kill();
			}
			NTL_EXEC_RANGE_END
			;
			/* Each ([i][0~batch)-th) column of encGrad[i] consists of the same value (gamma * encGrad[i][0~batch)) */
			/* CipherGD::encSigmoid(kdeg, encZData, encGrad, encIP, cnum, gamma, sBits, bBits, wBits, aBits); */

			/* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */
			NTL_EXEC_RANGE(cnum, first, last);
			for (long i = first; i < last; ++i) {

				if (encGrad[i].logp > encVData[i].logp) scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp - encVData[i].logp);
				if (encGrad[i].logp < encVData[i].logp) scheme.reScaleByAndEqual(encVData[i], encVData[i].logp - encGrad[i].logp);
				scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);

				Ciphertext ctmpw;
				scheme.add(ctmpw, encVData[i], encGrad[i]); // encGrad[i] has already self-multiplied with gamma
				// ctmpw = encVData[i] - encGrad[i]

				scheme.multByConst(encVData[i], ctmpw, 1. - eta, pBits);// encVData[i] = ( 1. - eta ) * ctmpw
				//scheme.reScaleByAndEqual(encVData[i], pBits-5);

				scheme.multByConstAndEqual(encWData[i], eta, pBits);// encWData[i] = eta * encWData[i]
				//scheme.reScaleByAndEqual(encWData[i], pBits-5);

				if (encWData[i].logq > encVData[i].logq) scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
				if (encWData[i].logq < encVData[i].logq) scheme.modDownToAndEqual(encVData[i], encWData[i].logq);
				if (encWData[i].logp != encVData[i].logp) {cout << "logp != logp"; exit(0);}

				scheme.addAndEqual(encVData[i], encWData[i]); // encVData[i] = encVData[i] + encWData[i]
				// encVData[i] = ( 1. - eta ) * ctmpw + eta * encWData[i]

				scheme.reScaleByAndEqual(encVData[i], pBits);
				encWData[i].copy(ctmpw);

				ctmpw.kill();
			}
			NTL_EXEC_RANGE_END

			;
			/* CipherGD::encNLGDstep(encWData, encVData, encGrad, eta, cnum, pBits); */

auto end_timeUpdatingWeights = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end_timeUpdatingWeights - start_timeUpdatingWeights;
        totalUpdatingWeightsTime += duration.count();  // Returns the elapsed time in seconds
        totalUpdatingWeights += 1;


			for (long i = 0; i < cnum; ++i)
				encGrad[i].kill();
			delete[] encGrad;
			/* cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits); */

			timeutils.stop("NesterovWithG : " + to_string(iter + 1)  + " -th iteration");

			openFileTIME << "," << timeutils.timeElapsed;
			openFileTIME.flush();
			openFileTIMELabel << ","
			                  << "NesterovWithGminBatch : " + to_string(iter + 1)
			                  + " -th epoch." + to_string(r + 1)
			                  + "-th iteration";
			openFileTIMELabel.flush();
			openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
			openFileCurrMEM.flush();
			openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
			openFilePeakMEM.flush();

			for (long i = 0; i < cnum; ++i)
				encZData[i].kill();
			delete[] encZData;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			cout << endl << "---------- TEST : THE " << iter * rnum + r + 1
			     << "-th ITERATION : Weights, AUC, MLE ----------" << endl;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);     */
			for (long i = 0; i < (cnum - 1); ++i) {
				complex<double> *dcvv = scheme.decrypt(secretKey, encVData[i]);
				for (long j = 0; j < batch; ++j) {
					cvData[batch * i + j] = dcvv[j].real();
				}
				delete[] dcvv;
			}
			complex<double> *dcvv = scheme.decrypt(secretKey,
			                                       encVData[cnum - 1]);
			long rest = factorDim - batch * (cnum - 1);
			for (long j = 0; j < rest; ++j) {
				cvData[batch * (cnum - 1) + j] = dcvv[j].real();
			}
			delete[] dcvv;
			/* cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits); */
			cout << "Current cWdata (encVData) : " << endl;
			for (long i = 0; i < factorDim; ++i)
				cout << setiosflags(ios::fixed) << setprecision(12) << cvData[i]
				     << ",\t";
			cout << endl;

			openFileTestAUC << ","
			                << MyTools::calculateAUC(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestAUC.flush();
			openFileTrainAUC << ","
			                 << MyTools::calculateAUC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainAUC.flush();

			openFileTestACC << ","
			                << MyTools::calculateACC(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestACC.flush();
			openFileTrainACC << ","
			                 << MyTools::calculateACC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainACC.flush();

			openFileTestMLE << ","
			                << MyTools::calculateMLE(zDataTest, cvData, factorDim,
			                        testSampleDim, enccor, encauc);
			openFileTestMLE.flush();
			openFileTrainMLE << ","
			                 << MyTools::calculateMLE(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc);
			openFileTrainMLE.flush();

			cout << endl << endl << endl << endl << endl;

			cout << "TrainMLE : "
			     << MyTools::calculateMLE(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << "TrainAUC : "
			     << MyTools::calculateAUC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << "TrainACC : "
			     << MyTools::calculateACC(zDataTrain, cvData, factorDim,
			                         trainSampleDim, enccor, encauc) << endl;
			cout << endl;
			cout << "Test MLE : "
			     << MyTools::calculateMLE(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << "Test AUC : "
			     << MyTools::calculateAUC(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << "Test ACC : "
			     << MyTools::calculateACC(zDataTest, cvData, factorDim,
			                              testSampleDim, enccor, encauc) << endl;
			cout << " ----------------- end of the " << (iter + 1)
			     << "-th epoch ----------------- " << endl;
			cout << " ----------------- end of the " << (r + 1)
			     << "-th iteration ----------------- " << endl << endl;
			cout
			        << "--------------------------------------------------------------------------------"
			        << endl;
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			cout << endl << endl << endl << endl << "encVData[0].logq = "
			     << encVData[0].logq << endl << endl << endl << endl;

			if (encVData[0].logq < 220 + encVData[0].logp) {
				//if ( encVData[0].logq <= 450 && iter < Epoch_Number-1 || encVData[0].logq < wBits && iter == Epoch_Number-1) {

				timeutils.start("Use Bootstrap To Recrypt Ciphertext");
				auto start_timeRefreshingWeights = std::chrono::high_resolution_clock::now();

			cout << endl << "before Bootstrap: encWData[0].logq = " << encWData[0].logq << endl;
				NTL_EXEC_RANGE(cnum, first, last);
				for (long i = first; i < last; ++i) {
					scheme.modDownToAndEqual(encWData[i], bootlogq);
					encWData[i].n = batch;
					scheme.bootstrapAndEqual(encWData[i], bootlogq, 990, logT, logI);
					encWData[i].n = slots;
				}
				NTL_EXEC_RANGE_END
				NTL_EXEC_RANGE(cnum, first, last);
				for (long i = first; i < last; ++i) {
					scheme.modDownToAndEqual(encVData[i], bootlogq);
					encVData[i].n = batch;
					scheme.bootstrapAndEqual(encVData[i], bootlogq, 990, logT, logI);
					encVData[i].n = slots;
				}
				NTL_EXEC_RANGE_END
			cout << endl << "after  Bootstrap: encWData[0].logq = " << encWData[0].logq << endl;
				auto end_timeRefreshingWeights = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duretion = end_timeRefreshingWeights - start_timeRefreshingWeights;
        totalRefreshingWeightsTime += duretion.count();  // Returns the elapsed time in seconds
        totalBootstrappingOperations += 1;



				timeutils.stop("Use Bootstrap To Recrypt Ciphertext");

				openFileTIME << "," << timeutils.timeElapsed;
				openFileTIME.flush();
				openFileTIMELabel << "," << "Bootstrapping";
				openFileTIMELabel.flush();
				openFileCurrMEM << "," << (MyTools::getCurrentRSS() >> 20);
				openFileCurrMEM.flush();
				openFilePeakMEM << "," << (MyTools::getPeakRSS() >> 20);
				openFilePeakMEM.flush();

			}
			/////////////////////////////////////////////////////////////////////////////
			//        BOOTSTRAPPING                                                    //
			//             Over and Out                                                //
			/////////////////////////////////////////////////////////////////////////////
			cout
			        << "--------------------------------------------------------------------------------"
			        << endl << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << endl << " !!! STOP " << iter + 1 << " ITERATION !!! "
			     << endl << endl;
			     			     	cout << endl << endl << endl << endl;
			     	cout << "totalUpdatingWeightsTime = " << totalUpdatingWeightsTime << endl;
			     	cout << "totalUpdatingWeights = " << totalUpdatingWeights << endl;
				cout << "totalRefreshingWeightsTime = " << totalRefreshingWeightsTime << endl;
				cout << "totalBootstrappingOperations = " << totalBootstrappingOperations << endl;
				cout << endl << endl << endl << endl;
		}
	}

	openFileTIME << endl;
	openFileTIME.flush();
	openFileTIMELabel << endl;
	openFileTIMELabel.flush();
	openFileTestAUC << endl;
	openFileTestAUC.flush();
	openFileTrainAUC << endl;
	openFileTrainAUC.flush();
	openFileTestACC << endl;
	openFileTestACC.flush();
	openFileTrainACC << endl;
	openFileTrainACC.flush();
	openFileTestMLE << endl;
	openFileTestMLE.flush();
	openFileTrainMLE << endl;
	openFileTrainMLE.flush();
	openFileCurrMEM << endl;
	openFileCurrMEM.flush();
	openFilePeakMEM << endl;
	openFilePeakMEM.flush();

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
