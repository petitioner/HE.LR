/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "MyTools.h"

#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>    // std::shuffle + std::copy
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

using namespace NTL;

long MyTools::suggestLogN(long lambda, long logQ) {
	long NBnd = ceil(logQ * (lambda + 110) / 3.6);
	double logNBnd = log2((double)NBnd);
	return (long)ceil(logNBnd);
}

bool MyTools::is_number(string& s)
{
	return( strspn( s.c_str(), "-.0123456789" ) == s.size() );
}

/**
 * Load data from train file, extract data to X and Y, transform Y{-0,+1} to Y{-1,+1};
 * @param  : factorDim : dimension of the raw data [X:y]
 * @param  : sampleDim : the number of rows in the data, except for the first row which is a interpretation.
 * @param  : isfirst : is Y at the beginning of each row ?
 * @return :
 * @author : no one
 */
double** MyTools::zDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst) {
	vector<vector<double>> zline;
	factorDim = 1; 	// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		
		while (getline(openFile, line)) {
			if(factorDim == 1)
				for (long i = 0; i < line.length(); ++i)	
					if (line[i] == ',')	factorDim++;

			size_t pos = line.find_first_of(',', 0);
			string first_word = line.substr(0, pos);
			if(!is_number(first_word))
				getline(openFile, line);

			size_t start, end;	
			vector<double> vecline;
			do {
				end = line.find_first_of(',', start);
				temp = line.substr(start, end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while (start);
			zline.push_back(vecline);
			sampleDim++;
		}

	} else {
		cout << "Error: cannot read file" << endl;
	}

	double** zData = new double*[sampleDim];
	if(isfirst) {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][0] - 1;               // change class label Y{0,1} to Y{-1,+1}
			for(long i = 1; i < factorDim; ++i){       // change data entry Y=-1       X{0,-1}
				zj[i] = zj[0] * zline[j][i];           // change data entry Y= 1       X{0, 1}
			}
			zData[j] = zj;
		}
	} else {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][factorDim - 1] - 1;
			for(long i = 1; i < factorDim; ++i){
				zj[i] = zj[0] * zline[j][i-1];
			}
			zData[j] = zj;
		}
	}
	return zData;
}
/**
 * To Extract X and Y from data.
 *      Step 2. Extract X and Y from data
 *              get X and Y as follows:
 *              X = | 1 X11 X12 X13 ... X1d|
 *                  | 1 X21 X22 X23 ... X2d|
 *                  | .  .   .   .  ...  . |
 *                  | .  .   .   .  ...  . |
 *                  | 1 Xn1 Xn2 Xn3 ... Xnd|
 *              Y = [ Y1 Y2  Y3  Y4 ... Yn ]
 *              turn y{+0,+1} to y{-1,+1}
 * @param  : sampleDim : the number of rows in the data
 * @return : zData = (Y{-1,+1}, Y@X)
 * @author : no one
 */
double** MyTools::dataFromFile(string& path, long& factorDim, long& sampleDim, double** &X, double* &Y){
	vector<vector<double>> zline;
	factorDim = 1; 	// dimension of (Y,X)
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if (openFile.is_open()) {
		string line, temp;
		
		while (getline(openFile, line)) {
			if(factorDim == 1)
				for (long i = 0; i < line.length(); ++i)	
					if (line[i] == ',')	factorDim++;

			size_t pos = line.find_first_of(',', 0);
			string first_word = line.substr(0, pos);
			if(!is_number(first_word))
				getline(openFile, line);

			size_t start, end;	
			vector<double> vecline;
			do {
				end = line.find_first_of(',', start);
				temp = line.substr(start, end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while (start);
			zline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}

	double** zData = new double*[sampleDim];
	X = new double*[sampleDim];
	Y = new double[sampleDim];

	bool flagY01 = true;
	for(int i=0;i<sampleDim;++i){
		if(zline[i][0]==-1){
			flagY01 = false;
			break;
		}
	}

	for (long j = 0; j < sampleDim; ++j) {
		double* zj = new double[factorDim];
		double* xj = new double[factorDim];
		Y[j] = zline[j][0];
		if(flagY01){
			zj[0] = 2 * zline[j][0] - 1;    // change class label Y{0,1} to Y{-1,+1}
			Y[j] = 2 * Y[j] - 1;
		}
		xj[0] = 1.0;
		for (long i = 1; i < factorDim; ++i) { // change data entry Y=-1       X{0,-1}
			zj[i] = zj[0] * zline[j][i]; // change data entry Y= 1       X{0, 1}

			xj[i] = zline[j][i];
		}
		zData[j] = zj;
		X[j] = xj;
	}

	return zData;
}
/**
 * Shuffle the order of each row in zData.
 * @param  : zData 
 * @param  : factorDim : dimension of the raw data [X:y]
 * @param  : sampleDim : the number of rows in the data, except for the first row which is a interpretation.
 * @param  : isfirst :
 * @return :
 * @author : no one
 */
void MyTools::shuffleZData(double** zData, long factorDim, long sampleDim) {
	srand(time(NULL));
	double* tmp = new double[factorDim];
	for (long i = 0; i < sampleDim; ++i) {
		long idx = i + rand() / (RAND_MAX / (sampleDim - i) + 1);
		copy(zData[i], zData[i] + factorDim, tmp);
		copy(zData[idx], zData[idx] + factorDim, zData[i]);
		copy(tmp, tmp + factorDim, zData[idx]);
	}
}
/**
 * Shuffle the order of each row in X and each element in Y in the same time.
 *
 * @param  : factorDim : dimension of the raw data X
 * @param  : sampleDim : the number of rows in the data.
 * @param  : isfirst :
 * @return :
 *
 * @author : no one
 */
void MyTools::shuffleDataSync(double** X, long factorDim, long sampleDim, double* Y){
	srand(time(NULL));
	double* tmp = new double[factorDim];
	double temp;
	for (long i = 0; i < sampleDim; ++i) {
		long idx = i + rand() / (RAND_MAX / (sampleDim - i) + 1);
		copy(X[i], X[i] + factorDim, tmp);      temp = Y[i];
		copy(X[idx], X[idx] + factorDim, X[i]); Y[i] = Y[idx];
		copy(tmp, tmp + factorDim, X[idx]);     Y[idx] = temp;
	}
}

/**
 * Shift the range of each element in the same column to the range [-1,+1];
 * In fact, it nomalize the data to the range [0, 1] since it combine y*X to X before.
 *
 * @param  : zData 
 * @param  : factorDim : dimension of the raw data [X:y]
 * @param  : sampleDim : the number of rows in the data.
 * @return :
 * @author : no one
 */
void MyTools::normalizeZData(double** zData, long factorDim, long sampleDim) {
	long i, j;
	double m;
	for (i = 0; i < factorDim; ++i) {
		m = 0.0;
		for (j = 0; j < sampleDim; ++j) {
			m = max(m, abs(zData[j][i]));
		}

		if(m < 1e-10) continue;

		for (j = 0; j < sampleDim; ++j) {
			zData[j][i] /= m;
		}
	}
}



double MyTools::trueIP(double* a, double* b, long size) {
    	double res = 0.0;
    	for(long i = 0; i < size; ++i) {
    		res += a[i] * b[i];
    	}
    	return res;
    }
double MyTools::calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
//	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
//		cout << wData[i] << ",";
	}
//	cout << endl;

	long TN = 0, FP = 0;

    vector<double> thetaTN;
    vector<double> thetaFP;

    for(int i = 0; i < sampleDim; ++i){
        if(zData[i][0] > 0){
            if(MyTools::trueIP(zData[i], wData, factorDim) < 0) TN++;
            thetaTN.push_back(zData[i][0] * MyTools::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(MyTools::trueIP(zData[i], wData, factorDim) < 0) FP++;
            thetaFP.push_back(zData[i][0] * MyTools::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        }
    }

    correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);

    if(thetaFP.size() == 0 || thetaTN.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        auc = 0.0;
    } else{
        auc = 0.0;
        for(long i = 0; i < thetaTN.size(); ++i){
            for(long j = 0; j < thetaFP.size(); ++j){
                if(thetaFP[j] <= thetaTN[i]) auc++;
            }
        }
        auc /= thetaTN.size() * thetaFP.size();
        cout << "AUC: " << auc << endl;
    }

    return auc;
}

double MyTools::calculateACC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
//	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
//		cout << wData[i] << ",";
	}
//	cout << endl;

	long TN = 0, FP = 0;

    vector<double> thetaTN;
    vector<double> thetaFP;

    for(int i = 0; i < sampleDim; ++i){
        if(zData[i][0] > 0){
            if(MyTools::trueIP(zData[i], wData, factorDim) < 0) TN++;
            thetaTN.push_back(zData[i][0] * MyTools::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(MyTools::trueIP(zData[i], wData, factorDim) < 0) FP++;
            thetaFP.push_back(zData[i][0] * MyTools::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        }
    }

    correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);

    if(thetaFP.size() == 0 || thetaTN.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        auc = 0.0;
    } else{
        auc = 0.0;
        for(long i = 0; i < thetaTN.size(); ++i){
            for(long j = 0; j < thetaFP.size(); ++j){
                if(thetaFP[j] <= thetaTN[i]) auc++;
            }
        }
        auc /= thetaTN.size() * thetaFP.size();
        //cout << "AUC: " << auc << endl;
    }

    return correctness;
}

double MyTools::calculateMLE(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {

	double* XW = new double[sampleDim]();
	for(int i=0;i<sampleDim;++i)
	{
		double res = 0;
		for(int j=0;j<factorDim;++j)
			res += zData[i][j]*wData[j];
		XW[i] = res;
	}
	double mle = 0.0;
	for(int i=0;i<sampleDim;++i)
	{
		double logsigm = -log(1+exp(-XW[i]));
		mle += logsigm;
	}

    return mle;
}



/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t MyTools::getPeakRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}


/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t MyTools::getCurrentRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}


