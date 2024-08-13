# DONT Delete This Repository or Change its NAME !!!
    It is the public source code of the first version of the HE.LR work including the mini-batch version.


```shell
Welcome to Ubuntu 24.04 LTS (GNU/Linux 6.8.0-39-generic x86_64)

 * Documentation:  https://help.ubuntu.com
 * Management:     https://landscape.canonical.com
 * Support:        https://ubuntu.com/pro

 System information as of Tue Aug 13 01:57:25 PM UTC 2024

  System load:             0.0
  Usage of /:              8.1% of 140.07GB
  Memory usage:            1%
  Swap usage:              0%
  Processes:               183
  Users logged in:         0
  IPv4 address for enp1s0: 107.191.42.70
  IPv6 address for enp1s0: 2001:19f0:5:5001:5400:5ff:fe0e:960b


Expanded Security Maintenance for Applications is not enabled.

0 updates can be applied immediately.

Enable ESM Apps to receive additional future security updates.
See https://ubuntu.com/esm or run: sudo pro status


The list of available updates is more than a week old.
To check for new updates run: sudo apt update


The programs included with the Ubuntu system are free software;
the exact distribution terms for each program are described in the
individual files in /usr/share/doc/*/copyright.

Ubuntu comes with ABSOLUTELY NO WARRANTY, to the extent permitted by
applicable law.

root@vultr:~# git clone https://github.com/petitioner/HE.LR
Cloning into 'HE.LR'...
remote: Enumerating objects: 1209, done.
remote: Counting objects: 100% (892/892), done.
remote: Compressing objects: 100% (651/651), done.
remote: Total 1209 (delta 312), reused 747 (delta 232), pack-reused 317 (from 1)
Receiving objects: 100% (1209/1209), 1.50 GiB | 110.52 MiB/s, done.
Resolving deltas: 100% (380/380), done.
Updating files: 100% (625/625), done.
root@vultr:~# ls
HE.LR
root@vultr:~# cd HE.LR/IDASH2017/lib/lib/
root@vultr:~/HE.LR/IDASH2017/lib/lib# unzip libntl.a.zip 
Archive:  libntl.a.zip
  inflating: libntl.a                
root@vultr:~/HE.LR/IDASH2017/lib/lib# cd ..
root@vultr:~/HE.LR/IDASH2017/lib# cd ..
root@vultr:~/HE.LR/IDASH2017# cd Default/
root@vultr:~/HE.LR/IDASH2017/Default# make clean
rm      ./src/MyIDASH2017.o ./src/MyMethods.o ./src/MyTools.o  ./HEAAN/HEAAN/src/BootContext.o ./HEAAN/HEAAN/src/Ciphertext.o ./HEAAN/HEAAN/src/EvaluatorUtils.o ./HEAAN/HEAAN/src/HEAAN.o ./HEAAN/HEAAN/src/Key.o ./HEAAN/HEAAN/src/Plaintext.o ./HEAAN/HEAAN/src/Ring.o ./HEAAN/HEAAN/src/RingMultiplier.o ./HEAAN/HEAAN/src/Scheme.o ./HEAAN/HEAAN/src/SchemeAlgo.o ./HEAAN/HEAAN/src/SecretKey.o ./HEAAN/HEAAN/src/SerializationUtils.o ./HEAAN/HEAAN/src/StringUtils.o ./HEAAN/HEAAN/src/TestScheme.o ./HEAAN/HEAAN/src/TimeUtils.o  ./src/MyIDASH2017.d ./src/MyMethods.d ./src/MyTools.d  ./HEAAN/HEAAN/src/BootContext.d ./HEAAN/HEAAN/src/Ciphertext.d ./HEAAN/HEAAN/src/EvaluatorUtils.d ./HEAAN/HEAAN/src/HEAAN.d ./HEAAN/HEAAN/src/Key.d ./HEAAN/HEAAN/src/Plaintext.d ./HEAAN/HEAAN/src/Ring.d ./HEAAN/HEAAN/src/RingMultiplier.d ./HEAAN/HEAAN/src/Scheme.d ./HEAAN/HEAAN/src/SchemeAlgo.d ./HEAAN/HEAAN/src/SecretKey.d ./HEAAN/HEAAN/src/SerializationUtils.d ./HEAAN/HEAAN/src/StringUtils.d ./HEAAN/HEAAN/src/TestScheme.d ./HEAAN/HEAAN/src/TimeUtils.d   CNNinference
rm: cannot remove 'CNNinference': No such file or directory
make: [makefile:75: clean] Error 1 (ignored)
 
root@vultr:~/HE.LR/IDASH2017/Default# make all
Building file: ../src/MyIDASH2017.cpp
Invoking: GCC C++ Compiler [CNNinference/Default/src/subdir.mk]
g++ -I/root/HE.LR/IDASH2017/lib/include -I/root/HE.LR/IDASH2017/HEAAN/HEAAN/src -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"src/MyIDASH2017.d" -MT"src/MyIDASH2017.o" -o "src/MyIDASH2017.o" "../src/MyIDASH2017.cpp"
../src/MyIDASH2017.cpp: In function ‘int main(int, char**)’:
../src/MyIDASH2017.cpp:49:18: warning: unused variable ‘zData’ [-Wunused-variable]
   49 |         double **zData = MyTools::dataFromFile(trainfile, trainfactorDim, trainSampleDim, traindataset, traindatalabel);
      |                  ^~~~~
../src/MyIDASH2017.cpp:50:18: warning: unused variable ‘zDate’ [-Wunused-variable]
   50 |         double **zDate = MyTools::dataFromFile(testfile, testfactorDim,         testSampleDim, testdataset, testdatalabel);
      |                  ^~~~~
../src/MyIDASH2017.cpp:86:14: warning: unused variable ‘testSampleNum’ [-Wunused-variable]
   86 |         long testSampleNum = testSampleDim;
      |              ^~~~~~~~~~~~~
In file included from ../src/MyMethods.h:4,
                 from ../src/MyIDASH2017.cpp:9:
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h: At global scope:
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:28:13: warning: ‘CONJUGATION’ defined but not used [-Wunused-variable]
   28 | static long CONJUGATION = 2;
      |             ^~~~~~~~~~~
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:27:13: warning: ‘MULTIPLICATION’ defined but not used [-Wunused-variable]
   27 | static long MULTIPLICATION  = 1;
      |             ^~~~~~~~~~~~~~
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:26:13: warning: ‘ENCRYPTION’ defined but not used [-Wunused-variable]
   26 | static long ENCRYPTION = 0;
      |             ^~~~~~~~~~
In file included from ../src/MyMethods.h:3:
../src/MyTools.h:81:15: warning: ‘degree7_36’ defined but not used [-Wunused-variable]
   81 | static double degree7_36[5] = {+0.5, -0.11634, +0.0011850, -0.0000062074, +0.000000011389};
      |               ^~~~~~~~~~
../src/MyTools.h:77:15: warning: ‘degree7_24’ defined but not used [-Wunused-variable]
   77 | static double degree7_24[5] = {+0.5, -0.15512, +0.0028088, -0.000026158,  +0.000000085324};
      |               ^~~~~~~~~~
../src/MyTools.h:74:15: warning: ‘degree7_12’ defined but not used [-Wunused-variable]
   74 | static double degree7_12[5] = {+0.5, -2.3514e-01, +1.2329e-02, -3.9345e-04, +4.7292e-06};  // 1 - poly(-yWTx)
      |               ^~~~~~~~~~
../src/MyTools.h:64:15: warning: ‘degree7’ defined but not used [-Wunused-variable]
   64 | static double degree7[5] = {+0.5, -0.10843500, +0.00102394287, -0.00000518229, +0.00000000934};  //|wTx| < 16
      |               ^~~~~~~
../src/MyTools.h:58:15: warning: ‘degree5’ defined but not used [-Wunused-variable]
   58 | static double degree5[4] = {+0.5, -0.19131, +0.0045963,   -0.0000412332};
      |               ^~~~~~~
../src/MyTools.h:50:15: warning: ‘degree3’ defined but not used [-Wunused-variable]
   50 | static double degree3[3] = {0.5,  -0.0843,          0.0002}; // MNIST : [-16, +16]
      |               ^~~~~~~
Finished building: ../src/MyIDASH2017.cpp
 
Building file: ../src/MyMethods.cpp
Invoking: GCC C++ Compiler [CNNinference/Default/src/subdir.mk]
g++ -I/root/HE.LR/IDASH2017/lib/include -I/root/HE.LR/IDASH2017/HEAAN/HEAAN/src -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"src/MyMethods.d" -MT"src/MyMethods.o" -o "src/MyMethods.o" "../src/MyMethods.cpp"
../src/MyMethods.cpp: In static member function ‘static double* MyMethods::testCryptoFullBatchNAGwithG(double**, double*, long int, long int, long int, double**, double*, long int, std::string)’:
../src/MyMethods.cpp:40:14: warning: unused variable ‘sdimBits’ [-Wunused-variable]
   40 |         long sdimBits = (long) ceil(log2(trainSampleDim));
      |              ^~~~~~~~
../src/MyMethods.cpp:199:17: warning: unused variable ‘cwData’ [-Wunused-variable]
  199 |         double *cwData = new double[factorDim]();
      |                 ^~~~~~
../src/MyMethods.cpp:202:21: warning: unused variable ‘encTrainData’ [-Wunused-variable]
  202 |         Ciphertext *encTrainData = new Ciphertext[rnum * cnum];
      |                     ^~~~~~~~~~~~
../src/MyMethods.cpp:203:21: warning: unused variable ‘encTrainLabel’ [-Wunused-variable]
  203 |         Ciphertext *encTrainLabel = new Ciphertext[rnum * cnum];
      |                     ^~~~~~~~~~~~~
../src/MyMethods.cpp:441:37: warning: unused variable ‘gamma’ [-Wunused-variable]
  441 |         double alpha0, alpha1, eta, gamma;
      |                                     ^~~~~
../src/MyMethods.cpp:442:32: warning: unused variable ‘truecor’ [-Wunused-variable]
  442 |         double enccor, encauc, truecor, trueauc;
      |                                ^~~~~~~
../src/MyMethods.cpp:442:41: warning: unused variable ‘trueauc’ [-Wunused-variable]
  442 |         double enccor, encauc, truecor, trueauc;
      |                                         ^~~~~~~
../src/MyMethods.cpp: In static member function ‘static double* MyMethods::testCryptoMiniBatchNAGwithG(double**, double*, long int, long int, long int, double**, double*, long int, std::string)’:
../src/MyMethods.cpp:953:14: warning: unused variable ‘sdimBits’ [-Wunused-variable]
  953 |         long sdimBits = (long) ceil(log2(trainSampleDim));
      |              ^~~~~~~~
../src/MyMethods.cpp:1112:17: warning: unused variable ‘cwData’ [-Wunused-variable]
 1112 |         double *cwData = new double[factorDim]();
      |                 ^~~~~~
../src/MyMethods.cpp:1115:21: warning: unused variable ‘encTrainData’ [-Wunused-variable]
 1115 |         Ciphertext *encTrainData = new Ciphertext[rnum * cnum];
      |                     ^~~~~~~~~~~~
../src/MyMethods.cpp:1116:21: warning: unused variable ‘encTrainLabel’ [-Wunused-variable]
 1116 |         Ciphertext *encTrainLabel = new Ciphertext[rnum * cnum];
      |                     ^~~~~~~~~~~~~
../src/MyMethods.cpp:1434:37: warning: unused variable ‘gamma’ [-Wunused-variable]
 1434 |         double alpha0, alpha1, eta, gamma;
      |                                     ^~~~~
../src/MyMethods.cpp:1435:32: warning: unused variable ‘truecor’ [-Wunused-variable]
 1435 |         double enccor, encauc, truecor, trueauc;
      |                                ^~~~~~~
../src/MyMethods.cpp:1435:41: warning: unused variable ‘trueauc’ [-Wunused-variable]
 1435 |         double enccor, encauc, truecor, trueauc;
      |                                         ^~~~~~~
../src/MyMethods.cpp: In static member function ‘static double* MyMethods::testCryptoFullBatchNAGwithG(double**, double*, long int, long int, long int, double**, double*, long int, std::string)’:
../src/MyMethods.cpp:934:1: warning: control reaches end of non-void function [-Wreturn-type]
  934 | }
      | ^
../src/MyMethods.cpp: In static member function ‘static double* MyMethods::testCryptoMiniBatchNAGwithG(double**, double*, long int, long int, long int, double**, double*, long int, std::string)’:
../src/MyMethods.cpp:1879:1: warning: control reaches end of non-void function [-Wreturn-type]
 1879 | }
      | ^
In file included from ../src/MyMethods.h:4,
                 from ../src/MyMethods.cpp:1:
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h: At global scope:
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:28:13: warning: ‘CONJUGATION’ defined but not used [-Wunused-variable]
   28 | static long CONJUGATION = 2;
      |             ^~~~~~~~~~~
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:27:13: warning: ‘MULTIPLICATION’ defined but not used [-Wunused-variable]
   27 | static long MULTIPLICATION  = 1;
      |             ^~~~~~~~~~~~~~
/root/HE.LR/IDASH2017/HEAAN/HEAAN/src/Scheme.h:26:13: warning: ‘ENCRYPTION’ defined but not used [-Wunused-variable]
   26 | static long ENCRYPTION = 0;
      |             ^~~~~~~~~~
In file included from ../src/MyMethods.h:3:
../src/MyTools.h:81:15: warning: ‘degree7_36’ defined but not used [-Wunused-variable]
   81 | static double degree7_36[5] = {+0.5, -0.11634, +0.0011850, -0.0000062074, +0.000000011389};
      |               ^~~~~~~~~~
../src/MyTools.h:77:15: warning: ‘degree7_24’ defined but not used [-Wunused-variable]
   77 | static double degree7_24[5] = {+0.5, -0.15512, +0.0028088, -0.000026158,  +0.000000085324};
      |               ^~~~~~~~~~
../src/MyTools.h:74:15: warning: ‘degree7_12’ defined but not used [-Wunused-variable]
   74 | static double degree7_12[5] = {+0.5, -2.3514e-01, +1.2329e-02, -3.9345e-04, +4.7292e-06};  // 1 - poly(-yWTx)
      |               ^~~~~~~~~~
../src/MyTools.h:64:15: warning: ‘degree7’ defined but not used [-Wunused-variable]
   64 | static double degree7[5] = {+0.5, -0.10843500, +0.00102394287, -0.00000518229, +0.00000000934};  //|wTx| < 16
      |               ^~~~~~~
../src/MyTools.h:58:15: warning: ‘degree5’ defined but not used [-Wunused-variable]
   58 | static double degree5[4] = {+0.5, -0.19131, +0.0045963,   -0.0000412332};
      |               ^~~~~~~
Finished building: ../src/MyMethods.cpp
 
Building file: ../src/MyTools.cpp
Invoking: GCC C++ Compiler [CNNinference/Default/src/subdir.mk]
g++ -I/root/HE.LR/IDASH2017/lib/include -I/root/HE.LR/IDASH2017/HEAAN/HEAAN/src -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"src/MyTools.d" -MT"src/MyTools.o" -o "src/MyTools.o" "../src/MyTools.cpp"
../src/MyTools.cpp: In static member function ‘static double** MyTools::zDataFromFile(std::string&, long int&, long int&, bool)’:
../src/MyTools.cpp:47:52: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::__cxx11::basic_string<char>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
   47 |                                 for (long i = 0; i < line.length(); ++i)
      |                                                  ~~^~~~~~~~~~~~~~~
../src/MyTools.cpp: In static member function ‘static double** MyTools::dataFromFile(std::string&, long int&, long int&, double**&, double*&)’:
../src/MyTools.cpp:118:52: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::__cxx11::basic_string<char>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
  118 |                                 for (long i = 0; i < line.length(); ++i)
      |                                                  ~~^~~~~~~~~~~~~~~
../src/MyTools.cpp: In static member function ‘static double MyTools::calculateAUC(double**, double*, long int, long int, double&, double&)’:
../src/MyTools.cpp:382:27: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::vector<double>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
  382 |         for(long i = 0; i < thetaTN.size(); ++i){
      |                         ~~^~~~~~~~~~~~~~~~
../src/MyTools.cpp:383:31: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::vector<double>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
  383 |             for(long j = 0; j < thetaFP.size(); ++j){
      |                             ~~^~~~~~~~~~~~~~~~
../src/MyTools.cpp: In static member function ‘static double MyTools::calculateACC(double**, double*, long int, long int, double&, double&)’:
../src/MyTools.cpp:423:27: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::vector<double>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
  423 |         for(long i = 0; i < thetaTN.size(); ++i){
      |                         ~~^~~~~~~~~~~~~~~~
../src/MyTools.cpp:424:31: warning: comparison of integer expressions of different signedness: ‘long int’ and ‘std::vector<double>::size_type’ {aka ‘long unsigned int’} [-Wsign-compare]
  424 |             for(long j = 0; j < thetaFP.size(); ++j){
      |                             ~~^~~~~~~~~~~~~~~~
In file included from ../src/MyTools.cpp:1:
../src/MyTools.h: At global scope:
../src/MyTools.h:81:15: warning: ‘degree7_36’ defined but not used [-Wunused-variable]
   81 | static double degree7_36[5] = {+0.5, -0.11634, +0.0011850, -0.0000062074, +0.000000011389};
      |               ^~~~~~~~~~
../src/MyTools.h:77:15: warning: ‘degree7_24’ defined but not used [-Wunused-variable]
   77 | static double degree7_24[5] = {+0.5, -0.15512, +0.0028088, -0.000026158,  +0.000000085324};
      |               ^~~~~~~~~~
../src/MyTools.h:74:15: warning: ‘degree7_12’ defined but not used [-Wunused-variable]
   74 | static double degree7_12[5] = {+0.5, -2.3514e-01, +1.2329e-02, -3.9345e-04, +4.7292e-06};  // 1 - poly(-yWTx)
      |               ^~~~~~~~~~
../src/MyTools.h:64:15: warning: ‘degree7’ defined but not used [-Wunused-variable]
   64 | static double degree7[5] = {+0.5, -0.10843500, +0.00102394287, -0.00000518229, +0.00000000934};  //|wTx| < 16
      |               ^~~~~~~
../src/MyTools.h:58:15: warning: ‘degree5’ defined but not used [-Wunused-variable]
   58 | static double degree5[4] = {+0.5, -0.19131, +0.0045963,   -0.0000412332};
      |               ^~~~~~~
../src/MyTools.h:50:15: warning: ‘degree3’ defined but not used [-Wunused-variable]
   50 | static double degree3[3] = {0.5,  -0.0843,          0.0002}; // MNIST : [-16, +16]
      |               ^~~~~~~
Finished building: ../src/MyTools.cpp
 
Building file: ../HEAAN/HEAAN/src/BootContext.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/BootContext.d" -MT"HEAAN/HEAAN/src/BootContext.o" -o "HEAAN/HEAAN/src/BootContext.o" "../HEAAN/HEAAN/src/BootContext.cpp"
Finished building: ../HEAAN/HEAAN/src/BootContext.cpp
 
Building file: ../HEAAN/HEAAN/src/Ciphertext.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/Ciphertext.d" -MT"HEAAN/HEAAN/src/Ciphertext.o" -o "HEAAN/HEAAN/src/Ciphertext.o" "../HEAAN/HEAAN/src/Ciphertext.cpp"
Finished building: ../HEAAN/HEAAN/src/Ciphertext.cpp
 
Building file: ../HEAAN/HEAAN/src/EvaluatorUtils.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/EvaluatorUtils.d" -MT"HEAAN/HEAAN/src/EvaluatorUtils.o" -o "HEAAN/HEAAN/src/EvaluatorUtils.o" "../HEAAN/HEAAN/src/EvaluatorUtils.cpp"
Finished building: ../HEAAN/HEAAN/src/EvaluatorUtils.cpp
 
Building file: ../HEAAN/HEAAN/src/HEAAN.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/HEAAN.d" -MT"HEAAN/HEAAN/src/HEAAN.o" -o "HEAAN/HEAAN/src/HEAAN.o" "../HEAAN/HEAAN/src/HEAAN.cpp"
Finished building: ../HEAAN/HEAAN/src/HEAAN.cpp
 
Building file: ../HEAAN/HEAAN/src/Key.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/Key.d" -MT"HEAAN/HEAAN/src/Key.o" -o "HEAAN/HEAAN/src/Key.o" "../HEAAN/HEAAN/src/Key.cpp"
Finished building: ../HEAAN/HEAAN/src/Key.cpp
 
Building file: ../HEAAN/HEAAN/src/Plaintext.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/Plaintext.d" -MT"HEAAN/HEAAN/src/Plaintext.o" -o "HEAAN/HEAAN/src/Plaintext.o" "../HEAAN/HEAAN/src/Plaintext.cpp"
Finished building: ../HEAAN/HEAAN/src/Plaintext.cpp
 
Building file: ../HEAAN/HEAAN/src/Ring.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/Ring.d" -MT"HEAAN/HEAAN/src/Ring.o" -o "HEAAN/HEAAN/src/Ring.o" "../HEAAN/HEAAN/src/Ring.cpp"
../HEAAN/HEAAN/src/Ring.cpp: In member function ‘void Ring::addBootContext(long int, long int)’:
../HEAAN/HEAAN/src/Ring.cpp:305:152: warning: ‘rp1’ may be used uninitialized [-Wmaybe-uninitialized]
  305 |          bootContextMap.insert(pair<long, BootContext*>(logSlots, new BootContext(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)));
      |                                                                                                                                                 ^

../HEAAN/HEAAN/src/Ring.cpp:167:27: note: ‘rp1’ was declared here
  167 |                 uint64_t* rp1;
      |                           ^~~
../HEAAN/HEAAN/src/Ring.cpp:305:152: warning: ‘rp2’ may be used uninitialized [-Wmaybe-uninitialized]
  305 |          bootContextMap.insert(pair<long, BootContext*>(logSlots, new BootContext(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)));
      |                                                                                                                                                 ^

../HEAAN/HEAAN/src/Ring.cpp:168:27: note: ‘rp2’ was declared here
  168 |                 uint64_t* rp2;
      |                           ^~~
../HEAAN/HEAAN/src/Ring.cpp:305:152: warning: ‘bnd1’ may be used uninitialized [-Wmaybe-uninitialized]
  305 |          bootContextMap.insert(pair<long, BootContext*>(logSlots, new BootContext(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)));
      |                                                                                                                                                 ^

../HEAAN/HEAAN/src/Ring.cpp:172:22: note: ‘bnd1’ was declared here
  172 |                 long bnd1;
      |                      ^~~~
../HEAAN/HEAAN/src/Ring.cpp:305:152: warning: ‘bnd2’ may be used uninitialized [-Wmaybe-uninitialized]
  305 |          bootContextMap.insert(pair<long, BootContext*>(logSlots, new BootContext(rpvec, rpvecInv, rp1, rp2, bndvec, bndvecInv, bnd1, bnd2, logp)));
      |                                                                                                                                                 ^

../HEAAN/HEAAN/src/Ring.cpp:173:22: note: ‘bnd2’ was declared here
  173 |                 long bnd2;
      |                      ^~~~
Finished building: ../HEAAN/HEAAN/src/Ring.cpp
 
Building file: ../HEAAN/HEAAN/src/RingMultiplier.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/RingMultiplier.d" -MT"HEAAN/HEAAN/src/RingMultiplier.o" -o "HEAAN/HEAAN/src/RingMultiplier.o" "../HEAAN/HEAAN/src/RingMultiplier.cpp"
../HEAAN/HEAAN/src/RingMultiplier.cpp: In member function ‘void RingMultiplier::multAndEqual(NTL::ZZ*, NTL::ZZ*, long int, const NTL::ZZ&)’:
../HEAAN/HEAAN/src/RingMultiplier.cpp:304:13: warning: unused variable ‘pHatnp’ [-Wunused-variable]
  304 |         ZZ* pHatnp = pHat[np - 1];
      |             ^~~~~~
../HEAAN/HEAAN/src/RingMultiplier.cpp:305:19: warning: unused variable ‘pHatInvModpnp’ [-Wunused-variable]
  305 |         uint64_t* pHatInvModpnp = pHatInvModp[np - 1];
      |                   ^~~~~~~~~~~~~
../HEAAN/HEAAN/src/RingMultiplier.cpp: In member function ‘void RingMultiplier::multNTTAndEqual(NTL::ZZ*, uint64_t*, long int, const NTL::ZZ&)’:
../HEAAN/HEAAN/src/RingMultiplier.cpp:334:13: warning: unused variable ‘pHatnp’ [-Wunused-variable]
  334 |         ZZ* pHatnp = pHat[np - 1];
      |             ^~~~~~
../HEAAN/HEAAN/src/RingMultiplier.cpp:335:19: warning: unused variable ‘pHatInvModpnp’ [-Wunused-variable]
  335 |         uint64_t* pHatInvModpnp = pHatInvModp[np - 1];
      |                   ^~~~~~~~~~~~~
../HEAAN/HEAAN/src/RingMultiplier.cpp: In member function ‘void RingMultiplier::square(NTL::ZZ*, NTL::ZZ*, long int, const NTL::ZZ&)’:
../HEAAN/HEAAN/src/RingMultiplier.cpp:365:13: warning: unused variable ‘pHatnp’ [-Wunused-variable]
  365 |         ZZ* pHatnp = pHat[np - 1];
      |             ^~~~~~
../HEAAN/HEAAN/src/RingMultiplier.cpp:366:19: warning: unused variable ‘pHatInvModpnp’ [-Wunused-variable]
  366 |         uint64_t* pHatInvModpnp = pHatInvModp[np - 1];
      |                   ^~~~~~~~~~~~~
Finished building: ../HEAAN/HEAAN/src/RingMultiplier.cpp
 
Building file: ../HEAAN/HEAAN/src/Scheme.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/Scheme.d" -MT"HEAAN/HEAAN/src/Scheme.o" -o "HEAAN/HEAAN/src/Scheme.o" "../HEAAN/HEAAN/src/Scheme.cpp"
Finished building: ../HEAAN/HEAAN/src/Scheme.cpp
 
Building file: ../HEAAN/HEAAN/src/SchemeAlgo.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/SchemeAlgo.d" -MT"HEAAN/HEAAN/src/SchemeAlgo.o" -o "HEAAN/HEAAN/src/SchemeAlgo.o" "../HEAAN/HEAAN/src/SchemeAlgo.cpp"
In file included from ../HEAAN/HEAAN/src/SchemeAlgo.h:18,
                 from ../HEAAN/HEAAN/src/SchemeAlgo.cpp:8:
../HEAAN/HEAAN/src/Scheme.h:28:13: warning: ‘CONJUGATION’ defined but not used [-Wunused-variable]
   28 | static long CONJUGATION = 2;
      |             ^~~~~~~~~~~
../HEAAN/HEAAN/src/Scheme.h:27:13: warning: ‘MULTIPLICATION’ defined but not used [-Wunused-variable]
   27 | static long MULTIPLICATION  = 1;
      |             ^~~~~~~~~~~~~~
../HEAAN/HEAAN/src/Scheme.h:26:13: warning: ‘ENCRYPTION’ defined but not used [-Wunused-variable]
   26 | static long ENCRYPTION = 0;
      |             ^~~~~~~~~~
Finished building: ../HEAAN/HEAAN/src/SchemeAlgo.cpp
 
Building file: ../HEAAN/HEAAN/src/SecretKey.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/SecretKey.d" -MT"HEAAN/HEAAN/src/SecretKey.o" -o "HEAAN/HEAAN/src/SecretKey.o" "../HEAAN/HEAAN/src/SecretKey.cpp"
Finished building: ../HEAAN/HEAAN/src/SecretKey.cpp
 
Building file: ../HEAAN/HEAAN/src/SerializationUtils.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/SerializationUtils.d" -MT"HEAAN/HEAAN/src/SerializationUtils.o" -o "HEAAN/HEAAN/src/SerializationUtils.o" "../HEAAN/HEAAN/src/SerializationUtils.cpp"
../HEAAN/HEAAN/src/SerializationUtils.cpp: In static member function ‘static Ciphertext* SerializationUtils::readCiphertext(std::string)’:
../HEAAN/HEAAN/src/SerializationUtils.cpp:56:16: warning: address of local variable ‘cipher’ returned [-Wreturn-local-addr]
   56 |         return &cipher;
      |                ^~~~~~~
../HEAAN/HEAAN/src/SerializationUtils.cpp:46:20: note: declared here
   46 |         Ciphertext cipher(logp, logq, n);
      |                    ^~~~~~
../HEAAN/HEAAN/src/SerializationUtils.cpp: In static member function ‘static Key* SerializationUtils::readKey(std::string)’:
../HEAAN/HEAAN/src/SerializationUtils.cpp:74:16: warning: address of local variable ‘key’ returned [-Wreturn-local-addr]
   74 |         return &key;
      |                ^~~~
../HEAAN/HEAAN/src/SerializationUtils.cpp:68:13: note: declared here
   68 |         Key key;
      |             ^~~
Finished building: ../HEAAN/HEAAN/src/SerializationUtils.cpp
 
Building file: ../HEAAN/HEAAN/src/StringUtils.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/StringUtils.d" -MT"HEAAN/HEAAN/src/StringUtils.o" -o "HEAAN/HEAAN/src/StringUtils.o" "../HEAAN/HEAAN/src/StringUtils.cpp"
Finished building: ../HEAAN/HEAAN/src/StringUtils.cpp
 
Building file: ../HEAAN/HEAAN/src/TestScheme.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/TestScheme.d" -MT"HEAAN/HEAAN/src/TestScheme.o" -o "HEAAN/HEAAN/src/TestScheme.o" "../HEAAN/HEAAN/src/TestScheme.cpp"
In file included from ../HEAAN/HEAAN/src/TestScheme.cpp:16:
../HEAAN/HEAAN/src/Scheme.h:28:13: warning: ‘CONJUGATION’ defined but not used [-Wunused-variable]
   28 | static long CONJUGATION = 2;
      |             ^~~~~~~~~~~
../HEAAN/HEAAN/src/Scheme.h:27:13: warning: ‘MULTIPLICATION’ defined but not used [-Wunused-variable]
   27 | static long MULTIPLICATION  = 1;
      |             ^~~~~~~~~~~~~~
../HEAAN/HEAAN/src/Scheme.h:26:13: warning: ‘ENCRYPTION’ defined but not used [-Wunused-variable]
   26 | static long ENCRYPTION = 0;
      |             ^~~~~~~~~~
Finished building: ../HEAAN/HEAAN/src/TestScheme.cpp
 
Building file: ../HEAAN/HEAAN/src/TimeUtils.cpp
Invoking: GCC C++ Compiler 4 Testing
g++ -I"/root/HE.LR/IDASH2017/HEAAN/HEAAN/src" -I"/root/HE.LR/IDASH2017/lib/include" -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"HEAAN/HEAAN/src/TimeUtils.d" -MT"HEAAN/HEAAN/src/TimeUtils.o" -o "HEAAN/HEAAN/src/TimeUtils.o" "../HEAAN/HEAAN/src/TimeUtils.cpp"
Finished building: ../HEAAN/HEAAN/src/TimeUtils.cpp
 
Building target: iDASH2017
Invoking: GCC C++ Linker4Testing
g++ -L/root/HE.LR/IDASH2017/lib/lib -L/root/HE.LR/IDASH2017/HEAAN/HEAAN/lib -pthread -o iDASH2017 ./src/MyIDASH2017.o ./src/MyMethods.o ./src/MyTools.o  ./HEAAN/HEAAN/src/BootContext.o ./HEAAN/HEAAN/src/Ciphertext.o ./HEAAN/HEAAN/src/EvaluatorUtils.o ./HEAAN/HEAAN/src/HEAAN.o ./HEAAN/HEAAN/src/Key.o ./HEAAN/HEAAN/src/Plaintext.o ./HEAAN/HEAAN/src/Ring.o ./HEAAN/HEAAN/src/RingMultiplier.o ./HEAAN/HEAAN/src/Scheme.o ./HEAAN/HEAAN/src/SchemeAlgo.o ./HEAAN/HEAAN/src/SecretKey.o ./HEAAN/HEAAN/src/SerializationUtils.o ./HEAAN/HEAAN/src/StringUtils.o ./HEAAN/HEAAN/src/TestScheme.o ./HEAAN/HEAAN/src/TimeUtils.o   -lntl -lgmp -lHEAAN -lm
/usr/bin/ld: /root/HE.LR/IDASH2017/lib/lib/libgmp.a(gcd_1.o): warning: relocation in read-only section `.text'
/usr/bin/ld: warning: creating DT_TEXTREL in a PIE
Finished building target: iDASH2017

```
# HE.LR
```shell

Ubuntu 2024.04  comes with ABSOLUTELY NO WARRANTY, to the extent permitted by
applicable law.

root@vultr:~# sudo swapoff /swapfile
root@vultr:~# sudo fallocate -l 16G /swapfile
root@vultr:~# sudo mkswap /swapfile
root@vultr:~# sudo swapon /swapfile
root@vultr:~# cd ..
root@vultr:/# ls
root@vultr:/# cd home/
root@vultr:/home# mkdir sunly
root@vultr:/home# cd sunly/
root@vultr:/home/sunly# mkdir Downloads
root@vultr:/home/sunly# cd Downloads/
root@vultr:/home/sunly/Downloads# mkdir IDASH2017-master
root@vultr:/home/sunly/Downloads# cd IDASH2017-master/
root@vultr:/home/sunly/Downloads/IDASH2017-master# ls
IDASH2017.zip
root@vultr:/home/sunly/Downloads/IDASH2017-master# unzip IDASH2017.zip 

root@vultr:/home/sunly/Downloads/IDASH2017-master# cd IDASH2017/
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017# cd lib/lib/
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# ls
README.md  libgmp.a  libgmpxx.a  libntl.a.zip
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# unzip libntl.a.zip 
Archive:  libntl.a.zip
  inflating: libntl.a                
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# cd ..
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib# cd ..
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017# cd Default/
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/Default# make clean 
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/Default# make all


66.135.23.150

B=g2Fn_n@Nxa}ME!

sftp://root@66.135.23.150
```


# How to run it:
```cpp        
Step 0: Clone the GitHub Repository
Execute the following command in the command line to clone the GitHub repository named HE.LR:
git clone https://github.com/petitioner/HE.LR

Step 1: Rename the Repository
Use the mv command to rename the cloned repository to IDASH2017-master:
mv HE.LR IDASH2017-master

Step 2: Place IDASH2017-master in the Correct Directory
Place the IDASH2017-master repository in the specified directory to align with /home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib.

Step 3: Unzip libntl.a.zip to IDASH2017/lib/lib/libntl.a




root@vultr:/# cd home/
root@vultr:/home# mkdir sunly
root@vultr:/home# cd sunly/
root@vultr:/home/sunly# mkdir Downloads
root@vultr:/home/sunly# cd Downloads/
root@vultr:/home/sunly/Downloads# git clone https://github.com/petitioner/HE.LR
root@vultr:/home/sunly/Downloads# mv HE.LR IDASH2017-master
root@vultr:/home/sunly/Downloads# cd IDASH2017-master/IDASH2017/lib/lib/
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# ls
libgmp.a  libgmpxx.a  libntl.a.zip  README.md
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# unzip libntl.a.zip 
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib/lib# cd ..
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/lib# cd ..
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017# cd Default/
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/Default# make clean 
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/Default# make all

Finished building target: iDASH2017
 
root@vultr:/home/sunly/Downloads/IDASH2017-master/IDASH2017/Default# ./iDASH2017 

```


---

# Privacy-Preserving Logistic Regression Training on Large Datasets 
        [a follow-up work for the first paper but on the mini-batch version]

---
## https://github.com/K-miran/HELR
## https://github.com/kimandrik/IDASH2017
## https://github.com/KyoohyungHan/HELR
