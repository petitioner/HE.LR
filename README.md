# DONT Delete This Repository or Change its NAME !!!
    It is the public source code of the first version of the HE.LR work including the mini-batch version.


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
