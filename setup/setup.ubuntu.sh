#
# This shell script will install all the dependent software
# necessary to run the ASR Pipeline
#
wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py -O - | python
sudo easy_install -U dendropy
sudo apt-get install python-biopython

mkdir common

# MUSCLE
curl -O http://drive5.com/muscle/downloads3.6/muscle3.6_src.tar.gz
mv muscle3.6_src.tar.gz /common/
cd common
gzip -d muscle3.6_src.tar.gz
tar -xvf muscle3.6_src.tar
cd muscle3.6_src/
make
ln -s /common/muscle3.6_src/muscle  /bin/muscle
cd /

#MSAPROBS
cp ~/asr-pipeline/apps/MSAProbs-0.9.7.tar /common/
cd /common
tar -xvf MSAProbs-0.9.7.tar 
cd MSAProbs-0.9.7/MSAProbs/
make
ln -s /common/MSAProbs-0.9.7/MSAProbs/msaprobs /bin/msaprobs




