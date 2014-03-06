cd libs

# Install bowtie2
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.1/bowtie2-2.2.1-source.zip
unzip bowtie2-2.2.1-source.zip
rm bowtie2-2.2.1-source.zip
cd bowtie2-2.2.1/
make
cd ..

# Install cd-hit-est
wget https://cdhit.googlecode.com/files/cd-hit-v4.6.1-2012-08-27.tgz
tar zxvf cd-hit-v4.6.1-2012-08-27.tgz
rm cd-hit-v4.6.1-2012-08-27.tgz
cd cd-hit-v4.6.1-2012-08-27
make openmp=yes
cd ..