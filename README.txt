# To compile on Linux:
make all

# To run program (automatically distributed on multiple threads using OpemMP):
./bin/nndistance config_test.txt
./bin/nndecluster config_test.txt


# To control explicitely the number of threads, use:
OMP_NUM_THREADS=5 ./bin/nndistance config_test.txt
OMP_NUM_THREADS=5 ./bin/nndecluster config_test.txt


# To compile on Windows 11:
g++.exe  src\strutils.cc src\arrayutils.cc src\fileutils.cc src\nnanalysis.cc src\nndistance.cc src\nnstats.cc -O3 -fopenmp -I.\include -IC:\Users\b94678\Documents\Gratuiciels\boost_1_82_0\boost -o .\bin\nndistance.exe

g++.exe  src\strutils.cc src\arrayutils.cc src\fileutils.cc src\nnanalysis.cc src\nndeclust.cc src\nnstats.cc -O3 -fopenmp -I.\include -IC:\Users\b94678\Documents\Gratuiciels\boost_1_82_0\boost -o .\bin\nndecluster.exe
