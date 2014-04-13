#! /bin/sh

# Initialize scripts
. /afs/cern.ch/sw/lcg/external/gcc/4.8.1/x86_64-slc6/setup.sh
. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

# Compiling program
echo "----------------------------------------"
echo "------------- COMPILING ----------------"
echo "----------------------------------------"

cd /afs/cern.ch/work/j/jeyserma/private/GEMSimulations/gain/
make

# Executing rules...
echo "----------------------------------------"
echo "-------------- EXECUTE -----------------"
echo "----------------------------------------"

#./gain 2_1
./gain 2_2
./gain 2_3
./gain 2_4
./gain 2_5
./gain 2_6
./gain 2_7
./gain 2_8
./gain 2_9
#./gain 2_10
