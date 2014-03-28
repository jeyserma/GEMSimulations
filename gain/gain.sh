#! /bin/sh

# Initialize scripts
. /afs/cern.ch/sw/lcg/external/gcc/4.8.1/x86_64-slc6/setup.sh
. /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.10/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

# Compiling program
echo "----------------------------------------"
echo "------------- COMPILING ----------------"
echo "----------------------------------------"

cd /afs/cern.ch/work/j/jeyserma/private/GEM/SingleGEM/gain/
make

# Executing rules...
echo "----------------------------------------"
echo "-------------- EXECUTE -----------------"
echo "----------------------------------------"

./gain 1_1
