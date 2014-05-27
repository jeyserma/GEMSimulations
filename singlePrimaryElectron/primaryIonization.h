/*
 * primaryIonization.h
 * Author: Jan Eysermans 
 * 2014
 */

// General libs
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdio.h>     
#include <stdlib.h> 
#include <math.h> 

// ROOT libs
#include <TROOT.h>

struct cluster {
    
    Double_t x0, y0, z0, t0, ec;
    Int_t nc;
    std::vector<Double_t> x1, y1, z1, t1, e1, vx1, vy1, vz1;
};


struct particle {
    
    TString type;
    Double_t x0, y0, z0, vx0, vy0, vz0, e0, t0;
    Int_t noClusters;
    Double_t lambda, stoppingPower;
    std::vector<cluster> clusters;
};
