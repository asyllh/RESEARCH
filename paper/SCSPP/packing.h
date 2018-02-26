/*--------------/
ALH
packing.h
12/10/17
/--------------*/
#ifndef PACKING_H
#define PACKING_H

#include <iostream>
#include <vector>
using namespace std;


int lowerBound(int maxStripWidth, double totalItemWidth);

void optimality(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int stripSize, int LB);

// FFD checking vicinal sum constraint for both sides of each item.
void basicFFD(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem, int maxItemWidth,
                         int maxStripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                         vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


// Packing each strip in turn, choosing smallest score width that meets vicinal sum constraint.
void pairSmallest(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores,
                           int maxStripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                           vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


// FFD including AHCA, instead of attempting to place item on end of strip, run AHCA to find feasible solution.
void FFDincAHCA(int instance, int tau, int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem, int maxItemWidth,
                        int maxStripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                        vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


// Initializing instance for AHCA using items on the strip from FFDincAHCA and the item to be packed.
void initializeInstance(int tau, int nScores, vector<vector<int> > &adjMat, vector<int> &scores, vector<int> &order, vector<int> &partnersX);


// Modified Maximum Cardinality Matching (MMCM) Algorithm.
void MMCM(int nScores, int &matchSize, vector<vector<int> > &adjMat, vector<int> &partnersX, vector<int> &matchList, vector<int> &cycleVertex);


// Matching-Partner Structure (MPS).
void MPS(int nScores, int &nCycles, vector<int> &partnersX, vector<int> &matchList, vector<vector<int> > &mpStructure);


// Bridge Recognition (BR) Algorithm.
void BR(int &qstar, int matchSize, vector<vector<int> > adjMat, vector<int> &matchList, vector<int> &cycleVertex, vector<vector<int> > &mpStructure,
        vector<vector<int> > &C, vector<vector<int> > &S);


// Connecting Procedure (CP).
void CP(int instance, int j1, int nScores, int nComp, bool &feasible, int qstar, int nCycles, vector<int> &partnersX, vector<int> &matchList, vector<int> &cycleVertex,
        vector<vector<int> > &C, vector<vector<int> > &S, vector<int> &altHam);


// Alternating Hamiltonian Construction Algorithm (AHCA).
void AHCA(int instance, int tau, int i1, int j1, bool &feasible, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
            vector<vector<int> > &itemWidths, vector<int> &itemDecrease, vector<int> &stripSum, vector<vector<int> > &strip);

#endif

