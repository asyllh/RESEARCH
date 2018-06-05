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


int LowerBound(double totalItemWidth, int stripWidth);


void MFFD(int numScores, int numItem, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners,
          vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


// Packing each strip in turn, choosing smallest score width that meets vicinal sum constraint.
void PairSmallest(int numScores, int stripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                  vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


// FFD including AHCA, instead of attempting to place item on end of strip, run AHCA to find feasible solution.
void MFFDPlus(int tau, int numScores, int numItem, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners,
              vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


void AHCA(int tau, int i1, int j1, int &feasible, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
           vector<vector<int> > &itemWidths, vector<int> &itemDecrease, vector<int> &stripSum, vector<vector<int> > &strip);

void InitInstance(int tau, int nScores, vector<vector<int> > &adjMat, vector<int> &scores, vector<int> &order, vector<int> &partnersX);


void MCM(int nScores, int &matchSize, vector<vector<int> > &adjMat, vector<int> &partnersX, vector<int> &matchList,
         vector<int> &cycleVertex);


void MPS(int nScores, int &nCycles, vector<int> &partnersX, vector<int> &matchList, vector<vector<int> > &mpStructure);


void BR(int &qstar, int matchSize, vector<vector<int> > &adjMat, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &edge,
        vector<vector<int> > &mpStructure, vector<vector<int> > &C, vector<vector<int> > &S);

void CP(int nScores, int nComp, int &feasible, int qstar, int nCycles, vector<int> &partnersX, vector<int> &matchList,
        vector<int> &cycleVertex, vector<int> &edge, vector<vector<int> > &adjMat, vector<vector<int> > &C, vector<vector<int> > &S, vector<int> &altHam);


#endif

