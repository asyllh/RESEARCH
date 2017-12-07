/*--------------/
ALH
packing.h
Evolutionary Algorithm with Local Search
05/12/2017
/--------------*/
#ifndef PACKING_H
#define PACKING_H

#include <iostream>
#include <vector>
using namespace std;

void swap(int &a, int &b);

int lowerBound(double totalBoxWidth, int maxStripWidth);

void FFDecreasing(int numScores, int numBox, int maxBoxWidth, vector<int> &mates, vector<vector<int> > &boxWidths, vector<int> &boxOrder);

void FFRandom(int numScores, int numBox, vector<int> &mates, vector<vector<int> > &boxWidths, vector<int> &boxOrder);

void FFShell(int numScores, int numBox, int maxBoxWidth, int maxStripWidth, vector<int> &mates,
         vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip, bool decrease);


void partialFFD(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &mates, vector<vector<int> > &adjMatrix,
                vector<vector<int> > &boxWidths, vector<int> &partialBoxes, vector<int> &partialSum, vector<vector<int> > &partialSol);


void createInitialPopulation(int numScores, int numBox, int maxBoxWidth, int maxStripWidth, vector<int> &allScores,
                             vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths,
                             vector<vector<int> > &populationSum, vector<vector<vector<int> > > &population);


void mutation(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &allScores, vector<int> &mates,
              vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip);


void localSearch(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
                 vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY);

void MBAHRA(int swapType, int moveType, int &feasible, int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores,
            vector<vector<int> > &boxWidths, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY);
#endif

