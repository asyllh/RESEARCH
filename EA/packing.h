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

int lowerBound(double totalItemWidth, int maxStripWidth);

double fitness (int maxStripWidth, vector<int> &stripSum, vector<vector<int> > &strip);

void FFDecreasing(int numScores, int numItem, int maxItemWidth, vector<int> &partners, vector<vector<int> > &itemWidths, vector<int> &itemOrder);

void FFRandom(int numScores, int numItem, vector<int> &partners, vector<vector<int> > &itemWidths, vector<int> &itemOrder);

void FFShell(int numScores, int numItem, int maxItemWidth, int maxStripWidth, vector<int> &partners,
         vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip, bool decrease);


void partialFFD(int numScores, int maxItemWidth, int maxStripWidth, vector<int> &partners, vector<vector<int> > &adjMatrix,
                vector<vector<int> > &itemWidths, vector<int> &partialItem, vector<int> &partialSum, vector<vector<int> > &partialSol);


void createInitialPopulation(int numScores, int numItem, int maxItemWidth, int maxStripWidth, vector<int> &allScores,
                             vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths,
                             vector<vector<int> > &populationSum, vector<vector<vector<int> > > &population);


void mutation(int numScores, int maxItemWidth, int maxStripWidth, vector<int> &allScores, vector<int> &partners,
              vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);


void localSearch(int numScores, int maxItemWidth, int maxStripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                 vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY);

void initAHCA(int swapType, int moveType, int &feasible, int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores,
              vector<vector<int> > &itemWidths, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY);

void AHCA(int swapType, int moveType, int &feasible, int i1, int a1, int b1, int j1, int c1, int d1,
          vector<int> &allScores,
          vector<vector<int> > &itemWidths, vector<int> &stripSumX, vector<vector<int> > &stripX,
          vector<int> &stripSumY, vector<vector<int> > &stripY);

void EA(int numScores, int maxItemWidth, int maxStripWidth, double &parent1cost, double &parent2cost, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
        vector<vector<int> > &itemWidths, vector<vector<int> > &populationSum, vector<vector<vector<int> > > &population);

void GGA(int numScores, int maxItemWidth, int maxStripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
         vector<vector<int> > &itemWidths, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY);
#endif

