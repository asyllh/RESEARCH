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

void swap(int &a, int &b);

int lowerBound(double totalBoxWidth, int maxStripWidth);

void packStripsFFD(int &totalCost, int numBox, int maxBoxWidth, int maxStripWidth, double totalBoxWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<int> &stripNumBoxes, vector<vector<int> > &strip, vector<vector<int> > &stripWidth);

int initCost(int &totalCost, int maxStripWidth, vector<int> &stripSum);

void localSearch(int maxStripWidth, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<int> &stripSumX, vector<int> &stripSumY, vector<vector<int> > &strip, vector<vector<int> > &stripX, vector<vector<int> > &stripY);

void MBAHRA(int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores, vector<vector<int> > &stripX, vector<vector<int> > &stripY);

#endif

