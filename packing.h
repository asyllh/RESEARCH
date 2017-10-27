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

void costEvalSwap(int i, int j, int k, int l, int &totalCost, int maxStripWidth, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum);

void swapBox(int i, int j, int k, int l, int &totalCost, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum);

void checkSwap(int &totalCost, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum);

void moveBox(int i, int j, int k, int l, int &moved, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum);

void checkMove(int moved, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum);


#endif

