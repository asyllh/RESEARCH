/*--------------------/
ALH
packing.h
18/08/2017
/--------------------*/
#ifndef PACKING_H
#define PACKING_H

#include <iostream>
#include <vector>
using namespace std;

void packStripsSmallest(int numScores, int numBox, int maxStripWidth, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths);

void packStripsMIS(int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &mateInduced, vector<vector<int> > &boxWidths);

void weakMatchPath(int vacant, int numScores, vector<int> &matchList, vector<int> &mates);

void packStripsBFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripsFFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripsNFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripsBFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripsFFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripsNFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripFFDScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

void packStripFFIScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths);

#endif