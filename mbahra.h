/*--------------------/
ALH
mbahra.h
18/08/2017
/--------------------*/
#ifndef MBAHRA_H
#define MBAHRA_H

#include <iostream>
#include <vector>
using namespace std;

void resetVectors(int vacant, int numScores, int numComp, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList, vector<int> &mates, vector<vector<int> > &S, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void clearVectors(vector<int> &allScores, vector<vector<int> > &mateInduced, vector<int> &lengthMateInduced, vector<vector<int> > &T, vector<int> &fullCycle, vector<int> &completePath);

void createInstance(int threshold, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, int numScores, int numBox, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void createInstanceUser(int threshold, int numScores, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void MTGMA(int vacant, int threshold, int numScores, int &matchSize, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList);

void MIS(int numScores, int &numCycles, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<int> &matchList, vector<vector<int> > &mateInduced, vector<int> &lengthMateInduced);

void FCA(int &qstar, int vacant, int matchSize, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList, vector<vector<int> > &mateInduced, vector<vector<int> > &S, vector<vector<int> > &T);

void oneTCyclePatch(int v, int full, int save, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &fullCycle, vector<vector<int> > &mateInduced, vector<vector<int> > &T);

void multipleTCyclePatch(int &u, int v, int save, int vacant, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &fullCycle, vector<vector<int> > &mateInduced, vector<vector<int> > &Tpatch, vector<int> &patchVertex);

void makePath(int numScores, vector<int> &fullCycle, vector<int> &completePath, vector<vector<int> > &boxWidths, vector<int> &allScores, vector<vector<int> > &allBoxes);

void patchGraph(int qstar, int vacant, int instance, int numScores, int numCycles, int &feasible, int &infeasible, int &fullT, int &splitT, int &noPatch, int &problem, vector<int> &matchList, vector<int> &cycleVertex, vector<vector<int> > &mateInduced, vector<vector<int> > &S, vector<vector<int> > &T, vector<int> &fullCycle, vector<int> &completePath, vector<vector<int> > &boxWidths, vector<int> &allScores, vector<vector<int> > &allBoxes);


#endif
