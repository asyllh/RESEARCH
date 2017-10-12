/*--------------/
ALH
base.h
12/10/17
/--------------*/

#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <vector>
using namespace std;

void resetVectors(int vacant, int numScores, int numComp, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList, vector<int> &mates, vector<vector<int> > &S, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void clearVectors(vector<int> &allScores, vector<vector<int> > &mateInduced, vector<int> &lengthMateInduced, vector<vector<int> > &T, vector<int> &fullCycle, vector<int> &completePath);

void createInstance(int threshold, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, int numScores, int numBox, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void createInstanceUser(int threshold, int numScores, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);



#endif
