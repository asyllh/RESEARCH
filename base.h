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

void resetVectors(int numScores, int numBox, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes, vector<int> &stripSum, vector<vector<int> > &strip);

void clearVectors(vector<int> &allScores, vector<int> &stripNumBoxes, vector<int> &stripSumX, vector<int> &stripSumY, vector<vector<int> > &stripX, vector<vector<int> > &stripY);

void createInstance(int threshold, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, int numScores, int numBox, double &totalBoxWidth, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

void createInstanceUser(int threshold, int numScores, double &totalBoxWidth, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);



#endif
