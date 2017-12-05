/*--------------/
ALH
base.h
Evolutionary Algorithm with Local Search
05/12/2017
/--------------*/

#ifndef BASE_H
#define BASE_H

#include <iostream>
#include <vector>
using namespace std;

void createInstance(int threshold, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, int numScores, int numBox, double &totalBoxWidth, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);




#endif
