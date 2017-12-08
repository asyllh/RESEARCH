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

void createInstance(int numScores, int numBox, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, double &totalBoxWidth,
                    vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);

#endif
