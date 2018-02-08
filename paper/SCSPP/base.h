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

void resetVectors(int numScores, int numBox, vector<int> &allScores, vector<int> &mates,
                  vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip);

void output(int opt, int opt90, int opt80, int opt70, int opt60, int opt50, int optLow, int numInstances);

void createInstance(int numScores, int numBox, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, double &totalBoxWidth,
                    vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths);

void createInstanceUser(int threshold, int numScores, double &totalBoxWidth, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes);



#endif
