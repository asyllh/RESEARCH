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

void resetVectors(int numScores, int numItem, vector<int> &allScores, vector<int> &partners,
                  vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);

void output(int opt, int opt90, int opt80, int opt70, int opt60, int opt50, int optLow, int numInstances);

void createInstance(int instance, int tau, int numScores, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, double &totalItemWidth,
                    vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths);



#endif
