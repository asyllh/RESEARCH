/*--------------/
ALH
base.h
Combined Program with Heuristics and EA
17/03/18
/--------------*/

#ifndef COMBINED_BASE_H
#define COMBINED_BASE_H

#include <iostream>
#include <vector>
using namespace std;

void ResetVar(int numScores, int numItem, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
              vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip);

void Output(int opt, int opt90, int opt80, int opt70, int opt60, int opt50, int optLow, int numInstances);

void CreateInstance(int tau, int numScores, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, double &totalItemWidth,
                    vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths);

#endif
