/*--------------------/
ALH
Rewritten Minimum Score Separation Problem
17/07/2017
/--------------------*/

#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <math.h>
#include <float.h>
#include <queue>
#include <limits.h>
#include <algorithm>
#include <iomanip> //header providing parametric manipulators
using namespace std;


int main(){

    //Variables
    int i, j, r;
    int randomSeed = 2;
    int numInstances = 100;
    int numBox = 5;
    int numScores = numBox * 2;
    int minWidth = 1;
    int maxWidth = 70;
    int threshold = 70;
    vector<int> allScores(numScores, 0);
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores ,0));
    srand(randomSeed);
    time_t startTime, endTime;
    startTime = clock();

    for(i = 0; i < numScores -2; ++i){
        allScores[i] = rand() % (maxWidth - minWidth + 1) + minWidth;
    }

    allScores[numScores - 2] = 71;
    allScores[numScores - 1] = 71;

    cout << "all scores:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << endl;
    }

    sort (allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    cout << "all scores:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << endl;
    }

    for(i = 0; i < allScores.size()-1; ++i){
        for(j = i+1; j< allScores.size(); ++j){
            if(allScores[i] + allScores[j] >= threshold){
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }
    cout << "adjacency matrix\n\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << "\t";
        }
        cout << endl;
    }








	endTime = clock();
	int totalTime = (int)(((endTime - startTime) / double(CLOCKS_PER_SEC)) * 100);
	cout << "CPU Time = " << totalTime << " milliseconds.\n";

}