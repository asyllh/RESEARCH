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
    unsigned int randomSeed = 2;
    int numInstances = 100; //number of instances of mssp, use in main for loop
    unsigned int numBox = 5; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    unsigned int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    int minWidth = 1; //minimum width of score (millimeters)
    int maxWidth = 70; //maximum width of score (millimeters)
    int threshold = 70; //adjacency threshold of scores, minimum knife distance
    vector<int> allScores(numScores, 0); //vector containing all score widths
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores ,0)); //adjacency matrix, 0 if width sum < threshold, 1 if width sum >= threshold, 2 if scores are mates (either side of same box)
    srand(randomSeed); //seed

    time_t startTime, endTime; //start clock
    startTime = clock();

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for(i = 0; i < numScores -2; ++i){
        allScores[i] = rand() % (maxWidth - minWidth + 1) + minWidth;
    }
    //add two dominating vertices with score widths = 71 (these scores will be either side of same box, mates)
    allScores[numScores - 2] = 71;
    allScores[numScores - 1] = 71;

    //Print out allScores vector
    cout << "all scores:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << endl;
    }

    //Sort all of the scores in the allScores vector in ascending order
    sort (allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    //Print out allScores vector (scores now in ascending order)
    cout << "all scores in increasing order:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << endl;
    }

    //Filling in adjacency matrix - if sum of two scores >= threshold (70), then insert 1 into the matrix, else leave as 0
    for(i = 0; i < allScores.size()-1; ++i){
        for(j = i+1; j< allScores.size(); ++j){
            if(allScores[i] + allScores[j] >= threshold){
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    //Print out adjacency matrix inc threshold
    cout << "adjacency matrix\n\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << "\t";
        }
        cout << endl;
    }

    //Create vector to be used to assign mates
    vector<int> randOrder(numScores, 0);

    //initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for(i = 0; i < numScores-2; ++i){
        randOrder[i] = i;
    }
    randOrder[numScores-2] = numScores -2;
    randOrder[numScores-1] = numScores -1;

    //Randomly shuffle all values in randOrder vector EXCEPT the last two values (dominating vertices, must stay as mates)
    random_shuffle(randOrder.begin(), randOrder.begin()+8);

    //Print out randOrder vector
    cout << "random order:\n";
    for(i = 0; i < randOrder.size(); ++i){
        cout << randOrder[i] << endl;
    }
    cout << endl;

    //Assign mates to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for(i = 0; i < numBox; ++i){
        adjMatrix[randOrder[2*i]][randOrder[2*i+1]] = 2;
        adjMatrix[randOrder[2*i+1]][randOrder[2*i]] = 2;
    }

    //Print out adjacency matrix inc threshold and mates
    cout << "adjacency matrix with twin:\n\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << "\t";
        }
        cout << endl;
    }

    //MATCHING ALGORITHM











	endTime = clock();
	int totalTime = (int)(((endTime - startTime) / double(CLOCKS_PER_SEC)) * 100);
	cout << "CPU Time = " << totalTime << " milliseconds.\n";

}