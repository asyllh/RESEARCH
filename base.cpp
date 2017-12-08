/*--------------/
ALH
base.cpp
Evolutionary Algorithm with Local Search
05/12/2017
/--------------*/
#include <algorithm>
#include <iomanip>
#include "base.h"
using namespace std;

void createInstance(int numScores, int numBox, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, double &totalBoxWidth,
                    vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes){

    int i, j, k;
    int threshold = 70;
    int count = 1;
    vector<int> randOrder;
    vector<int> checkBox(numScores, 0);
    totalBoxWidth = 0.0;

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for (i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }


    //Sort all of the scores in the allScores vector in ascending order
    sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    //cout << "All scores:\n";
    for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl;

    //Filling in adjacency matrix - if sum of two scores >= threshold (70), then insert 1 into the matrix, else leave as 0
    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= threshold) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    //Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for (i = 0; i < numScores; ++i) {
        randOrder.push_back(i);
    }

    //Randomly shuffle all values in randOrder vector
    random_shuffle(randOrder.begin(), randOrder.end());

    //Assign mates to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for (i = 0; i < numBox; ++i) {
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    /*cout << "AdjMatrix:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;*/



    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                mates[i] = j;
                break;
            }
        }
    }
    /*cout << "Mates Vector:\n";
    for(i = 0; i < mates.size(); ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;*/

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && boxWidths[i][j] == 0){
                boxWidths[i][j] = rand() % (maxBoxWidth - minBoxWidth + 1) + minBoxWidth;
                boxWidths[j][i] = boxWidths[i][j];
                break;
            }

        }
    }

    k = 1;
    for(i = 0; i < numScores; ++i){
        for(j = i+1; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                allBoxes[i][j] = k;
                allBoxes[j][i] = k * -1;
                ++k;
                break;
            }
        }
    }
    //cout << endl;

    cout << right << setw(5) << "Box#" << setw(12) << "Scores" << setw(12) << "Mates" << setw(12) << "Width\n";
    for(i = 0; i < numScores; ++i){
        if(checkBox[i] == 1){
            continue;
        }
        cout << setw(5) << count << setw(10) << allScores[i] << "-" << allScores[mates[i]] << setw(10) << i  << "-" << mates[i] << setw(10) << boxWidths[i][mates[i]] << endl;
        totalBoxWidth += boxWidths[i][mates[i]];
        checkBox[i] = 1;
        checkBox[mates[i]] = 1;
        ++count;

    }

    cout << "Total Box Widths: " << totalBoxWidth << endl << endl;

    /*cout << "allBoxes:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allBoxes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}
