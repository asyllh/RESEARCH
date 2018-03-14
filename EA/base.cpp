/*--------------/
ALH
base.cpp
Evolutionary Algorithm with Local Search
05/12/2017
06/03/2018
/--------------*/
#include <algorithm>
#include <iomanip>
#include "base.h"
using namespace std;

void createInstance(int tau, int numScores, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, double &totalItemWidth,
                    vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths,
                    vector<vector<int> > &allItems){

    int i, j, k;
    int count = 1;
    vector<int> randOrder;
    vector<int> checkItem(numScores, 0);
    totalItemWidth = 0.0;

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for (i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }


    //Sort all of the scores in the allScores vector in ascending order
    sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    //cout << "All scores:\n";
    /*for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl;*/

    //Filling in adjacency matrix - if sum of two scores >= tau (70), then insert 1 into the matrix, else leave as 0
    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= tau) {
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

    //Assign partners to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for (i = 0; i < numItem; ++i) {
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
                partners[i] = j;
                break;
            }
        }
    }
    /*cout << "Mates Vector:\n";
    for(i = 0; i < partners.size(); ++i){
        cout << partners[i] << " ";
    }
    cout << endl << endl;*/

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && itemWidths[i][j] == 0){
                itemWidths[i][j] = rand() % (maxItemWidth - minItemWidth + 1) + minItemWidth;
                itemWidths[j][i] = itemWidths[i][j];
                break;
            }

        }
    }

    k = 1;
    for(i = 0; i < numScores; ++i){
        for(j = i+1; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                allItems[i][j] = k;
                allItems[j][i] = k * -1;
                ++k;
                break;
            }
        }
    }
    //cout << endl;

    cout << right << setw(5) << "Box#" << setw(12) << "Scores" << setw(12) << "Mates" << setw(12) << "Width\n";
    for(i = 0; i < numScores; ++i){
        if(checkItem[i] == 1){
            continue;
        }
        cout << setw(5) << count << setw(10) << allScores[i] << "-" << allScores[partners[i]] << setw(10) << i  << "-" << partners[i] << setw(10) << itemWidths[i][partners[i]] << endl;
        totalItemWidth += itemWidths[i][partners[i]];
        checkItem[i] = 1;
        checkItem[partners[i]] = 1;
        ++count;

    }

    cout << "Total Item Widths: " << totalItemWidth << endl << endl;

    /*cout << "allItems:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allItems[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}
