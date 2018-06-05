/*--------------/
ALH
base.cpp
12/10/17
/--------------*/
#include <algorithm>
#include <iomanip>
#include "base.h"
using namespace std;

void ResetVar(int numScores, int numItem, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
              vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){
    int i, j;

    for(i = 0; i < numScores; ++i){
        partners[i] = 0;
        for(j = 0; j < numScores; ++j){
            adjMatrix[i][j] = 0;
            itemWidths[i][j] = 0;
        }
    }

    strip.clear();
    strip.resize(numItem);
    stripSum.clear();
    stripSum.resize(numItem, 0);
    allScores.clear();

} //End ResetVar



void CreateInstance(int tau, int numScores, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, double &totalItemWidth,
                    vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths){

    int i, j;
    int count = 1;
    vector<int> randOrder;
    vector<int> checkItem(numScores, 0);
    totalItemWidth = 0.0;

    for (i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }

    sort(allScores.begin(), allScores.end());

    //int ones = 0;
    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= tau) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
                //ones += 2;
            }
        }

    }

    //ones -= numScores;
    //cout << tau << ": " << ones << endl;

    for (i = 0; i < numScores; ++i) {
        randOrder.push_back(i);
    }

    random_shuffle(randOrder.begin(), randOrder.end());

    for (i = 0; i < numItem; ++i) {
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                partners[i] = j;
                break;
            }
        }
    }

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && itemWidths[i][j] == 0){
                itemWidths[i][j] = rand() % (maxItemWidth - minItemWidth + 1) + minItemWidth;
                itemWidths[j][i] = itemWidths[i][j];
                break;
            }

        }
    }


    //cout << "Box#" << setw(10) << "Mates" << setw(10) << "Scores" << setw(10) << "Width\n";
    for(i = 0; i < numScores; ++i){
        if(checkItem[i] == 1){
            continue;
        }
        //cout << count << setw(10)  << i << "-" << mates[i] << setw(10) << allScores[i] << "-" << allScores[mates[i]] << setw(8) << itemWidths[i][mates[i]] << endl;
        totalItemWidth += itemWidths[i][partners[i]];
        checkItem[i] = 1;
        checkItem[partners[i]] = 1;
        ++count;

    }

    //cout << "TotalItemWidth: " << totalItemWidth << endl;

} //End CreateInstance