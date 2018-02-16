/*--------------
ALH
Fixing bug in mbahra from SCSPP
16/02/2018
---------------*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

int main() {

    int i, j, k;
    int threshold = 70;
    int vacant = 999;
    vector<int> scores = { 1, 69, 8, 3, 70, 13, 68, 10, 68, 6, 70, 30, 62, 45, 57, 5, 19, 43, 70, 70 };
    vector<int> order;
    vector<int> final;

    int nScores = scores.size();
    int nBox = nScores /2;
    int nComp = (nBox + (nBox % 2)) / 2;
    vector<int> invOrder(nScores);
    vector<vector<int> > adjMat(nScores, vector<int>(nScores, 0));
    vector<int> mates(nScores, vacant);

    for(k = 0; k < nScores; ++k){
        order.push_back(k);
    }

    for (i = 1; i < nScores; ++i) {
        for (j = i - 1; j >= 0; --j) {
            if (scores[i] < scores[order[j]]) {
                order[j + 1] = order[j];
                order[j] = i;
            }
        }
    }

    cout << "Order:\n";
    for(int v : order){
        cout << v << " ";
    }
    cout << endl;

    for (k = 0; k < nScores; ++k) {
        invOrder[order[k]] = k;
    }

    cout << "invOrder:\n";
    for(int v : invOrder){
        cout << v << " ";
    }
    cout << endl;

    for(i = 0; i < nScores; i+=2){
        adjMat[invOrder[i]][invOrder[i+1]] = 2;
        adjMat[invOrder[i+1]][invOrder[i]] = 2;
    }

    sort(scores.begin(), scores.end());

    for(i = 0; i < scores.size() - 1; ++i){
        for(j = i + 1; j < scores.size(); ++j){
            if(scores[i] + scores[j] >= threshold && adjMat[i][j] != 2){
                adjMat[i][j] = 1;
                adjMat[j][i] = 1;
            }
        }
    }


    for(i = 0; i < nScores; ++i){
        if(mates[i] == vacant) {
            for (j = 0; j < nScores; ++j) {
                if (adjMat[i][j] == 2) {
                    mates[i] = j;
                    mates[j] = i;
                    break;
                }
            }
        }
    }

    cout << " Mates:\n";
    for(int v : mates){
        cout << v << " ";
    }
    cout << endl;


    //MTGMA
    int vacantFlag = 0;
    int matchSize = 0;
    int lastMatch = vacant;
    int mateMatch = vacant;
    vector<int> cycleVertex(nScores, 1);
    vector<int> matchList(nScores, vacant);

    for(i = 0; i < nScores; ++i){
        vacantFlag = 0;
        if(matchList[i] == vacant){
            for(j = nScores -1; j > i; --j){
                if(adjMat[i][j] == 1 && matchList[j] == vacant){
                    matchList[i] = j;
                    matchList[j] = i;
                    lastMatch = i;
                    ++matchSize;
                    if(vacantFlag == 1){
                        cycleVertex[i] = vacant;
                        cycleVertex[j] = vacant;
                    }
                    break;
                }
                else if(adjMat[i][j] == 2 && matchList[j] == vacant){
                    vacantFlag = 1;
                }
            }
            if(matchList[i] == vacant){
                mateMatch = mates[i]; //check if matches line 588 in SCSPP.
            }
        }
    }






















}