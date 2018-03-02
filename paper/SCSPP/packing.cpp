/*--------------/
ALH
packing.cpp
12/10/17
/--------------*/
#include <algorithm>
#include "packing.h"
using namespace std;

int lowerBound(int maxStripWidth, double totalItemWidth){
    int lBound;
    lBound = ceil(totalItemWidth/maxStripWidth);
    return lBound;
} //End lowerBound

void optimality(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int stripSize, int LB){

    double c = static_cast<double>(LB) / stripSize;
    if(stripSize == LB){
        ++opt;
    }
    else if(c >= 0.9){
        ++opt90;
    }
    else if(c >= 0.8){
        ++opt80;
    }
    else if(c >= 0.7){
        ++opt70;
    }
    else if(c >= 0.6){
        ++opt60;
    }
    else if(c >= 0.5){
        ++opt50;
    }
    else{
        ++optLow;
    }
} //End optimality

// FFD checking vicinal sum constraint for both sides of each item.
void basicFFD(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem, int maxItemWidth,
              int maxStripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
              vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){
    int i, j, mini;
    int min = 0;
    int max = maxItemWidth;
    vector<int> itemDecrease;
    vector<int> checked(numScores, 0);

    while(itemDecrease.size() < numItem) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (itemWidths[i][partners[i]] > min && itemWidths[i][partners[i]] <= max) {
                min = itemWidths[i][partners[i]];
                mini = i;
            }
        }
        itemDecrease.push_back(mini);
        checked[mini] = 1;
        checked[partners[mini]] = 1;
        max = min;
        min = 0;
    }

    strip[0].push_back(itemDecrease[0]);
    strip[0].push_back(partners[itemDecrease[0]]);
    stripSum[0] += itemWidths[itemDecrease[0]][partners[itemDecrease[0]]];

    for(j = 1; j < itemDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][itemDecrease[j]] == 1){
                        strip[i].push_back(itemDecrease[j]);
                        strip[i].push_back(partners[itemDecrease[j]]);
                        stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][partners[itemDecrease[j]]] == 1){
                        strip[i].push_back(partners[itemDecrease[j]]);
                        strip[i].push_back(itemDecrease[j]);
                        stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(itemDecrease[j]);
                strip[i].push_back(partners[itemDecrease[j]]);
                stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                break;
            }
        }
    }

    int k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    int l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalItemWidth);

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);

} //End basicFFD


// Packing each strip in turn, choosing smallest score width that meets vicinal sum constraint.
void pairSmallest(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int maxStripWidth,
                  double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths,
                  vector<int> &stripSum, vector<vector<int> > &strip){
    int i, j;
    int count = 0;
    vector<int> scoreIncrease;
    vector<int> checked(numScores, 0);

    for(i = 0; i < numScores; ++i){
        scoreIncrease.push_back(i);
    }
    for(i = 0; i < scoreIncrease.size() - 1; ++i){
        for(j = i+1; j < scoreIncrease.size(); ++j){
            if(allScores[scoreIncrease[i]] == allScores[scoreIncrease[j]]){
                if(allScores[partners[scoreIncrease[i]]] < allScores[partners[scoreIncrease[j]]]){
                    swap(scoreIncrease[i], scoreIncrease[j]);
                }
            }
            else{
                break;
            }
        }
    }

    j = 0;
    while (count < numScores){
        for(i = 0; i < scoreIncrease.size(); ++i){
            if(checked[scoreIncrease[i]] == 1){
                continue;
            }
            if(strip[j].empty()){
                strip[j].push_back(scoreIncrease[i]);
                strip[j].push_back(partners[scoreIncrease[i]]);
                stripSum[j] += itemWidths[scoreIncrease[i]][partners[scoreIncrease[i]]];
                checked[scoreIncrease[i]] = 1;
                checked[partners[scoreIncrease[i]]] = 1;
                count += 2;
                i = -1;
                continue;
            }
            else if(!strip[j].empty()){
                if(adjMatrix[strip[j].back()][scoreIncrease[i]] == 1 && stripSum[j] + itemWidths[scoreIncrease[i]][partners[scoreIncrease[i]]] <= maxStripWidth){
                    strip[j].push_back(scoreIncrease[i]);
                    strip[j].push_back(partners[scoreIncrease[i]]);
                    stripSum[j] += itemWidths[scoreIncrease[i]][partners[scoreIncrease[i]]];
                    checked[scoreIncrease[i]] = 1;
                    checked[partners[scoreIncrease[i]]] = 1;
                    count += 2;
                    i = -1;
                    continue;
                }
            }
        }
        ++j;
    }

    int k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    int l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalItemWidth);

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);

} //End pairSmallest


// FFD including AHCA, instead of attempting to place item on end of strip, run AHCA to find feasible solution.
void FFDincAHCA(int &cp, int &na, int &type0, int &type1, int &type2, int &type3, int instance, int tau, int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem, int maxItemWidth,
                int maxStripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, mini;
    int min = 0;
    int max = maxItemWidth;
    bool feasible;
    vector<int> itemDecrease;
    vector<int> checked(numScores, 0);

    while(itemDecrease.size() < numItem) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (itemWidths[i][partners[i]] > min && itemWidths[i][partners[i]] <= max) {
                min = itemWidths[i][partners[i]];
                mini = i;
            }
        }
        itemDecrease.push_back(mini);
        checked[mini] = 1;
        checked[partners[mini]] = 1;
        max = min;
        min = 0;
    }

    strip[0].push_back(itemDecrease[0]);
    strip[0].push_back(partners[itemDecrease[0]]);
    stripSum[0] += itemWidths[itemDecrease[0]][partners[itemDecrease[0]]];

    for(j = 1; j < itemDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= maxStripWidth){
                    feasible = false;
                    AHCA(cp, na, type0, type1, type2, type3, instance, tau, i, j, feasible, allScores, partners, adjMatrix, itemWidths, itemDecrease, stripSum, strip);
                    if(feasible){
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(itemDecrease[j]);
                strip[i].push_back(partners[itemDecrease[j]]);
                stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                break;
            }
        }
    }

    int k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    int l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalItemWidth);

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);

} //End FFDincAHCA


// Initializing instance for AHCA using items on the strip from FFDincAHCA and the item to be packed.
void initializeInstance(int tau, int nScores, vector<vector<int> > &adjMat, vector<int> &scores, vector<int> &order, vector<int> &partnersX){

    int i, j;
    int vacant = 999;
    vector<int> invOrder(nScores);

    for (i = 0; i < nScores; ++i) {
        order.push_back(i);
    }

    for (i = 1; i < nScores; ++i) {
        for (j = i - 1; j >= 0; --j) {
            if (scores[i] < scores[order[j]]) {
                order[j + 1] = order[j];
                order[j] = i;
            }
        }
    }

    for (i = 0; i < nScores; ++i) {
        invOrder[order[i]] = i;
    }

    for (i = 0; i < nScores - 1; i += 2) {
        adjMat[invOrder[i]][invOrder[i + 1]] = 2;
        adjMat[invOrder[i + 1]][invOrder[i]] = 2;
    }

    sort(scores.begin(), scores.end());

    for (i = 0; i < scores.size() - 1; ++i) {
        for (j = i + 1; j < scores.size(); ++j) {
            if (scores[i] + scores[j] >= tau && adjMat[i][j] != 2) {
                adjMat[i][j] = 1;
                adjMat[j][i] = 1;
            }
        }

    }

    for (i = 0; i < nScores; ++i) {
        if(partnersX[i] == vacant) {
            for (j = 0; j < nScores; ++j) {
                if (adjMat[i][j] == 2) {
                    partnersX[i] = j;
                    partnersX[j] = i;
                    break;
                }
            }
        }
    }

} //End initializeInstance


// Modified Maximum Cardinality Matching (MMCM) Algorithm.
void MMCM(int nScores, int &matchSize, vector<vector<int> > &adjMat, vector<int> &partnersX, vector<int> &matchList, vector<int> &cycleVertex){

    int i, j;
    int vacant = 999;
    int vacantFlag = 0;
    int lastMatch = vacant;

    for (i = 0; i < nScores; ++i) {
        vacantFlag = 0;
        if (matchList[i] == vacant) {
            for (j = nScores - 1; j > i; --j) {
                if (adjMat[i][j] == 1 && matchList[j] == vacant) {
                    matchList[i] = j;
                    matchList[j] = i;
                    lastMatch = i;
                    ++matchSize;
                    if (vacantFlag == 1) {
                        cycleVertex[i] = vacant;
                        cycleVertex[j] = vacant;
                    }
                    break;
                }
                else if (adjMat[i][j] == 2 && matchList[j] == vacant) {
                    vacantFlag = 1;
                }
            }
            if (matchList[i] == vacant) {
                if ((matchList[partnersX[i]] == vacant) && (lastMatch != vacant)
                    && (partnersX[i] > i) && (adjMat[lastMatch][partnersX[i]] == 1)) {
                    matchList[i] = matchList[lastMatch];
                    matchList[lastMatch] = partnersX[i];
                    matchList[partnersX[i]] = lastMatch;
                    matchList[matchList[i]] = i;
                    cycleVertex[lastMatch] = vacant;
                    cycleVertex[partnersX[i]] = vacant;
                    lastMatch = i;
                    ++matchSize;
                }
            }
        }
    }

} //End MMCM


// Matching-Partner Structure (MPS).
void MPS(int nScores, int &nCycles, vector<int> &partnersX, vector<int> &matchList, vector<vector<int> > &mpStructure){

    int i, current;
    int smallest = nScores - 2;
    vector<int> temp;
    vector<int> checked(nScores, 0);

    do {
        current = smallest;
        do {
            temp.push_back(current);
            checked[current] = 1;
            temp.push_back(partnersX[current]);
            checked[partnersX[current]] = 1;
            current = matchList[partnersX[current]];
        } while (current != smallest);

        mpStructure.push_back(temp);
        temp.clear();

        for (i = 0; i < nScores; ++i) {
            if (checked[i] == 0) {
                smallest = i;
                break;
            }
        }
    } while (smallest != current);

    nCycles = mpStructure.size();

} //End MPS


// Bridge Recognition (BR) Algorithm.
void BR(int &qstar, int matchSize, vector<vector<int> > adjMat, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &edge, vector<vector<int> > &mpStructure,
        vector<vector<int> > &C, vector<vector<int> > &S){

    int i, j, k, nEdges;
    int vacant = 999;
    vector<int> temp;

    for (i = 0; i < mpStructure.size(); ++i) {
        for (j = 0; j < mpStructure[i].size(); ++j) {
            if (cycleVertex[mpStructure[i][j]] != vacant) {
                cycleVertex[mpStructure[i][j]] = i;
            }
        }
    }

    for (i = 0; i < matchSize; ++i) {
        while (cycleVertex[i] == vacant) {
            ++i;
        }
        edge.push_back(i);
    }
    nEdges = edge.size();

    qstar = -1;
    k = 0;
    do {
        while (k < nEdges - 2 && (adjMat[edge[k]][matchList[edge[k + 1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k + 1]])) {
            ++k;
        }
        if (adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k + 1]]) {
            ++qstar;
            temp.push_back(edge[k]);
            S[qstar][cycleVertex[edge[k]]] = 1;
            while (k < nEdges - 1 && adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && S[qstar][cycleVertex[edge[k + 1]]] == 0) {
                ++k;
                temp.push_back(edge[k]);
                S[qstar][cycleVertex[edge[k]]] = 1;
            }
            C.push_back(temp);
            temp.clear();
        }
        ++k;
    } while (k < nEdges - 1);

} //End BR


// Connecting Procedure (CP).
void CP(int &cp, int &na, int &type0, int &type1, int &type2, int &type3, int instance, int j1, int nScores, int nComp, bool &feasible, int qstar, int nCycles, vector<int> &partnersX, vector<int> &matchList,
        vector<int> &cycleVertex, vector<int> &edge, vector<vector<int> > &adjMat, vector<vector<int> > &C, vector<vector<int> > &S, vector<int> &altHam){

    int a, i, j, k, q, u, v, SSum, SqIntS;
    int v1 = 0;
    int v2 = 0;
    int vacant = 999;
    int full = vacant;
    int current = nScores - 2;
    int maxRowSize = 0;
    int maxRow;
    int type;
    //int nEdges = edge.size();
    vector<int> temp;
    vector<int> SSet2;
    vector<int> SSet3;
    vector<int> edgeCopy;
    vector<int> connectML; //was patchML;
    vector<int> QSet(nComp, 0);
    vector<int> inCycle(nScores, 0);
    vector<int> connectCycle(nComp, vacant); //was patchCycleX
    vector<vector<int> > Cconnect;

    ++cp;

    for (i = 0; i < C.size(); ++i) {
        if (C[i].size() == nCycles) {
            full = i;
            break;
        }
    }


    if (full != vacant) {

        copy(matchList.begin(), matchList.end(), back_inserter(connectML));

        for (v = 0; v < C[full].size() - 1; ++v) {
            connectML[C[full][v]] = matchList[C[full][v + 1]];
            connectML[matchList[C[full][v + 1]]] = C[full][v];
        }
        connectML[C[full][C[full].size() - 1]] = matchList[C[full][0]];
        connectML[matchList[C[full][0]]] = C[full][C[full].size() - 1];

        do {
            altHam.push_back(current);
            inCycle[current] = 1;
            altHam.push_back(partnersX[current]);
            inCycle[partnersX[current]] = 1;
            if (inCycle[connectML[partnersX[current]]] == 0) {
                current = connectML[partnersX[current]];
            }
            else {
                current = matchList[partnersX[current]];
            }
        } while (altHam.size() < nScores);

        altHam.erase(altHam.begin(), altHam.begin() + 2);
        feasible = true;
        ++type0;

    }


    else {
        type = 0;

        //region TYPE 1
        //TYPE 1: Searching for two C-cycles that connect all MPS cycles, S rows only intersect once.
        SSum = 0;
        for(a = 0; a < C.size() - 1; ++a){
            for(q = a+1; q < C.size(); ++q){
                SqIntS = 0;
                for(i = 0; i < nCycles; ++i){
                    if(S[a][i] + S[q][i] == 0){
                        SqIntS = 0;
                        break;
                    }
                    if(S[a][i] + S[q][i] == 2){
                        ++SqIntS;
                    }
                }
                if(SqIntS == 1){
                    v1 = a;
                    v2 = q;
                    SSum = nCycles;
                    break;
                }
            }
            if(SSum == nCycles){
                type = 1;
                break;
            }
        } //End Type 1
        //endregion

        //region TYPE 2
        //TYPE 2: Search for longest C-cycle, then find other C-cycles that intersect once with the longest cycle and cover cycles not yet connected.
        if(type == 0){
            temp.clear();
            SSum = 0;
            copy(edge.begin(), edge.end(), back_inserter(edgeCopy));
            for(i = 0; i < C.size(); ++i){
                if(C[i].size() > maxRowSize){
                    maxRowSize = C[i].size();
                    maxRow = i;
                }
            }

            for(j = 0; j < C[maxRow].size(); ++j){
                temp.push_back(C[maxRow][j]);
            }
            Cconnect.push_back(temp);
            temp.clear();

            for(j = 0; j < nCycles; ++j){
                SSet2.push_back(S[maxRow][j]);
            }
            for(j = 0; j < nCycles; ++j) {
                SSum = SSum + SSet2[j];
            }

            for(i = 0; i < Cconnect[0].size(); ++i){
                for(j = 0; j < edgeCopy.size(); ++j){
                    if(edgeCopy[j] == Cconnect[0][i]){
                        edgeCopy.erase(edgeCopy.begin()+j);
                        break;
                    }
                }
            }

            int nEdgesC = edgeCopy.size();
            k = 0;

            do{
                while(k < nEdgesC - 2 && (adjMat[edgeCopy[k]][matchList[edgeCopy[k+1]]] != 1 || cycleVertex[edgeCopy[k]] == cycleVertex[edgeCopy[k+1]])){
                    ++k;
                }
                if(adjMat[edgeCopy[k]][matchList[edgeCopy[k+1]]] == 1 && cycleVertex[edgeCopy[k]] != cycleVertex[edgeCopy[k+1]]
                   && ((SSet2[cycleVertex[edgeCopy[k]]] == 0 && SSet2[cycleVertex[edgeCopy[k+1]]] == 1)
                   || (SSet2[cycleVertex[edgeCopy[k]]] == 1 && SSet2[cycleVertex[edgeCopy[k+1]]] == 0))){
                    temp.push_back(edgeCopy[k]);
                    temp.push_back(edgeCopy[k+1]);
                    SSet2[cycleVertex[edgeCopy[k]]] = 1;
                    SSet2[cycleVertex[edgeCopy[k+1]]] = 1;
                    ++SSum;
                    if(SSum < nCycles){
                        ++k;
                        while(k < nEdgesC - 1 && SSet2[cycleVertex[edgeCopy[k+1]]] == 0 && adjMat[edgeCopy[k]][matchList[edgeCopy[k+1]]] == 1){
                            ++k;
                            temp.push_back(edge[k]);
                            SSet2[cycleVertex[edgeCopy[k]]] = 1;
                            ++SSum;
                        }
                    }
                    Cconnect.push_back(temp);
                    temp.clear();
                }
                ++k;
            } while (k < nEdgesC - 1 && SSum < nCycles);
            if(SSum == nCycles){
                type = 2;
            }

        }// End find overlaps
        //endregion

        //region TYPE 3
        //TYPE 3: Original
        if(type == 0){
            q = 0;
            QSet[0] = 1;
            SSum = 0;
            for (i = 0; i < nCycles; ++i) {
                SSet3.push_back(S[q][i]);
            }
            for (i = 0; i < nCycles; ++i) {
                SSum = SSum + SSet3[i];
            }
            if (SSum >= 1) {
                connectCycle[q] = 1;
            }
            while (q <= qstar && SSum < nCycles) {
                do {
                    ++q;
                    SqIntS = vacant;
                    if (q <= qstar) {
                        for (j = 0; j < nCycles; ++j) {
                            if (S[q][j] == 1 && SSet3[j] == 1) {
                                SqIntS = 1;
                                break;
                            }
                        }
                    }
                } while (q < qstar + 1 && (QSet[q] == 1 || SqIntS == vacant));

                if (q <= qstar) {
                    for (i = 0; i < nCycles; ++i) {
                        if (SSet3[i] == 0 && S[q][i] == 1) {
                            SSet3[i] = 1;
                            ++SSum;
                            connectCycle[q] = 1;
                        }
                    }
                    QSet[q] = 1;
                    q = 0;
                }
            }
            if(SSum == nCycles){
                type = 3;
            }
        } //End Type 3
        //endregion

        if (SSum == nCycles) {
            temp.clear();
            if(type == 1){
                ++type1;
                for (j = 0; j < C[v1].size(); ++j) {
                    temp.push_back(C[v1][j]);
                }
                Cconnect.push_back(temp);
                temp.clear();

                for (j = 0; j < C[v2].size(); ++j) {
                    temp.push_back(C[v2][j]);
                }
                Cconnect.push_back(temp);
                temp.clear();
            }
            else if(type == 2){
                ++type2;
            }
            else if(type == 3){
                ++type3;
                for (i = 0; i < connectCycle.size(); ++i) {
                    if (connectCycle[i] == 1) {
                        for (j = 0; j < C[i].size(); ++j) {
                            temp.push_back(C[i][j]);
                        }
                        Cconnect.push_back(temp);
                        temp.clear();
                    }
                }

            }
            else{
                cout << "[ERROR]: NO TYPE.\n";
                exit(1);
            }


            copy(matchList.begin(), matchList.end(), back_inserter(connectML));

            for (u = 0; u < Cconnect.size(); ++u) {
                for (v = 0; v < Cconnect[u].size() - 1; ++v) {
                    connectML[Cconnect[u][v]] = matchList[Cconnect[u][v + 1]];
                    connectML[matchList[Cconnect[u][v + 1]]] = Cconnect[u][v];
                }
                connectML[Cconnect[u][Cconnect[u].size() - 1]] = matchList[Cconnect[u][0]];
                connectML[matchList[Cconnect[u][0]]] = Cconnect[u][Cconnect[u].size() - 1];
            }

            current = nScores - 2;
            do {
                altHam.push_back(current);
                inCycle[current] = 1;
                altHam.push_back(partnersX[current]);
                inCycle[partnersX[current]] = 1;
                if (inCycle[connectML[partnersX[current]]] == 0) {
                    current = connectML[partnersX[current]];
                }
                else {
                    current = matchList[partnersX[current]];
                }
            } while (altHam.size() < nScores);

            altHam.erase(altHam.begin(), altHam.begin() + 2);
            feasible = true;
        }

        else {
            ++na;
            feasible = false;
        }
    }

} //End CP


// Alternating Hamiltonian Construction Algorithm (AHCA).
void AHCA(int &cp, int &na, int &type0, int &type1, int &type2, int &type3, int instance, int tau, int i1, int j1, bool &feasible, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
          vector<vector<int> > &itemWidths, vector<int> &itemDecrease, vector<int> &stripSum, vector<vector<int> > &strip){

    int k;
    vector<int> scores;
    vector<int> original;
    feasible = false;

    //Creating scores vector
    for(k = 0; k < strip[i1].size(); ++k){
        scores.push_back(allScores[strip[i1][k]]);
        original.push_back(strip[i1][k]);
    }
    scores.push_back(allScores[itemDecrease[j1]]);
    original.push_back(itemDecrease[j1]);
    scores.push_back(allScores[partners[itemDecrease[j1]]]);
    original.push_back(partners[itemDecrease[j1]]);
    scores.push_back(tau);
    scores.push_back(tau);

    //Variables
    int i, j, qstar;
    int once = 0;
    int vacant = 999;
    int matchSize = 0;
    int nCycles = 0;
    int nScores = scores.size();
    int nItem = scores.size() / 2;
    int nComp = (nItem + (nItem % 2)) / 2;
    vector<int> order;
    vector<int> altHam;
    vector<int> final;
    vector<int> cycleVertex(nScores, 1);
    vector<int> partnersX(nScores, vacant);
    vector<int> matchList(nScores, vacant);
    vector<int> edge;
    vector<vector<int> > C;
    vector<vector<int> > mpStructure;
    vector<vector<int> > S(nComp, vector<int>(nComp, 0));
    vector<vector<int> > adjMat(nScores, vector<int>(nScores, 0));

    do {
        once = 1;
        initializeInstance(tau, nScores, adjMat, scores, order, partnersX);

        MMCM(nScores, matchSize, adjMat, partnersX, matchList, cycleVertex);
        if (matchSize < nItem) {
            feasible = false;
            break;
        }

        MPS(nScores, nCycles, partnersX, matchList, mpStructure);
        if (mpStructure[0].size() == nScores) {
            for (j = 0; j < mpStructure[0].size(); ++j) {
                altHam.push_back(mpStructure[0][j]);
            }

            altHam.erase(altHam.begin(), altHam.begin() + 2);

            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            feasible = true;
            break;
        }

        BR(qstar, matchSize, adjMat, matchList, cycleVertex, edge, mpStructure, C, S);
        if (qstar == -1) {
            feasible = false;
            break;
        }

        CP(cp, na, type0, type1, type2, type3, instance, j1, nScores, nComp, feasible, qstar, nCycles, partnersX, matchList, cycleVertex, edge, adjMat, C, S, altHam);
        if (feasible) {
            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            break;
        }

    } while (once == 0);

    if(feasible){
        stripSum[i1] += itemWidths[itemDecrease[j1]][partners[itemDecrease[j1]]];
        strip[i1].swap(final);
    }

} //End AHCA



























































































