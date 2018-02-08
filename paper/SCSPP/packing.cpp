/*--------------/
ALH
packing.cpp
12/10/17
/--------------*/
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "packing.h"
using namespace std;

void swap(int &a, int &b){
    int temp = a;
    a = b;
    b = temp;
}

int lowerBound(int maxStripWidth, double totalBoxWidth){
    int lBound;

    lBound = ceil(totalBoxWidth/maxStripWidth);
    return lBound;
}

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

}


void packStripsFFDApprox(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numBox, int maxBoxWidth,
                         int maxStripWidth, double totalBoxWidth, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
                         vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, mini, k, l;
    int min = 0;
    int max = maxBoxWidth;
    vector<int> boxDecrease;
    vector<int> checked(numScores, 0);



    while(boxDecrease.size() < numBox) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] <= max) {
                min = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxDecrease.push_back(mini);
        checked[mini] = 1;
        checked[mates[mini]] = 1;
        max = min;
        min = 0;
    }

    /*cout << "Box decrease:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;*/

    strip[0].push_back(boxDecrease[0]);
    strip[0].push_back(mates[boxDecrease[0]]);
    stripSum[0] += boxWidths[boxDecrease[0]][mates[boxDecrease[0]]];


    for(j = 1; j < boxDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[boxDecrease[j]][mates[boxDecrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][boxDecrease[j]] == 1){
                        strip[i].push_back(boxDecrease[j]);
                        strip[i].push_back(mates[boxDecrease[j]]);
                        stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][mates[boxDecrease[j]]] == 1){
                        strip[i].push_back(mates[boxDecrease[j]]);
                        strip[i].push_back(boxDecrease[j]);
                        stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(boxDecrease[j]);
                strip[i].push_back(mates[boxDecrease[j]]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                break;
            }
        }
    }

    k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalBoxWidth);

    /*cout << "After FFD: " << stripSize << " strips\n";

    cout << "Lower Bound: " << LB << " strips\n";*/

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);

/*
    cout << "Strips FFD (scores):\n";
    for(i = 0; i < strip.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < strip[i].size(); ++j){
            cout << strip[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Strip" << setw(8) << "Width" << setw(12) << "Residual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << setw(9) <<  stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;*/
}

void packStripsFFDSmallest(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numBox, int maxBoxWidth,
                           int maxStripWidth, double totalBoxWidth, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
                           vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, k, l;
    int count = 0;
    vector<int> scoreInc;
    vector<int> checked(numScores, 0);

    for(i = 0; i < numScores; ++i){
        scoreInc.push_back(i);
    }

    for(i = 0; i < scoreInc.size() - 1; ++i){
        for(j = i+1; j < scoreInc.size(); ++j){
            if(allScores[scoreInc[i]] == allScores[scoreInc[j]]){
                if(allScores[mates[scoreInc[i]]] < allScores[mates[scoreInc[j]]]){
                    swap(scoreInc[i], scoreInc[j]);
                }
            }
            else{
                break;
            }
        }
    }

    /*cout << "Scores Increase:\n";
    for(auto number : scoreInc){
        cout << number << " ";
    }
    cout << endl;*/


    j = 0;
    while (count < numScores){
        for(i = 0; i < scoreInc.size(); ++i){
            if(checked[scoreInc[i]] == 1){
                continue;
            }
            if(strip[j].empty()){
                strip[j].push_back(scoreInc[i]);
                strip[j].push_back(mates[scoreInc[i]]);
                stripSum[j] += boxWidths[scoreInc[i]][mates[scoreInc[i]]];
                checked[scoreInc[i]] = 1;
                checked[mates[scoreInc[i]]] = 1;
                count += 2;
                i = -1;
                continue;
            }
            else if(!strip[j].empty()){
                if(adjMatrix[strip[j].back()][scoreInc[i]] == 1 && stripSum[j] + boxWidths[scoreInc[i]][mates[scoreInc[i]]] <= maxStripWidth){
                    strip[j].push_back(scoreInc[i]);
                    strip[j].push_back(mates[scoreInc[i]]);
                    stripSum[j] += boxWidths[scoreInc[i]][mates[scoreInc[i]]];
                    checked[scoreInc[i]] = 1;
                    checked[mates[scoreInc[i]]] = 1;
                    count += 2;
                    i = -1;
                    continue;
                }
            }
        }
        ++j;
    }

    k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalBoxWidth);

    /*cout << "After FFD: " << stripSize << " strips\n";

    cout << "Lower Bound: " << LB << " strips\n";*/

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);


    //cout << "Strips FFD (scores):\n";
    /*for(i = 0; i < strip.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < strip[i].size(); ++j){
            cout << strip[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

    //cout << "Strip" << setw(8) << "Width" << setw(12) << "Residual\n";
    /*for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << setw(9) <<  stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;*/


}

void packStripsFFDExact(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numBox, int maxBoxWidth,
                        int maxStripWidth, double totalBoxWidth, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
                        vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, mini, k, l;
    int min = 0;
    int max = maxBoxWidth;
    int feasible;
    vector<int> boxDecrease;
    vector<int> checked(numScores, 0);


    while(boxDecrease.size() < numBox) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] <= max) {
                min = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxDecrease.push_back(mini);
        checked[mini] = 1;
        checked[mates[mini]] = 1;
        max = min;
        min = 0;
    }

    /*cout << "Box decrease:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;*/

    strip[0].push_back(boxDecrease[0]);
    strip[0].push_back(mates[boxDecrease[0]]);
    stripSum[0] += boxWidths[boxDecrease[0]][mates[boxDecrease[0]]];

    //cout << "hello\n";
    for(j = 1; j < boxDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[boxDecrease[j]][mates[boxDecrease[j]]] <= maxStripWidth){
                    feasible = 0;
                    MBAHRA(i, j, feasible, allScores, mates, adjMatrix, boxWidths, boxDecrease, stripSum, strip);
                    if(feasible == 1){
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(boxDecrease[j]);
                strip[i].push_back(mates[boxDecrease[j]]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                break;
            }
        }
    }




    k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    l = stripSum.size() -1;
    while(stripSum[l] == 0){
        stripSum.pop_back();
        --l;
    }

    int stripSize = strip.size();
    int LB = lowerBound(maxStripWidth, totalBoxWidth);

    /*cout << "After FFD: " << stripSize << " strips\n";

    cout << "Lower Bound: " << LB << " strips\n";*/

    optimality(opt, opt90, opt80, opt70,opt60, opt50, optLow, stripSize, LB);


    //cout << "Strips FFD (scores):\n";
    /*for(i = 0; i < strip.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < strip[i].size(); ++j){
            cout << strip[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

    //cout << "Strip" << setw(8) << "Width" << setw(12) << "Residual\n";
    /*for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << setw(9) <<  stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;*/


}

void MBAHRA(int i1, int j1, int &feasible, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
            vector<vector<int> > &boxWidths, vector<int> &boxDecrease, vector<int> &stripSum, vector<vector<int> > &strip){

    feasible = 0;
    int i, j, k;
    int threshold = 70;
    int vacant = 999;
    vector<int> scoresX;
    vector<int> orderX;
    vector<int> originalX;
    vector<int> finalX;

    //region Creating scoresX vector
    for(k = 0; k < strip[i1].size(); ++k){
        scoresX.push_back(allScores[strip[i1][k]]);
        originalX.push_back(strip[i1][k]);
    }
    scoresX.push_back(allScores[boxDecrease[j1]]);
    originalX.push_back(boxDecrease[j1]);
    scoresX.push_back(allScores[mates[boxDecrease[j1]]]);
    originalX.push_back(mates[boxDecrease[j1]]);
    scoresX.push_back(70);
    scoresX.push_back(70);

    //endregion

    //region Variables
    int nScoresX = scoresX.size();
    int nBoxX = scoresX.size() / 2;
    int nCompX = (nBoxX + (nBoxX % 2)) / 2;
    vector<int> invOrderX(nScoresX);
    vector<vector<int> > adjMatX(nScoresX, vector<int>(nScoresX, 0));
    vector<int> matesX(nScoresX, 0);

    //MTGMA
    int lastMatchX = vacant;
    int mateMatchX = vacant;
    int vacantFlagX = 0;
    int matchSizeX = 0;
    vector<int> cycleVertexX(nScoresX, 1);
    vector<int> matchListX(nScoresX, vacant);

    //MIS
    int numCyclesX = 0;
    int smallestVertexX;
    int currentVertexX;
    vector<vector<int> > mateInducedX;
    vector<int> lengthMateInducedX;
    vector<int> tempMISX;
    vector<int> checkedX(nScoresX, 0);

    //FCA
    int qstarX;
    int numEdgesX;
    vector<vector<int> > SX(nCompX, vector<int>(nCompX, 0));
    vector<vector<int> > TX;
    vector<int> edgeX;
    vector<int> tX;

    //PatchGraph
    int q, u, v, saveX, SSumX, SqIntSX;
    int fullX = vacant;
    vector<int> SSetX;
    vector<int> tempPGX;
    vector<int> patchVertexX(nScoresX, vacant);
    vector<vector<int> > TpatchX;

    //MakePath
    vector<int> completePathX;
    vector<int> fullCycleX;
    //endregion

    if(feasible == 0) {
        //region Creating Instance

        for (k = 0; k < nScoresX; ++k) {
            orderX.push_back(k);
        }

        /*cout << "Scores:\n";
        for (k = 0; k < nScoresX; ++k) {
            cout << scoresX[k] << " ";
        }
        cout << endl << endl;*/

        /*cout << "Order:\n";
        for (k = 0; k < nScoresX; ++k) {
            cout << orderX[k] << " ";
        }
        cout << endl << endl;*/

        for (i = 1; i < nScoresX; ++i) {
            for (j = i - 1; j >= 0; --j) {
                if (scoresX[i] < scoresX[orderX[j]]) {
                    orderX[j + 1] = orderX[j];
                    orderX[j] = i;
                }
            }
        }

        /*cout << "Order:\n";
        for (k = 0; k < nScoresX; ++k) {
            cout << orderX[k] << " ";
        }
        cout << endl << endl;*/


        for (k = 0; k < nScoresX; ++k) {
            invOrderX[orderX[k]] = k;
        }

        /*cout << "Inverse Order:\n";
        for (k = 0; k < nScoresX; ++k) {
            cout << invOrderX[k] << " ";
        }
        cout << endl << endl;*/

        for (i = 0; i < nScoresX - 1; i += 2) {
            adjMatX[invOrderX[i]][invOrderX[i + 1]] = 2;
            adjMatX[invOrderX[i + 1]][invOrderX[i]] = 2;
        }

        sort(scoresX.begin(), scoresX.end());

        for (i = 0; i < scoresX.size() - 1; ++i) {
            for (j = i + 1; j < scoresX.size(); ++j) {
                if (scoresX[i] + scoresX[j] >= threshold && adjMatX[i][j] != 2) {
                    adjMatX[i][j] = 1;
                    adjMatX[j][i] = 1;
                }
            }

        }

        /*cout << "Adjacency Matrix\n";
        for (i = 0; i < nScoresX; ++i) {
            for (j = 0; j < nScoresX; ++j) {
                cout << adjMatX[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;*/

        for (i = 0; i < nScoresX; ++i) {
            for (j = 0; j < nScoresX; ++j) {
                if (adjMatX[i][j] == 2) {
                    matesX[i] = j;
                    break;
                }
            }
        }
        /*cout << "Mates Vector:\n";
        for (i = 0; i < matesX.size(); ++i) {
            cout << matesX[i] << " ";
        }
        cout << endl << endl;*/
        //endregion

        //region MTGMA
        /**MTGMA**/
        for (i = 0; i < nScoresX; ++i) { //check all vertices
            vacantFlagX = 0;
            if (matchListX[i] == vacant) { //if vertex has not yet been matched
                for (j = nScoresX - 1; j >
                                       i; --j) { //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
                    if (adjMatX[i][j] == 1 && matchListX[j] ==
                                              vacant) { //if vertices i and j are adjacent, and if vertex j has not yet been matched
                        matchListX[i] = j;
                        matchListX[j] = i;
                        lastMatchX = i;
                        ++matchSizeX;
                        if (vacantFlagX ==
                            1) { //delete edge for FCA if matching was not with highest vertex due to the highest vertex being its mate
                            cycleVertexX[i] = vacant;
                            cycleVertexX[j] = vacant;
                        }
                        break;
                    }
                    else if (adjMatX[i][j] == 2 && matchListX[j] == vacant) { //if potential match == mate
                        vacantFlagX = 1;
                    }
                }//end for j
                if (matchListX[i] == vacant) { //if vertex has still not been matched
                    for (k = 0; k < nScoresX - 2; ++k) {
                        if (adjMatX[i][k] == 2) { //if vertex i and vertex k are mates
                            mateMatchX = k;
                            break;
                        }
                    }
                    if ((scoresX[i] + scoresX[mateMatchX] >= threshold) //match with mate?
                        && (matchListX[mateMatchX] == vacant) //is mate unmatched?
                        && (lastMatchX != vacant) //has the previous vertex been matched?
                        && (mateMatchX >
                            i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
                        && (scoresX[lastMatchX] + scoresX[mateMatchX] >=
                            threshold)) { //can mate be matched with last matched vertex?
                        // if so, then swap mates
                        matchListX[i] = matchListX[lastMatchX];
                        matchListX[lastMatchX] = mateMatchX;
                        matchListX[mateMatchX] = lastMatchX;
                        matchListX[matchListX[i]] = i;
                        cycleVertexX[lastMatchX] = vacant; //edge from mate swap will not count for FCA
                        cycleVertexX[mateMatchX] = vacant; //edge from mate swap will not count for FCA
                        lastMatchX = i;
                        ++matchSizeX;
                    }
                }//end if matchList == vacant
            }//end if matchList[i] == i
        }//end for i


        /*cout << "Cycle Vertex vector after MTGMA:\n";
        for(i = 0; i < cycleVertexX.size(); ++i){
            cout << cycleVertexX[i] << " ";
        }
        cout << endl << endl;

        cout << "Matching List:\n";
        for(i = 0; i < matchListX.size(); ++i){
            cout << matchListX[i] << " ";
        }
        cout << endl << endl;*/

        //endregion
        if (matchSizeX < nBoxX) {
            //NOT ENOUGH MATCHING EDGES
            feasible = 0;
            goto End;
        }

        //region MIS
        /**MIS**/
        //find the smallest vertex not yet checked for mate-induced structure - start with this vertex
        for (i = 0; i < nScoresX; ++i) {
            if (checkedX[i] == 0) {
                smallestVertexX = i;
                break;
            }
        }

        //Building the mate-induced structure
        do {
            currentVertexX = smallestVertexX;
            do {
                tempMISX.push_back(currentVertexX);
                checkedX[currentVertexX] = 1;
                tempMISX.push_back(matesX[currentVertexX]);
                checkedX[matesX[currentVertexX]] = 1;
                currentVertexX = matchListX[matesX[currentVertexX]];
            } while (currentVertexX != smallestVertexX);

            mateInducedX.push_back(tempMISX);
            tempMISX.clear();

            for (i = 0; i < nScoresX; ++i) {
                if (checkedX[i] == 0) {
                    smallestVertexX = i;
                    break;
                }
            }


        } while (smallestVertexX != currentVertexX);

        tempMISX.clear(); //clear cycle vector again for next instance

        numCyclesX = mateInducedX.size(); //number of cycles in the mate-induced structure

        /*cout << "Mate-Induced Structure:\n";
        for(i = 0; i < mateInducedX.size(); ++i){
            for(j = 0; j < mateInducedX[i].size(); ++j){
                cout << mateInducedX[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;*/
        //cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

        for (i = 0; i < mateInducedX.size(); ++i) {
            lengthMateInducedX.push_back(mateInducedX[i].size());
        }

        //endregion
        if (lengthMateInducedX[0] == nScoresX) {
            for (j = 0; j < mateInducedX[0].size(); ++j) {
                fullCycleX.push_back(mateInducedX[0][j]);
            }
            //region MAKEPATH
            /**MAKEPATH FUNCTION**/
            for (i = 0; i < fullCycleX.size() - 1; ++i) {
                if ((fullCycleX[i] == nScoresX - 1 && fullCycleX[i + 1] == nScoresX - 2) ||
                    (fullCycleX[i] == nScoresX - 2 && fullCycleX[i + 1] == nScoresX - 1)) {
                    if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                        for (j = 2; j < fullCycleX.size(); ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }

                    else if (i == fullCycleX.size() -
                                  2) { //if the dominating vertices are at the end of the fullCycle vector
                        for (j = 0; j < fullCycleX.size() - 2; ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }
                    else { //if the dominating vertices are in the middle of the fullCycle vector
                        for (j = i + 2; j < fullCycleX.size(); ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;

                    }

                }
            }

            for (i = 0; i < completePathX.size(); ++i) {
                finalX.push_back(originalX[orderX[completePathX[i]]]);
            }
            feasible = 1;
            goto End;

            //END MAKE PATH FUNCTION
            //endregion
        }

        //region FCA
        /**FCA**/
        //create list cycleVertex that contains for each vertex the cycle that each edge belongs to
        for (i = 0; i < mateInducedX.size(); ++i) {
            for (j = 0; j < mateInducedX[i].size(); ++j) {
                if (cycleVertexX[mateInducedX[i][j]] != vacant) { //if edge is not deleted for FCA
                    cycleVertexX[mateInducedX[i][j]] = i;
                }
            }
        }

        /*cout << "Cycle Vertex:\n";
        for (i = 0; i < cycleVertexX.size(); ++i) {
            cout << cycleVertexX[i] << " ";
        }
        cout << endl << endl;*/

        //create list of edges without empty edges (those generated by mate swap)
        for (i = 0; i < matchSizeX; ++i) {
            while (cycleVertexX[i] == vacant) {
                ++i;
            }
            edgeX.push_back(i);
        }
        numEdgesX = edgeX.size();

        /*cout << "Edges vector:\n";
        for (i = 0; i < edgeX.size(); ++i) {
            cout << edgeX[i] << " ";
        }
        cout << endl << endl;*/

        //cout << "Number of Edges: " << numEdges << endl;
        //FCA Algorithm
        qstarX = -1;
        k = 0; //edge from matching that is under consideration

        do {
            while (k < numEdgesX - 2 && (adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] != 1 ||
                                         cycleVertexX[edgeX[k]] == cycleVertexX[edgeX[k + 1]])) {
                ++k;
            }
            if (adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] == 1 &&
                cycleVertexX[edgeX[k]] != cycleVertexX[edgeX[k + 1]]) {
                ++qstarX;
                tX.push_back(edgeX[k]);
                SX[qstarX][cycleVertexX[edgeX[k]]] = 1;
                while (k < numEdgesX - 1 && adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] == 1 &&
                       SX[qstarX][cycleVertexX[edgeX[k + 1]]] == 0) { //add more edges to current T-cycle
                    ++k;
                    tX.push_back(edgeX[k]);
                    SX[qstarX][cycleVertexX[edgeX[k]]] = 1;
                }
                TX.push_back(tX);
                tX.clear();
            } // end if
            ++k;
        } while (k < numEdgesX - 1);

        tX.clear();

        /*cout << "T matrix:\n";
        for(i = 0; i < TX.size(); ++i){
            for(j = 0; j < TX[i].size(); ++j){
                cout << TX[i][j] << "  ";
            }
            cout << endl;
        }
        cout << endl << endl;*/

        /*cout << "S Matrix:\n";
        for(i = 0; i < SX.size(); ++i){
            for(j = 0; j < SX[i].size(); ++j){
                cout << SX[i][j] << "  ";
            }
            cout << endl;
        }
        cout << endl;*/

        //cout << "qstar: " << qstarX << endl << endl;

        //endregion
        if (qstarX == -1) {
            //INFEASIBLE, EXIT
            feasible = 0;
            goto End;
        }

        //region PATCHGRAPH
        /**PATCHGRAPH**/
        vector<int> QSetX(qstarX, 0);
        vector<int> patchCycleX(qstarX, vacant);
        for (i = 0; i < TX.size(); ++i) {
            if (TX[i].size() == numCyclesX) {
                fullX = i;
                break;
            }
        }

        if (fullX != vacant) {
            saveX = 0;
            //cout << "Full: " << full << endl;
            for (v = 0; v < TX[fullX].size(); ++v) {
                for (j = 0; j < mateInducedX[cycleVertexX[TX[fullX][v]]].size(); ++j) {
                    if (mateInducedX[cycleVertexX[TX[fullX][v]]][j] == matchListX[TX[fullX][v]]) {
                        saveX = j;
                        break;
                    }
                }
                //region oneTCyclePatch
                /****************************oneTCyclePatch algorithm: ***********************************/
                //CASE ONE: if element matchList[T[full][v]] is before element T[full][v] in the mateInduced cycle
                //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save + 1" is T[full][v]
                if (mateInducedX[cycleVertexX[TX[fullX][v]]][saveX + 1] == TX[fullX][v]) {
                    for (i = saveX + 1; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                    for (i = mateInducedX[cycleVertexX[TX[fullX][v]]].size();
                         i-- > saveX + 1;) { //from end of cycle to element at position save+1
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                }

                    //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
                    //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
                else if (mateInducedX[cycleVertexX[TX[fullX][v]]][saveX - 1] == TX[fullX][v]) {
                    for (i = saveX; i < mateInducedX[cycleVertexX[TX[fullX][v]]].size(); ++i) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                    for (i = 0; i < saveX; ++i) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                }

                    //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
                    //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
                else if (saveX == 0 &&
                         mateInducedX[cycleVertexX[TX[fullX][v]]][mateInducedX[cycleVertexX[TX[fullX][v]]].size() -
                                                                  1] == TX[fullX][v]) {
                    for (i = 0; i < mateInducedX[cycleVertexX[TX[fullX][v]]].size(); ++i) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                }

                    //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
                    //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
                else if (saveX == mateInducedX[cycleVertexX[TX[fullX][v]]].size() - 1 &&
                         mateInducedX[cycleVertexX[TX[fullX][v]]][0] == TX[fullX][v]) {
                    for (i = mateInducedX[cycleVertexX[TX[fullX][v]]].size(); i-- > 0;) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                    }
                }
                //END ONETCYCLEPATCH FUNCTION
                //endregion
            }

            //region MAKEPATH
            /**MAKEPATH FUNCTION**/
            for (i = 0; i < fullCycleX.size() - 1; ++i) {
                if ((fullCycleX[i] == nScoresX - 1 && fullCycleX[i + 1] == nScoresX - 2) ||
                    (fullCycleX[i] == nScoresX - 2 && fullCycleX[i + 1] == nScoresX - 1)) {
                    if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                        for (j = 2; j < fullCycleX.size(); ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }

                    else if (i == fullCycleX.size() -
                                  2) { //if the dominating vertices are at the end of the fullCycle vector
                        for (j = 0; j < fullCycleX.size() - 2; ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }
                    else { //if the dominating vertices are in the middle of the fullCycle vector
                        for (j = i + 2; j < fullCycleX.size(); ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;

                    }

                }
            }

            for (i = 0; i < completePathX.size(); ++i) {
                finalX.push_back(originalX[orderX[completePathX[i]]]);
            }
            feasible = 1;
            goto End;
            //END MAKE PATH FUNCTION
            //endregion

        }

        else {
            q = 0; //Start with first Tq-cycle
            QSetX[0] = 1;
            SSumX = 0; //number of MIS-cycles that have been included
            for (i = 0; i < numCyclesX; ++i) {
                SSetX.push_back(SX[q][i]); // ==1 if MIS cycle i has been included
            }
            for (i = 0; i < numCyclesX; ++i) {
                SSumX = SSumX + SSetX[i];
            }

            if (SSumX >= 1) {
                patchCycleX[q] = 1;
            }

            //Start connectivity check
            while (q <= qstarX && SSumX < numCyclesX) {
                do {
                    ++q;
                    SqIntSX = vacant;
                    if (q <= qstarX) {
                        for (j = 0; j < numCyclesX; ++j) { //is there a j such that S[q][j] = 1 and SSet[j] = 1?
                            if (SX[q][j] == 1 && SSetX[j] == 1) {
                                SqIntSX = 1;
                                //break here? no need to check all other j indices once one has been found such that S[q][j] =1 and SSet[j] = 1
                            }
                        }
                    }
                } while (q < qstarX + 1 && (QSetX[q] == 1 || SqIntSX == vacant));

                if (q <= qstarX) { //if Tq-cyce for enlargement has been found
                    for (i = 0; i < numCyclesX; ++i) {
                        if (SSetX[i] == 0 && SX[q][i] == 1) {
                            SSetX[i] = 1;
                            ++SSumX;
                            patchCycleX[q] = 1;
                        }
                    }
                    QSetX[q] = 1;
                    q = 0;
                }
            }//end while


            //If patching graph is connected, then instance is feasible, else infeasible
            if (SSumX == numCyclesX) {
                for (i = 0; i < patchCycleX.size(); ++i) {
                    if (patchCycleX[i] == 1) {
                        for (j = 0; j < TX[i].size(); ++j) {
                            tempPGX.push_back(TX[i][j]);
                        }
                        TpatchX.push_back(tempPGX);
                        tempPGX.clear();
                    }
                }
                tempPGX.clear();

                for (i = 0; i < TpatchX.size(); ++i) {
                    for (j = 0; j < TpatchX[i].size(); ++j) {
                        patchVertexX[TpatchX[i][j]] = i;
                        patchVertexX[matchListX[TpatchX[i][j]]] = i;
                    }
                }

                u = 0;
                v = 0;
                saveX = 0;
                for (j = 0; j < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++j) {
                    if (mateInducedX[cycleVertexX[TpatchX[u][v]]][j] == matchListX[TpatchX[u][v]]) {
                        saveX = j;
                        break;
                    }
                }
                //region multipleTCyclePatch
                /***************************** MULTIPLE T CYCLE PATCH FUNCTION ***********************/
                int x = 0;

                //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
                if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == TpatchX[u][v]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = 0; i < saveX; ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == matchListX[TpatchX[u][v]]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = 0; i < saveX; ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == matchListX[TpatchX[u][v]]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = saveX; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size();
                             i-- > saveX + 1;) { //from end of cycle to element at position save+1
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == TpatchX[u][v]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = saveX; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size();
                             i-- > saveX + 1;) { //from end of cycle to element at position save+1
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                    //the last element in the cycle is Tpatch[u][v]
                else if (saveX == 0 &&
                         mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size() -
                                                                   1] == TpatchX[u][v]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                    //the first element in the cycle is Tpatch[u][v]
                else if (saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1 &&
                         mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == TpatchX[u][v]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1; i-- > 0;) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                    //and the last element in the cycle is matchList[Tpatch[u][v]]
                else if (saveX == 0 &&
                         mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size() -
                                                                   1] == matchListX[TpatchX[u][v]]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                    //the first element in the cycle is matchList[Tpatch[u][v]]
                else if (saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1 &&
                         mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == matchListX[TpatchX[u][v]]) {
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1; i-- > 0;) {
                        if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }
                //END MULTIPLE T CYCLE PATCH FUNCTION
                //endregion

                while (fullCycleX.size() < nScoresX) {
                    saveX = 0;
                    for (i = 0; i < TpatchX[u].size(); ++i) {
                        if (TpatchX[u][i] == fullCycleX.back()) {
                            if (i == TpatchX[u].size() - 1) {
                                v = 0;
                            }
                            else {
                                v = ++i;
                            }
                            for (j = 0; j < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++j) {
                                if (mateInducedX[cycleVertexX[TpatchX[u][v]]][j] == matchListX[TpatchX[u][v]]) {
                                    saveX = j;
                                    break;
                                }
                            }
                            break;
                        }
                        else if (matchListX[TpatchX[u][i]] == fullCycleX.back()) {
                            if (i == 0) {
                                v = TpatchX[u].size() - 1;
                            }
                            else {
                                v = --i;
                            }
                            for (j = 0; j < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++j) {
                                if (mateInducedX[cycleVertexX[TpatchX[u][v]]][j] == TpatchX[u][v]) {
                                    saveX = j;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    //region multipleTCyclePatch
                    /***************************** MULTIPLE T CYCLE PATCH FUNCTION ***********************/
                    int x = 0;

                    //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
                    if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == TpatchX[u][v]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = 0; i < saveX; ++i) {
                                if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                }
                                else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                    u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                    else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == matchListX[TpatchX[u][v]]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = 0; i < saveX; ++i) {
                                if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                }
                                else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                    u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                    else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == matchListX[TpatchX[u][v]]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = saveX; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size();
                                 i-- > saveX + 1;) { //from end of cycle to element at position save+1
                                if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                }
                                else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                    u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                    else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == TpatchX[u][v]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = saveX; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size();
                                 i-- > saveX + 1;) { //from end of cycle to element at position save+1
                                if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                }
                                else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                    u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                        //the last element in the cycle is Tpatch[u][v]
                    else if (saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][
                                                   mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1] ==
                                           TpatchX[u][v]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                        //the first element in the cycle is Tpatch[u][v]
                    else if (saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1 &&
                             mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == TpatchX[u][v]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1; i-- > 0;) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                        //and the last element in the cycle is matchList[Tpatch[u][v]]
                    else if (saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][
                                                   mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1] ==
                                           matchListX[TpatchX[u][v]]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                        //the first element in the cycle is matchList[Tpatch[u][v]]
                    else if (saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1 &&
                             mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == matchListX[TpatchX[u][v]]) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size() - 1; i-- > 0;) {
                            if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if (patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                    //END MULTIPLE T CYCLE PATCH FUNCTION
                    //endregion

                }

                //region MAKEPATH
                /**MAKEPATH FUNCTION**/
                for (i = 0; i < fullCycleX.size() - 1; ++i) {
                    if ((fullCycleX[i] == nScoresX - 1 && fullCycleX[i + 1] == nScoresX - 2) ||
                        (fullCycleX[i] == nScoresX - 2 && fullCycleX[i + 1] == nScoresX - 1)) {
                        if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                            for (j = 2; j < fullCycleX.size(); ++j) {
                                completePathX.push_back(fullCycleX[j]);
                            }
                            break;
                        }

                        else if (i == fullCycleX.size() -
                                      2) { //if the dominating vertices are at the end of the fullCycle vector
                            for (j = 0; j < fullCycleX.size() - 2; ++j) {
                                completePathX.push_back(fullCycleX[j]);
                            }
                            break;
                        }
                        else { //if the dominating vertices are in the middle of the fullCycle vector
                            for (j = i + 2; j < fullCycleX.size(); ++j) {
                                completePathX.push_back(fullCycleX[j]);
                            }
                            for (j = 0; j < i; ++j) {
                                completePathX.push_back(fullCycleX[j]);
                            }
                            break;

                        }

                    }
                }

                for (i = 0; i < completePathX.size(); ++i) {
                    finalX.push_back(originalX[orderX[completePathX[i]]]);
                }
                feasible = 1;
                goto End;
                //END MAKE PATH FUNCTION
                //endregion

            }
            else if (SSumX < numCyclesX) {
                feasible = 0;
                goto End;

            }
            else {
                feasible = 0;
                goto End;
            }

        }

        //endregion
    }

    //END MBAHRA STRIPX


    //region End
    End:
    if(feasible == 1){
        stripSum[i1] += boxWidths[boxDecrease[j1]][mates[boxDecrease[j1]]];
        strip[i1].swap(finalX);
    }
    //endregion

}//end void MBAHRA



























































































