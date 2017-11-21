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

int lowerBound(double totalBoxWidth, int maxStripWidth){
    int lBound;

    lBound = ceil(totalBoxWidth/maxStripWidth);
    return lBound;
}

void packStripsFFD(int &totalCost, int numBox, int maxBoxWidth, int maxStripWidth, double totalBoxWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<int> &stripNumBoxes, vector<vector<int> > &strip, vector<vector<int> > &stripWidth){

    int i, j, mini, k, l, m;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> boxDecrease;

    while(boxDecrease.size() < numBox) {
        for (i = 0; i < mates.size(); ++i) {
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] < max) { //what happens if two boxes have the same width?
                min = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxDecrease.push_back(mini);
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
    stripWidth[0].push_back(boxWidths[boxDecrease[0]][mates[boxDecrease[0]]]);


    for(j = 1; j < boxDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[boxDecrease[j]][mates[boxDecrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][boxDecrease[j]] == 1){
                        strip[i].push_back(boxDecrease[j]);
                        strip[i].push_back(mates[boxDecrease[j]]);
                        stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][mates[boxDecrease[j]]] == 1){
                        strip[i].push_back(mates[boxDecrease[j]]);
                        strip[i].push_back(boxDecrease[j]);
                        stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(boxDecrease[j]);
                strip[i].push_back(mates[boxDecrease[j]]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                break;
            }
        }
    }

    k = strip.size() - 1;
    while(strip[k].empty()){
        strip.pop_back();
        --k;
    }

    l = stripWidth.size() -1;
    while(stripWidth[l].empty()){
        stripWidth.pop_back();
        --l;
    }

    m = stripSum.size() -1;
    while(stripSum[m] == 0){
        stripSum.pop_back();
        --m;
    }

    for(i = 0; i < strip.size(); ++i){
        stripNumBoxes.push_back(strip[i].size() / 2);
        ++numStrips;
    }


    cout << "FFD: " << numStrips << " strips\n";

    cout << "Lower Bound: " << lowerBound(totalBoxWidth, maxStripWidth) << " strips\n---------------\n";


    cout << "Strips FFD (scores):\n";
    for(i = 0; i < strip.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < strip[i].size(); ++j){
            cout << strip[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    /*cout << "Strips FFD (boxWidths):\n";
    for(i = 0; i < stripWidth.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < stripWidth[i].size(); ++j){
            cout << stripWidth[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

    cout << "Strip" << setw(8) << "#Boxes" << setw(8) << "Width" << setw(12) << "Residual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << setw(9) << stripNumBoxes[i] << setw(10) <<  stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    int costFFD = initCost(totalCost, maxStripWidth, stripSum);
    cout << "\nInitial Cost (cost from FFD): " << costFFD << "\n-------------------------------\n\n";


}

int initCost(int &totalCost, int maxStripWidth, vector<int> &stripSum){

    int i;
    totalCost = 0;

    for(i = 0; i < stripSum.size(); ++i){
        totalCost += pow((maxStripWidth - stripSum[i]), 2);
    }

    return totalCost;

}

void localSearch(int maxStripWidth, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<int> &stripSumX, vector<int> &stripSumY, vector<vector<int> > &strip, vector<vector<int> > &stripX, vector<vector<int> > &stripY){

    int a, b, c, d, i, j, k, l, half, pairSizeX, pairSizeY;

    /* Initially:
     * We split the strips into two sets of strips, stripX and stripY
     * If the total number of strips is even, then we put half the strips in each set
     * If the total number of strips is odd, then stripX will have one more strip than stripY
     * We also have to create two new vectors to hold the total strip widths, stripSumX and stripSumY
     */

    /*if(strip.size() % 2 == 1){
        half = (strip.size() / 2) + 1;
    }
    else{
        half = strip.size() / 2;
    }

    for(i = 0; i < half; ++i){
        stripX.push_back(strip[i]);
    }
    for(i = half; i < strip.size(); ++i){
        stripY.push_back(strip[i]);
    }
    for(i = 0; i < half; ++i){
        stripSumX.push_back(stripSum[i]);
    }
    for(i = half; i < strip.size(); ++i){
        stripSumY.push_back(stripSum[i]);
    }*/

    if(strip.size() % 2 == 1){
        for(i = 0; i < strip.size() - 1; i+=2){
            stripX.push_back(strip[i]);
            stripSumX.push_back(stripSum[i]);
            stripY.push_back(strip[i+1]);
            stripSumY.push_back(stripSum[i+1]);
        }
        stripX.push_back(strip[strip.size()-1]);
        stripSumX.push_back(stripSum[stripSum.size()-1]);
    }
    else{
        for(i = 0; i < strip.size(); i+=2){
            stripX.push_back(strip[i]);
            stripSumX.push_back(stripSum[i]);
            stripY.push_back(strip[i+1]);
            stripSumY.push_back(stripSum[i+1]);
        }
    }


    cout << "stripX:\n";
    for(i = 0; i < stripX.size(); ++i){
        for(j = 0; j < stripX[i].size(); ++j){
            cout << stripX[i][j] << " ";
        }
        cout << endl;
    }
    cout << "stripSumX: ";
    for(i = 0; i < stripSumX.size(); ++i){
        cout << stripSumX[i] << " ";
    }
    cout << endl << endl;

    cout << "stripY:\n";
    for(i = 0; i < stripY.size(); ++i){
        for(j = 0; j < stripY[i].size(); ++j){
            cout << stripY[i][j] << " ";
        }
        cout << endl;
    }
    cout << "stripSumY: ";
    for(i = 0; i < stripSumY.size(); ++i){
        cout << stripSumY[i] << " ";
    }
    cout << endl << endl;


    /*SWAPPING A PAIR OF BOXES FROM EACH SET*/
    for(i = 0; i < stripX.size(); ++i){
        if(stripX[i].size() >= 4){
            for(a = 0; a < stripX[i].size()-3; a+=2){
                for(b = a+2; b < stripX[i].size()-1; b+=2){
                    pairSizeX = boxWidths[stripX[i][a]][stripX[i][a+1]] + boxWidths[stripX[i][b]][stripX[i][b+1]];
                    for(j = 0; j < stripY.size(); ++j){
                        if(stripY[j].size() >= 4){
                            for(c = 0; c < stripY[j].size()-3; c+=2){
                                for(d = c+2; d < stripY[j].size()-1; d+=2){
                                    pairSizeY = boxWidths[stripY[j][c]][stripY[j][c+1]] + boxWidths[stripY[j][d]][stripY[j][d+1]];
                                    if(pairSizeX < pairSizeY && stripSumX[i] - pairSizeX + pairSizeY <= maxStripWidth){
                                        //MBAHRA
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "NO PAIR PAIR SWAP PERFORMED.\n";


    /*SWAPPING A PAIR OF BOXES FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        if(stripX[i].size() >= 4){ //If there are at least 2 boxes on stripX[i]
            for(a = 0; a < stripX[i].size()-3; a+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                for(b = a+2; b < stripX[i].size()-1; b+=2){ //Starting from the first score on the second box until the first score on the last box
                    pairSizeX = boxWidths[stripX[i][a]][stripX[i][a+1]] + boxWidths[stripX[i][b]][stripX[i][b+1]]; //Sum box widths
                    //Check if there exists a box on a strip in set stripY whose width is larger than pairSizeX
                    for(j = 0; j < stripY.size(); ++j){
                        //Go through each box on stripY[j]
                        for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box unil the first score on the last box
                            //Check if pairSizeX < width of box in stripY, and that box can fit onto strip
                            if(pairSizeX <= boxWidths[stripY[j][c]][stripY[j][c+1]] && stripSumX[i] - pairSizeX + boxWidths[stripY[j][c]][stripY[j][c+1]] <= maxStripWidth){
                                //Check if boxes meet mssc on opposite strips
                                //If so, perform swap, update stripSums, and move onto SinSin
                                //If not, continue
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "NO PAIR SIN SWAP PERFORMED.\n";

    /*SWAPPING ONE BOX FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        for(a = 0; a < stripX[i].size()-1; a+=2){ // Starting from the first score on the first box until the first score on the last box
            for(j = 0; j < stripY.size(); ++j){ // For each strip in the set stripY
                for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
                    //Check if boxwidth[a] < boxWidth[c] and that box can fit on strip
                    if(boxWidths[stripX[i][a]][stripX[i][a+1]] < boxWidths[stripY[j][c]][stripY[j][c+1]]
                       && stripSumX[i] - boxWidths[stripX[i][a]][stripX[i][a+1]] + boxWidths[stripY[j][c]][stripY[j][c+1]] <= maxStripWidth){
                        //Check if boxes meet mssc on opposite strips
                        //If so, perform swap, update stripSums, and move onto SinSin
                        //If not, continue
                    }
                }
            }
        }
    }
    cout << "NO SINGLE SINGLE SWAP PERFORMED.\n";

    /*MOVING ONE BOX FROM SET STRIPY TO SET STRIPX*/
    for(j = 0; j < stripY.size(); ++j){ //For each strip in the set stripY
        for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
            for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
                if(stripSumX[i] + boxWidths[stripY[j][c]][stripY[j][c+1]] <= maxStripWidth){
                    //Check if box c meets mssc on strip in set stripX
                    //If so, perform move, update stripSums, go back to PairPair
                    //If not, continue onto next stripX[i]
                }
            }
        }
    }
    cout << "NO SINGLE MOVE PERFORMED.\n";

    /*SWAPPING A PAIR OF BOXES FROM EACH SET*/
    /*for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        if(stripX[i].size() >= 4){ //If there are at least 2 boxes on stripX[i] (note that each element represents a score, so 4 elements = 2 boxes)
            //Go through each pair of boxes on stripX[i]
            for(a = 0; a < stripX[i].size()-3; a+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                for(b = a+2; b < stripX[i].size()-1; b+=2){ //Starting from the first score on the second box until the first score on the last box
                    pairSizeX = boxWidths[stripX[i][a]][stripX[i][a+1]] + boxWidths[stripX[i][b]][stripX[i][b+1]]; //Sum box widths
                    //Check if there exists a pair of boxes on a strip in set stripY that have a combined width larger than pairSizeX
                    for(j = 0; j < stripY.size(); ++j){ //For each strip in the set stripY
                        if(stripY[j].size() >= 4){ //If there are at least 2 boxes on stripY[j]
                            //Go through each pair of boxes on stripY[j]
                            for(c = 0; c < stripY[j].size()-3; c+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                                for(d = c+2; d < stripY[j].size()-1; d+=2){ //Starting from the first score on the second box until the first score on the last box
                                    pairSizeY = boxWidths[stripY[j][c]][stripY[j][c+1]] + boxWidths[stripY[j][d]][stripY[j][d+1]]; //Sum box widths
                                    //Check if pairSizeX < pairSizeY and that boxes can fit onto strip
                                    if(pairSizeX < pairSizeY && stripSumX[i] - pairSizeX + pairSizeY <= maxStripWidth){
                                        //Check if boxes meet mssc on opposite strips
                                        //If so, perform swap, update stripSums, and move onto PairSin
                                        //If not, continue
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "NO PAIR PAIR SWAP PERFORMED.\n";*/


}

void MBAHRA(int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores, vector<vector<int> > &stripX, vector<vector<int> > &stripY){

    int i, j, k, nScoresX, nBoxX, nCompX;
    int threshold = 70;
    int vacant = 999;
    vector<int> scoresX;
    vector<int> orderX;
    vector<int> originalX;
    vector<int> finalX;


    /**FOR NEW STRIPX PAIRPAIR**/
    for(k = 0; k < stripX[i1].size(); ++k){
        if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1){
            continue;
        }
        scoresX.push_back(allScores[stripX[i1][k]]);
        originalX.push_back(stripX[i1][k]);

    }
    scoresX.push_back(allScores[stripY[j1][c1]]);
    originalX.push_back(stripY[j1][c1]);
    scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
    originalX.push_back(stripY[j1][c1+1]);
    scoresX.push_back(allScores[stripY[j1][d1]]);
    originalX.push_back(stripY[j1][d1]);
    scoresX.push_back(allScores[stripY[j1][d1 + 1]]);
    originalX.push_back(stripY[j1][d1+1]);
    scoresX.push_back(70);
    scoresX.push_back(70);

    //region Initialisation
    nScoresX = scoresX.size();
    nBoxX = scoresX.size() /2;
    nCompX = (nBoxX + (nBoxX % 2)) / 2;
    vector<int> invOrderX(nScoresX);
    vector<vector<int> > adjMatX(nScoresX, vector<int>(nScoresX, 0));
    vector<int> matesX(nScoresX, 0);
    vector<int> completePathX;

    for(k = 0; k < nScoresX; ++k){
        orderX.push_back(k);
    }

    cout << "Scores:\n";
    for(k = 0; k < nScoresX; ++k){
        cout << scoresX[k] << " ";
    }
    cout << endl << endl;

    cout << "Order:\n";
    for(k = 0; k < nScoresX; ++k){
        cout << orderX[k] << " ";
    }
    cout << endl << endl;

    for(i = 1; i < nScoresX; ++i){
        for(j = i-1; j >= 0; --j){
            if(scoresX[i] < scoresX[orderX[j]]){
                orderX[j+1] = orderX[j];
                orderX[j] = i;
            }
        }
    }

    cout << "Order:\n";
    for(k = 0; k < nScoresX; ++k){
        cout << orderX[k] << " ";
    }
    cout << endl << endl;


    for(k = 0; k < nScoresX; ++k){
        invOrderX[orderX[k]] = k;
    }

    cout << "Inverse Order:\n";
    for(k = 0; k < nScoresX; ++k){
        cout << invOrderX[k] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < nScoresX - 1; i+=2){
        adjMatX[invOrderX[i]][invOrderX[i+1]] = 2;
        adjMatX[invOrderX[i+1]][invOrderX[i]] = 2;
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

    cout << "Adjacency Matrix\n";
    for(i = 0; i < nScoresX; ++i){
        for(j = 0; j < nScoresX; ++j){
            cout << adjMatX[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;


    for (i = 0; i < nScoresX; ++i) {
        for (j = 0; j < nScoresX; ++j) {
            if (adjMatX[i][j] == 2) {
                matesX[i] = j;
                break;
            }
        }
    }
    cout << "Mates Vector:\n";
    for(i = 0; i < matesX.size(); ++i){
        cout << matesX[i] << " ";
    }
    cout << endl << endl;
    //endregion


    //region MTGMA
    /**MTGMA**/
    int lastMatchX = vacant;
    int mateMatchX = vacant;
    int vacantFlagX = 0;
    int matchSizeX = 0;
    vector<int> cycleVertexX(nScoresX, 1);
    vector<int> matchListX(nScoresX, vacant);

    for (i = 0; i < nScoresX; ++i) { //check all vertices
        vacantFlagX = 0;
        if (matchListX[i] == vacant) { //if vertex has not yet been matched
            for (j = nScoresX - 1; j > i; --j) { //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
                if (adjMatX[i][j] == 1 && matchListX[j] == vacant) { //if vertices i and j are adjacent, and if vertex j has not yet been matched
                    matchListX[i] = j;
                    matchListX[j] = i;
                    lastMatchX = i;
                    ++matchSizeX;
                    if (vacantFlagX == 1) { //delete edge for FCA if matching was not with highest vertex due to the highest vertex being its mate
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
                    && (mateMatchX > i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
                    && (scoresX[lastMatchX] + scoresX[mateMatchX] >= threshold)) { //can mate be matched with last matched vertex?
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
    if(matchSizeX < nBoxX){
        //NOT ENOUGH MATCHING EDGES
        /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
    }

    //region MIS
    /**MIS**/
    int numCyclesX = 0;
    int smallestVertexX;
    int currentVertexX;
    vector<vector<int> > mateInducedX;
    vector<int> lengthMateInducedX;
    vector<int> tempMISX;
    vector<int> checkedX(nScoresX, 0);
    vector<int> fullCycleX;

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
    if(lengthMateInducedX[0] == nScoresX){
        for(j = 0; j < mateInducedX[0].size(); ++j){
            fullCycleX.push_back(mateInducedX[0][j]);
        }
        //region MAKEPATH
        /******************** MAKEPATH FUNCTION *************************************/
        for(i = 0; i < fullCycleX.size()-1; ++i){
            if((fullCycleX[i] == nScoresX - 1 && fullCycleX[i+1] == nScoresX - 2) || (fullCycleX[i] == nScoresX - 2 && fullCycleX[i+1] == nScoresX - 1)){
                if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                    for(j = 2; j < fullCycleX.size(); ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;
                }

                else if(i == fullCycleX.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                    for(j = 0; j < fullCycleX.size()-2; ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;
                }
                else{ //if the dominating vertices are in the middle of the fullCycle vector
                    for(j = i+2; j < fullCycleX.size(); ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    for(j = 0; j < i; ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;

                }

            }
        }

        for(i = 0; i < completePathX.size(); ++i){
            finalX.push_back(originalX[orderX[completePathX[i]]]);
        }

        //END MAKE PATH FUNCTION
        //endregion
        /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/
    }

    //region FCA
    /**FCA**/
    int qstarX;
    int numEdgesX;
    vector<vector<int> > SX(nCompX, vector<int>(nCompX, 0));
    vector<vector<int> > TX;
    vector<int> edgeX;
    vector<int> tX;

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
        while (k < numEdgesX - 2 && (adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] != 1 || cycleVertexX[edgeX[k]] == cycleVertexX[edgeX[k + 1]])) {
            ++k;
        }
        if (adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] == 1 && cycleVertexX[edgeX[k]] != cycleVertexX[edgeX[k + 1]]) {
            ++qstarX;
            tX.push_back(edgeX[k]);
            SX[qstarX][cycleVertexX[edgeX[k]]] = 1;
            while (k < numEdgesX - 1 && adjMatX[edgeX[k]][matchListX[edgeX[k + 1]]] == 1 && SX[qstarX][cycleVertexX[edgeX[k + 1]]] == 0) { //add more edges to current T-cycle
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
    if(qstarX == -1){
        //INFEASIBLE, EXIT
        /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
    }

    //region PATCHGRAPH
    /**PATCHGRAPH**/
    int q, u, v, saveX, SSumX, SqIntSX;
    int fullX = vacant;
    vector<int> QSetX(qstarX, 0);
    vector<int> SSetX;
    vector<int> patchCycleX(qstarX, vacant);
    vector<int> tempPGX;
    vector<int> patchVertexX(nScoresX, vacant);
    vector<vector<int> > TpatchX;

    for (i = 0; i < TX.size(); ++i) {
        if (TX[i].size() == numCyclesX) {
            fullX = i;
            break;
        }
    }

    if(fullX != vacant){
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
            if(mateInducedX[cycleVertexX[TX[fullX][v]]][saveX + 1] == TX[fullX][v]){
                for(i = saveX + 1; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
                for(i = mateInducedX[cycleVertexX[TX[fullX][v]]].size(); i-- > saveX +1;){ //from end of cycle to element at position save+1
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
            }

                //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
                //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
            else if (mateInducedX[cycleVertexX[TX[fullX][v]]][saveX - 1] == TX[fullX][v]){
                for(i = saveX; i < mateInducedX[cycleVertexX[TX[fullX][v]]].size(); ++i){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
                for(i = 0; i < saveX; ++i){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
            }

                //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
                //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
            else if (saveX == 0 && mateInducedX[cycleVertexX[TX[fullX][v]]][mateInducedX[cycleVertexX[TX[fullX][v]]].size()-1] == TX[fullX][v]){
                for(i = 0; i < mateInducedX[cycleVertexX[TX[fullX][v]]].size(); ++i){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
            }

                //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
                //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
            else if(saveX == mateInducedX[cycleVertexX[TX[fullX][v]]].size()-1 && mateInducedX[cycleVertexX[TX[fullX][v]]][0] == TX[fullX][v]){
                for(i = mateInducedX[cycleVertexX[TX[fullX][v]]].size(); i-- > 0; ){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TX[fullX][v]]][i]);
                }
            }
            //END ONETCYCLEPATCH FUNCTION
            //endregion
        }

        //region MAKEPATH
        /********************************* MAKEPATH FUNCTION **************************************/
        for(i = 0; i < fullCycleX.size()-1; ++i){
            if((fullCycleX[i] == nScoresX - 1 && fullCycleX[i+1] == nScoresX - 2) || (fullCycleX[i] == nScoresX - 2 && fullCycleX[i+1] == nScoresX - 1)){
                if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                    for(j = 2; j < fullCycleX.size(); ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;
                }

                else if(i == fullCycleX.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                    for(j = 0; j < fullCycleX.size()-2; ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;
                }
                else{ //if the dominating vertices are in the middle of the fullCycle vector
                    for(j = i+2; j < fullCycleX.size(); ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    for(j = 0; j < i; ++j){
                        completePathX.push_back(fullCycleX[j]);
                    }
                    break;

                }

            }
        }

        for(i = 0; i < completePathX.size(); ++i){
            finalX.push_back(originalX[orderX[completePathX[i]]]);
        }
        //END MAKE PATH FUNCTION
        //endregion
        /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/
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
            if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == TpatchX[u][v]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = 0; i < saveX; ++i) {
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
            else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == matchListX[TpatchX[u][v]]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = 0; i < saveX; ++i) {
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
            else if(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == matchListX[TpatchX[u][v]]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = saveX; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); i-- > saveX + 1;) { //from end of cycle to element at position save+1
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
            else if(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == TpatchX[u][v]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = saveX; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); i-- > saveX + 1;) { //from end of cycle to element at position save+1
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                //the last element in the cycle is Tpatch[u][v]
            else if(saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1] == TpatchX[u][v]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                //the first element in the cycle is Tpatch[u][v]
            else if(saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1 && mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == TpatchX[u][v]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1; i-- > 0; ){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                //and the last element in the cycle is matchList[Tpatch[u][v]]
            else if(saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1] == matchListX[TpatchX[u][v]]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                //the first element in the cycle is matchList[Tpatch[u][v]]
            else if(saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1 && mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == matchListX[TpatchX[u][v]]){
                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                for(i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1; i-- > 0; ){
                    if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                        fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                    }
                    else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
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
                        } else {
                            v = ++i;
                        }
                        for (j = 0; j < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++j) {
                            if (mateInducedX[cycleVertexX[TpatchX[u][v]]][j] == matchListX[TpatchX[u][v]]) {
                                saveX = j;
                                break;
                            }
                        }
                        break;
                    } else if (matchListX[TpatchX[u][i]] == fullCycleX.back()) {
                        if (i == 0) {
                            v = TpatchX[u].size() - 1;
                        } else {
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
                if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == TpatchX[u][v]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = 0; i < saveX; ++i) {
                            if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                else if (mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX - 1] == matchListX[TpatchX[u][v]]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = saveX + 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = 0; i < saveX; ++i) {
                            if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                else if(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == matchListX[TpatchX[u][v]]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = saveX; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); i-- > saveX + 1;) { //from end of cycle to element at position save+1
                            if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                else if(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX + 1] == TpatchX[u][v]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = saveX; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); i-- > saveX + 1;) { //from end of cycle to element at position save+1
                            if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            }
                            else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                                fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                                u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                    //the last element in the cycle is Tpatch[u][v]
                else if(saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1] == TpatchX[u][v]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                    //the first element in the cycle is Tpatch[u][v]
                else if(saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1 && mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == TpatchX[u][v]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1; i-- > 0; ){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                    //and the last element in the cycle is matchList[Tpatch[u][v]]
                else if(saveX == 0 && mateInducedX[cycleVertexX[TpatchX[u][v]]][mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1] == matchListX[TpatchX[u][v]]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = 1; i < mateInducedX[cycleVertexX[TpatchX[u][v]]].size(); ++i){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                            u = patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                    //the first element in the cycle is matchList[Tpatch[u][v]]
                else if(saveX == mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1 && mateInducedX[cycleVertexX[TpatchX[u][v]]][0] == matchListX[TpatchX[u][v]]){
                    fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][saveX]);
                    for(i = mateInducedX[cycleVertexX[TpatchX[u][v]]].size()-1; i-- > 0; ){
                        if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] == vacant) {
                            fullCycleX.push_back(mateInducedX[cycleVertexX[TpatchX[u][v]]][i]);
                        }
                        else if(patchVertexX[mateInducedX[cycleVertexX[TpatchX[u][v]]][i]] != vacant){
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
            /******************** MAKEPATH FUNCTION *************************************/
            for(i = 0; i < fullCycleX.size()-1; ++i){
                if((fullCycleX[i] == nScoresX - 1 && fullCycleX[i+1] == nScoresX - 2) || (fullCycleX[i] == nScoresX - 2 && fullCycleX[i+1] == nScoresX - 1)){
                    if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                        for(j = 2; j < fullCycleX.size(); ++j){
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }

                    else if(i == fullCycleX.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                        for(j = 0; j < fullCycleX.size()-2; ++j){
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;
                    }
                    else{ //if the dominating vertices are in the middle of the fullCycle vector
                        for(j = i+2; j < fullCycleX.size(); ++j){
                            completePathX.push_back(fullCycleX[j]);
                        }
                        for(j = 0; j < i; ++j){
                            completePathX.push_back(fullCycleX[j]);
                        }
                        break;

                    }

                }
            }

            for(i = 0; i < completePathX.size(); ++i){
                finalX.push_back(originalX[orderX[completePathX[i]]]);
            }
            //END MAKE PATH FUNCTION
            //endregion

            /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/

        }
        else if (SSumX < numCyclesX) {
            //cout << instance << ": Infeasible SSum < numCycles\n\n";
            /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/

        }
        else {
            //cout << instance << ": Problem.\n\n";
            /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
        }


    }
    //endregion

    //END MBAHRA STRIPX

    /*********************************************************************************************************/

    /**FOR NEW STRIPY PAIRPAIR**/

    int nScoresY, nBoxY, nCompY;
    vector<int> scoresY;
    vector<int> orderY;
    vector<int> originalY;
    vector<int> finalY;

    /**FOR NEW STRIPY PAIRPAIR**/
    for(k = 0; k < stripY[j1].size(); ++k){
        if(k == c1 || k == c1 + 1 || k == d1 || k == d1 + 1){
            continue;
        }
        scoresY.push_back(allScores[stripY[j1][k]]);
        originalY.push_back(stripY[j1][k]);

    }
    scoresY.push_back(allScores[stripX[i1][a1]]);
    originalY.push_back(stripX[i1][a1]);
    scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
    originalY.push_back(stripX[i1][a1+1]);
    scoresY.push_back(allScores[stripX[i1][b1]]);
    originalY.push_back(stripX[i1][b1]);
    scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
    originalY.push_back(stripX[i1][b1+1]);
    scoresY.push_back(70);
    scoresY.push_back(70);

    //region Initialisation
    nScoresY = scoresY.size();
    nBoxY = scoresY.size() /2;
    nCompY = (nBoxY + (nBoxY % 2)) / 2;
    vector<int> invOrderY(nScoresY);
    vector<vector<int> > adjMatY(nScoresY, vector<int>(nScoresY, 0));
    vector<int> matesY(nScoresY, 0);
    vector<int> completePathY;

    for(k = 0; k < nScoresY; ++k){
        orderY.push_back(k);
    }

    cout << "Scores:\n";
    for(k = 0; k < nScoresY; ++k){
        cout << scoresY[k] << " ";
    }
    cout << endl << endl;

    cout << "Order:\n";
    for(k = 0; k < nScoresY; ++k){
        cout << orderY[k] << " ";
    }
    cout << endl << endl;

    for(i = 1; i < nScoresY; ++i){
        for(j = i-1; j >= 0; --j){
            if(scoresY[i] < scoresY[orderY[j]]){
                orderY[j+1] = orderY[j];
                orderY[j] = i;
            }
        }
    }

    cout << "Order:\n";
    for(k = 0; k < nScoresY; ++k){
        cout << orderY[k] << " ";
    }
    cout << endl << endl;


    for(k = 0; k < nScoresY; ++k){
        invOrderY[orderY[k]] = k;
    }

    cout << "Inverse Order:\n";
    for(k = 0; k < nScoresY; ++k){
        cout << invOrderY[k] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < nScoresY - 1; i+=2){
        adjMatY[invOrderY[i]][invOrderY[i+1]] = 2;
        adjMatY[invOrderY[i+1]][invOrderY[i]] = 2;
    }

    sort(scoresY.begin(), scoresY.end());

    for (i = 0; i < scoresY.size() - 1; ++i) {
        for (j = i + 1; j < scoresY.size(); ++j) {
            if (scoresY[i] + scoresY[j] >= threshold && adjMatY[i][j] != 2) {
                adjMatY[i][j] = 1;
                adjMatY[j][i] = 1;
            }
        }

    }

    cout << "Adjacency Matrix\n";
    for(i = 0; i < nScoresY; ++i){
        for(j = 0; j < nScoresY; ++j){
            cout << adjMatY[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;


    for (i = 0; i < nScoresY; ++i) {
        for (j = 0; j < nScoresY; ++j) {
            if (adjMatY[i][j] == 2) {
                matesY[i] = j;
                break;
            }
        }
    }
    cout << "Mates Vector:\n";
    for(i = 0; i < matesY.size(); ++i){
        cout << matesY[i] << " ";
    }
    cout << endl << endl;
    //endregion


    //region MTGMA
    /**MTGMA**/
    int lastMatchY = vacant;
    int mateMatchY = vacant;
    int vacantFlagY = 0;
    int matchSizeY = 0;
    vector<int> cycleVertexY(nScoresY, 1);
    vector<int> matchListY(nScoresY, vacant);

    for (i = 0; i < nScoresY; ++i) { //check all vertices
        vacantFlagY = 0;
        if (matchListY[i] == vacant) { //if vertex has not yet been matched
            for (j = nScoresY - 1; j > i; --j) { //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
                if (adjMatY[i][j] == 1 && matchListY[j] == vacant) { //if vertices i and j are adjacent, and if vertex j has not yet been matched
                    matchListY[i] = j;
                    matchListY[j] = i;
                    lastMatchY = i;
                    ++matchSizeY;
                    if (vacantFlagY == 1) { //delete edge for FCA if matching was not with highest vertex due to the highest vertex being its mate
                        cycleVertexY[i] = vacant;
                        cycleVertexY[j] = vacant;
                    }
                    break;
                }
                else if (adjMatY[i][j] == 2 && matchListY[j] == vacant) { //if potential match == mate
                    vacantFlagY = 1;
                }
            }//end for j
            if (matchListY[i] == vacant) { //if vertex has still not been matched
                for (k = 0; k < nScoresY - 2; ++k) {
                    if (adjMatY[i][k] == 2) { //if vertex i and vertex k are mates
                        mateMatchY = k;
                        break;
                    }
                }
                if ((scoresY[i] + scoresY[mateMatchY] >= threshold) //match with mate?
                    && (matchListY[mateMatchY] == vacant) //is mate unmatched?
                    && (lastMatchY != vacant) //has the previous vertex been matched?
                    && (mateMatchY > i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
                    && (scoresY[lastMatchY] + scoresY[mateMatchY] >= threshold)) { //can mate be matched with last matched vertex?
                    // if so, then swap mates
                    matchListY[i] = matchListY[lastMatchY];
                    matchListY[lastMatchY] = mateMatchY;
                    matchListY[mateMatchY] = lastMatchY;
                    matchListY[matchListY[i]] = i;
                    cycleVertexY[lastMatchY] = vacant; //edge from mate swap will not count for FCA
                    cycleVertexY[mateMatchY] = vacant; //edge from mate swap will not count for FCA
                    lastMatchY = i;
                    ++matchSizeY;
                }
            }//end if matchList == vacant
        }//end if matchList[i] == i
    }//end for i


    /*cout << "Cycle Vertex vector after MTGMA:\n";
    for(i = 0; i < cycleVertexY.size(); ++i){
        cout << cycleVertexY[i] << " ";
    }
    cout << endl << endl;*/

    /*cout << "Matching List:\n";
    for(i = 0; i < matchListY.size(); ++i){
        cout << matchListY[i] << " ";
    }
    cout << endl << endl;*/

    //endregion
    if(matchSizeY < nBoxY){
        //NOT ENOUGH MATCHING EDGES
        /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
    }

    //region MIS
    /**MIS**/
    int numCyclesY = 0;
    int smallestVertexY;
    int currentVertexY;
    vector<vector<int> > mateInducedY;
    vector<int> lengthMateInducedY;
    vector<int> tempMISY;
    vector<int> checkedY(nScoresY, 0);
    vector<int> fullCycleY;

    //find the smallest vertex not yet checked for mate-induced structure - start with this vertex
    for (i = 0; i < nScoresY; ++i) {
        if (checkedY[i] == 0) {
            smallestVertexY = i;
            break;
        }
    }

    //Building the mate-induced structure
    do {
        currentVertexY = smallestVertexY;
        do {
            tempMISY.push_back(currentVertexY);
            checkedY[currentVertexY] = 1;
            tempMISY.push_back(matesY[currentVertexY]);
            checkedY[matesY[currentVertexY]] = 1;
            currentVertexY = matchListY[matesY[currentVertexY]];
        } while (currentVertexY != smallestVertexY);

        mateInducedY.push_back(tempMISY);
        tempMISY.clear();

        for (i = 0; i < nScoresY; ++i) {
            if (checkedY[i] == 0) {
                smallestVertexY = i;
                break;
            }
        }


    } while (smallestVertexY != currentVertexY);

    tempMISY.clear(); //clear cycle vector again for next instance

    numCyclesY = mateInducedY.size(); //number of cycles in the mate-induced structure

   /* cout << "Mate-Induced Structure:\n";
    for(i = 0; i < mateInducedY.size(); ++i){
        for(j = 0; j < mateInducedY[i].size(); ++j){
            cout << mateInducedY[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;*/
    //cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

    for (i = 0; i < mateInducedY.size(); ++i) {
        lengthMateInducedY.push_back(mateInducedY[i].size());
    }

    //endregion
    if(lengthMateInducedY[0] == nScoresY){
        for(j = 0; j < mateInducedY[0].size(); ++j){
            fullCycleY.push_back(mateInducedY[0][j]);
        }
        //region MAKEPATH
        /******************** MAKEPATH FUNCTION *************************************/
        for(i = 0; i < fullCycleY.size()-1; ++i){
            if((fullCycleY[i] == nScoresY - 1 && fullCycleY[i+1] == nScoresY - 2) || (fullCycleY[i] == nScoresY - 2 && fullCycleY[i+1] == nScoresY - 1)){
                if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                    for(j = 2; j < fullCycleY.size(); ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;
                }

                else if(i == fullCycleY.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                    for(j = 0; j < fullCycleY.size()-2; ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;
                }
                else{ //if the dominating vertices are in the middle of the fullCycle vector
                    for(j = i+2; j < fullCycleY.size(); ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    for(j = 0; j < i; ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;

                }

            }
        }
        for(i = 0; i < completePathY.size(); ++i){
            finalY.push_back(originalY[orderY[completePathY[i]]]);
        }


        //END MAKE PATH FUNCTION
        //endregion
        /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/
    }

    //region FCA
    /**FCA**/
    int qstarY;
    int numEdgesY;
    vector<vector<int> > SY(nCompY, vector<int>(nCompY, 0));
    vector<vector<int> > TY;
    vector<int> edgeY;
    vector<int> tY;

    //create list cycleVertex that contains for each vertex the cycle that each edge belongs to
    for (i = 0; i < mateInducedY.size(); ++i) {
        for (j = 0; j < mateInducedY[i].size(); ++j) {
            if (cycleVertexY[mateInducedY[i][j]] != vacant) { //if edge is not deleted for FCA
                cycleVertexY[mateInducedY[i][j]] = i;
            }
        }
    }

    /*cout << "Cycle Vertex:\n";
    for (i = 0; i < cycleVertexY.size(); ++i) {
        cout << cycleVertexY[i] << " ";
    }
    cout << endl << endl;*/

    //create list of edges without empty edges (those generated by mate swap)
    for (i = 0; i < matchSizeY; ++i) {
        while (cycleVertexY[i] == vacant) {
            ++i;
        }
        edgeY.push_back(i);
    }
    numEdgesY = edgeY.size();

    /*cout << "Edges vector:\n";
    for (i = 0; i < edgeY.size(); ++i) {
        cout << edgeY[i] << " ";
    }
    cout << endl << endl;*/

    //cout << "Number of Edges: " << numEdges << endl;
    //FCA Algorithm
    qstarY = -1;
    k = 0; //edge from matching that is under consideration

    do {
        while (k < numEdgesY - 2 && (adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] != 1 || cycleVertexY[edgeY[k]] == cycleVertexY[edgeY[k + 1]])) {
            ++k;
        }
        if (adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] == 1 && cycleVertexY[edgeY[k]] != cycleVertexY[edgeY[k + 1]]) {
            ++qstarY;
            tY.push_back(edgeY[k]);
            SY[qstarY][cycleVertexY[edgeY[k]]] = 1;
            while (k < numEdgesY - 1 && adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] == 1 && SY[qstarY][cycleVertexY[edgeY[k + 1]]] == 0) { //add more edges to current T-cycle
                ++k;
                tY.push_back(edgeY[k]);
                SY[qstarY][cycleVertexY[edgeY[k]]] = 1;
            }
            TY.push_back(tY);
            tY.clear();
        } // end if
        ++k;
    } while (k < numEdgesY - 1);

    tY.clear();

    /*cout << "T matrix:\n";
    for(i = 0; i < TY.size(); ++i){
        for(j = 0; j < TY[i].size(); ++j){
            cout << TY[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl << endl;*/

    /*cout << "S Matrix:\n";
    for(i = 0; i < SY.size(); ++i){
        for(j = 0; j < SY[i].size(); ++j){
            cout << SY[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;*/

    //cout << "qstar: " << qstarY << endl << endl;

    //endregion
    if(qstarY == -1){
        //INFEASIBLE, EXIT
        /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
    }

    //region PATCHGRAPH
    /**PATCHGRAPH**/
    int saveY, SSumY, SqIntSY;
    int fullY = vacant;
    vector<int> QSetY(qstarY, 0);
    vector<int> SSetY;
    vector<int> patchCycleY(qstarY, vacant);
    vector<int> tempPGY;
    vector<int> patchVertexY(nScoresY, vacant);
    vector<vector<int> > TpatchY;

    for (i = 0; i < TY.size(); ++i) {
        if (TY[i].size() == numCyclesY) {
            fullY = i;
            break;
        }
    }

    if(fullY != vacant){
        saveY = 0;
        //cout << "Full: " << full << endl;
        for (v = 0; v < TY[fullY].size(); ++v) {
            for (j = 0; j < mateInducedY[cycleVertexY[TY[fullY][v]]].size(); ++j) {
                if (mateInducedY[cycleVertexY[TY[fullY][v]]][j] == matchListY[TY[fullY][v]]) {
                    saveY = j;
                    break;
                }
            }
            //region oneTCyclePatch
            /****************************oneTCyclePatch algorithm: ***********************************/
            //CASE ONE: if element matchList[T[full][v]] is before element T[full][v] in the mateInduced cycle
            //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save + 1" is T[full][v]
            if(mateInducedY[cycleVertexY[TY[fullY][v]]][saveY + 1] == TY[fullY][v]){
                for(i = saveY + 1; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
                for(i = mateInducedY[cycleVertexY[TY[fullY][v]]].size(); i-- > saveY +1;){ //from end of cycle to element at position save+1
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
            }

                //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
                //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
            else if (mateInducedY[cycleVertexY[TY[fullY][v]]][saveY - 1] == TY[fullY][v]){
                for(i = saveY; i < mateInducedY[cycleVertexY[TY[fullY][v]]].size(); ++i){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
                for(i = 0; i < saveY; ++i){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
            }

                //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
                //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
            else if (saveY == 0 && mateInducedY[cycleVertexY[TY[fullY][v]]][mateInducedY[cycleVertexY[TY[fullY][v]]].size()-1] == TY[fullY][v]){
                for(i = 0; i < mateInducedY[cycleVertexY[TY[fullY][v]]].size(); ++i){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
            }

                //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
                //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
            else if(saveY == mateInducedY[cycleVertexY[TY[fullY][v]]].size()-1 && mateInducedY[cycleVertexY[TY[fullY][v]]][0] == TY[fullY][v]){
                for(i = mateInducedY[cycleVertexY[TY[fullY][v]]].size(); i-- > 0; ){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                }
            }
            //END ONETCYCLEPATCH FUNCTION
            //endregion
        }

        //region MAKEPATH
        /******************** MAKEPATH FUNCTION *************************************/
        for(i = 0; i < fullCycleY.size()-1; ++i){
            if((fullCycleY[i] == nScoresY - 1 && fullCycleY[i+1] == nScoresY - 2) || (fullCycleY[i] == nScoresY - 2 && fullCycleY[i+1] == nScoresY - 1)){
                if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                    for(j = 2; j < fullCycleY.size(); ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;
                }

                else if(i == fullCycleY.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                    for(j = 0; j < fullCycleY.size()-2; ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;
                }
                else{ //if the dominating vertices are in the middle of the fullCycle vector
                    for(j = i+2; j < fullCycleY.size(); ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    for(j = 0; j < i; ++j){
                        completePathY.push_back(fullCycleY[j]);
                    }
                    break;

                }

            }
        }
        for(i = 0; i < completePathY.size(); ++i){
            finalY.push_back(originalY[orderY[completePathY[i]]]);
        }
        //END MAKE PATH FUNCTION
        //endregion
        /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/
    }

    else {
        q = 0; //Start with first Tq-cycle
        QSetY[0] = 1;
        SSumY = 0; //number of MIS-cycles that have been included
        for (i = 0; i < numCyclesY; ++i) {
            SSetY.push_back(SY[q][i]); // ==1 if MIS cycle i has been included
        }
        for (i = 0; i < numCyclesY; ++i) {
            SSumY = SSumY + SSetY[i];
        }

        if (SSumY >= 1) {
            patchCycleY[q] = 1;
        }

        //Start connectivity check
        while (q <= qstarY && SSumY < numCyclesY) {
            do {
                ++q;
                SqIntSY = vacant;
                if (q <= qstarY) {
                    for (j = 0; j < numCyclesY; ++j) { //is there a j such that S[q][j] = 1 and SSet[j] = 1?
                        if (SY[q][j] == 1 && SSetY[j] == 1) {
                            SqIntSY = 1;
                            //break here? no need to check all other j indices once one has been found such that S[q][j] =1 and SSet[j] = 1
                        }
                    }
                }
            } while (q < qstarY + 1 && (QSetY[q] == 1 || SqIntSY == vacant));

            if (q <= qstarY) { //if Tq-cyce for enlargement has been found
                for (i = 0; i < numCyclesY; ++i) {
                    if (SSetY[i] == 0 && SY[q][i] == 1) {
                        SSetY[i] = 1;
                        ++SSumY;
                        patchCycleY[q] = 1;
                    }
                }
                QSetY[q] = 1;
                q = 0;
            }
        }//end while


        //If patching graph is connected, then instance is feasible, else infeasible
        if (SSumY == numCyclesY) {
            for (i = 0; i < patchCycleY.size(); ++i) {
                if (patchCycleY[i] == 1) {
                    for (j = 0; j < TY[i].size(); ++j) {
                        tempPGY.push_back(TY[i][j]);
                    }
                    TpatchY.push_back(tempPGY);
                    tempPGY.clear();
                }
            }
            tempPGY.clear();

            for (i = 0; i < TpatchY.size(); ++i) {
                for (j = 0; j < TpatchY[i].size(); ++j) {
                    patchVertexY[TpatchY[i][j]] = i;
                    patchVertexY[matchListY[TpatchY[i][j]]] = i;
                }
            }

            u = 0;
            v = 0;
            saveY = 0;
            for (j = 0; j < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++j) {
                if (mateInducedY[cycleVertexY[TpatchY[u][v]]][j] == matchListY[TpatchY[u][v]]) {
                    saveY = j;
                    break;
                }
            }
            //region multipleTCyclePatch
            /***************************** MULTIPLE T CYCLE PATCH FUNCTION ***********************/
            int x = 0;

            //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
            if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == TpatchY[u][v]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = 0; i < saveY; ++i) {
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
            else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == matchListY[TpatchY[u][v]]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = 0; i < saveX; ++i) {
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
            else if(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == matchListY[TpatchY[u][v]]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = saveY; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); i-- > saveY + 1;) { //from end of cycle to element at position save+1
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
            else if(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == TpatchY[u][v]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = saveY; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        x = 1;
                        break;
                    }
                }
                if(x == 0) {
                    for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); i-- > saveY + 1;) { //from end of cycle to element at position save+1
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }
            }

                //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                //the last element in the cycle is Tpatch[u][v]
            else if(saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1] == TpatchY[u][v]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                //the first element in the cycle is Tpatch[u][v]
            else if(saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1 && mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == TpatchY[u][v]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1; i-- > 0; ){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                //and the last element in the cycle is matchList[Tpatch[u][v]]
            else if(saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1] == matchListY[TpatchY[u][v]]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        break;
                    }
                }
            }

                //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                //the first element in the cycle is matchList[Tpatch[u][v]]
            else if(saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1 && mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == matchListY[TpatchY[u][v]]){
                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                for(i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1; i-- > 0; ){
                    if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                    }
                    else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                        break;
                    }
                }
            }
            //END MULTIPLE T CYCLE PATCH FUNCTION
            //endregion

            while (fullCycleY.size() < nScoresY) {
                saveY = 0;
                for (i = 0; i < TpatchY[u].size(); ++i) {
                    if (TpatchY[u][i] == fullCycleY.back()) {
                        if (i == TpatchY[u].size() - 1) {
                            v = 0;
                        } else {
                            v = ++i;
                        }
                        for (j = 0; j < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++j) {
                            if (mateInducedY[cycleVertexY[TpatchY[u][v]]][j] == matchListY[TpatchY[u][v]]) {
                                saveY = j;
                                break;
                            }
                        }
                        break;
                    }
                    else if (matchListY[TpatchY[u][i]] == fullCycleY.back()) {
                        if (i == 0) {
                            v = TpatchY[u].size() - 1;
                        } else {
                            v = --i;
                        }
                        for (j = 0; j < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++j) {
                            if (mateInducedY[cycleVertexY[TpatchY[u][v]]][j] == TpatchY[u][v]) {
                                saveY = j;
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
                if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == TpatchY[u][v]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = 0; i < saveY; ++i) {
                            if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == matchListY[TpatchY[u][v]]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = 0; i < saveX; ++i) {
                            if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                else if(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == matchListY[TpatchY[u][v]]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = saveY; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); i-- > saveY + 1;) { //from end of cycle to element at position save+1
                            if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                else if(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == TpatchY[u][v]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = saveY; i-- > 0;){ //from element at position 'save' to the first element in the cycle
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if(x == 0) {
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); i-- > saveY + 1;) { //from end of cycle to element at position save+1
                            if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                    //the last element in the cycle is Tpatch[u][v]
                else if(saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1] == TpatchY[u][v]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                    //the first element in the cycle is Tpatch[u][v]
                else if(saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1 && mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == TpatchY[u][v]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1; i-- > 0; ){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                    //and the last element in the cycle is matchList[Tpatch[u][v]]
                else if(saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1] == matchListY[TpatchY[u][v]]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                    //the first element in the cycle is matchList[Tpatch[u][v]]
                else if(saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1 && mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == matchListY[TpatchY[u][v]]){
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for(i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size()-1; i-- > 0; ){
                        if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if(patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant){
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }
                //END MULTIPLE T CYCLE PATCH FUNCTION
                //endregion

            }


            //region MAKEPATH
            /******************** MAKEPATH FUNCTION *************************************/
            for(i = 0; i < fullCycleY.size()-1; ++i){
                if((fullCycleY[i] == nScoresY - 1 && fullCycleY[i+1] == nScoresY - 2) || (fullCycleY[i] == nScoresY - 2 && fullCycleY[i+1] == nScoresY - 1)){
                    if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                        for(j = 2; j < fullCycleY.size(); ++j){
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }

                    else if(i == fullCycleY.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                        for(j = 0; j < fullCycleY.size()-2; ++j){
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }
                    else{ //if the dominating vertices are in the middle of the fullCycle vector
                        for(j = i+2; j < fullCycleY.size(); ++j){
                            completePathY.push_back(fullCycleY[j]);
                        }
                        for(j = 0; j < i; ++j){
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;

                    }

                }
            }
            for(i = 0; i < completePathY.size(); ++i){
                finalY.push_back(originalY[orderY[completePathY[i]]]);
            }
            //END MAKE PATH FUNCTION
            //endregion

            /*********************************** OUTPUT STRIPX, MOVE ONTO STRIP Y ********************************/

        }
        else if (SSumY < numCyclesY) {
            //cout << instance << ": Infeasible SSum < numCycles\n\n";
            /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/

        }
        else {
            //cout << instance << ": Problem.\n\n";
            /*************************************** EXIT- NO FEASIBLE SOLUTION **********************************/
        }


    }
    //endregion

    //END MBAHRA STRIPY





}//end void MBAHRA



























































































