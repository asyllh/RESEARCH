/*--------------/
ALH
packing.cpp
Evolutionary Algorithm with Local Search
05/12/2017
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

void FFD(int numScores, int numBox, int maxBoxWidth, int maxStripWidth, vector<int> &mates,
         vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, mini, k, l;
    int min = 0;
    int max = maxBoxWidth;
    vector<int> boxOrder;
    vector<int> checked(numScores, 0);


    while(boxOrder.size() < numBox) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] <= max) {
                min = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxOrder.push_back(mini);
        checked[mini] = 1;
        checked[mates[mini]] = 1;
        max = min;
        min = 0;
    }


    cout << "Box decrease:\n";
    for(i = 0; i < boxOrder.size(); ++i){
        cout << boxOrder[i] << " ";
    }
    cout << endl << endl;

    strip[0].push_back(boxOrder[0]);
    strip[0].push_back(mates[boxOrder[0]]);
    stripSum[0] += boxWidths[boxOrder[0]][mates[boxOrder[0]]];


    for(j = 1; j < boxOrder.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[boxOrder[j]][mates[boxOrder[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][boxOrder[j]] == 1){
                        strip[i].push_back(boxOrder[j]);
                        strip[i].push_back(mates[boxOrder[j]]);
                        stripSum[i] += boxWidths[boxOrder[j]][mates[boxOrder[j]]];
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][mates[boxOrder[j]]] == 1){
                        strip[i].push_back(mates[boxOrder[j]]);
                        strip[i].push_back(boxOrder[j]);
                        stripSum[i] += boxWidths[boxOrder[j]][mates[boxOrder[j]]];
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(boxOrder[j]);
                strip[i].push_back(mates[boxOrder[j]]);
                stripSum[i] += boxWidths[boxOrder[j]][mates[boxOrder[j]]];
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



    cout << "After FFD: " << strip.size() << " strips\n";

    //cout << "Lower Bound: " << lowerBound(totalBoxWidth, maxStripWidth) << " strips\n";


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
            cout << i << setw(9) << stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;



}

void partialFFD(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &mates, vector<vector<int> > &adjMatrix,
                vector<vector<int> > &boxWidths, vector<int> &partialBoxes, vector<int> &partialSum, vector<vector<int> > &partialSol){

    int i, j, mini, k, l;
    int min = 0;
    int max = maxBoxWidth;
    vector<int> boxDecrease;
    vector<int> checked(numScores, 0);


    while(boxDecrease.size() < partialBoxes.size()/2) {
        for (i = 0; i < partialBoxes.size(); ++i) {
            if(checked[partialBoxes[i]] == 1){
                continue;
            }
            if (boxWidths[partialBoxes[i]][mates[partialBoxes[i]]] > min && boxWidths[partialBoxes[i]][mates[partialBoxes[i]]] <= max) {
                min = boxWidths[partialBoxes[i]][mates[partialBoxes[i]]];
                mini = partialBoxes[i];
            }
        }
        boxDecrease.push_back(mini);
        checked[mini] = 1;
        checked[mates[mini]] = 1;
        max = min;
        min = 0;
    }

    cout << "\nBox decrease:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;

    partialSol[0].push_back(boxDecrease[0]);
    partialSol[0].push_back(mates[boxDecrease[0]]);
    partialSum[0] += boxWidths[boxDecrease[0]][mates[boxDecrease[0]]];


    for(j = 1; j < boxDecrease.size(); ++j){
        for(i = 0; i < partialSol.size(); ++i){
            if(!partialSol[i].empty()){
                if(partialSum[i] + boxWidths[boxDecrease[j]][mates[boxDecrease[j]]] <= maxStripWidth){
                    if(adjMatrix[partialSol[i].back()][boxDecrease[j]] == 1){
                        partialSol[i].push_back(boxDecrease[j]);
                        partialSol[i].push_back(mates[boxDecrease[j]]);
                        partialSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        break;
                    }
                    else if (adjMatrix[partialSol[i].back()][mates[boxDecrease[j]]] == 1){
                        partialSol[i].push_back(mates[boxDecrease[j]]);
                        partialSol[i].push_back(boxDecrease[j]);
                        partialSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                        break;
                    }
                }
            }
            else if (partialSol[i].empty()){
                partialSol[i].push_back(boxDecrease[j]);
                partialSol[i].push_back(mates[boxDecrease[j]]);
                partialSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                break;
            }
        }
    }

    k = partialSol.size() - 1;
    while(partialSol[k].empty()){
        partialSol.pop_back();
        --k;
    }

    l = partialSum.size() -1;
    while(partialSum[l] == 0){
        partialSum.pop_back();
        --l;
    }




    cout << "After REPAIR FFD: " << partialSol.size() << " strips\n";


    cout << "REPAIR Strips FFD (scores):\n";
    for(i = 0; i < partialSol.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < partialSol[i].size(); ++j){
            cout << partialSol[i][j] << " ";
        }
        cout << endl;
    }


    cout << "stripSum: ";
    for(i = 0; i < partialSum.size(); ++i){
        cout << partialSum[i] << " ";
    }
    cout << endl << endl;

}

void createInitialPopulation(int numScores, int numBox, int maxBoxWidth, int maxStripWidth, vector<int> &allScores,
                             vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths,
                             vector<vector<int> > &populationSum, vector<vector<vector<int> > > &population){

    int i, j;
    vector<vector<int> > strip(numBox);
    vector<int> stripSum(numBox, 0);

    FFD(numScores, numBox, maxBoxWidth, maxStripWidth, mates, adjMatrix, boxWidths, stripSum, strip);

    mutation(numScores, maxBoxWidth, maxStripWidth, allScores, mates, adjMatrix, boxWidths, stripSum, strip);

    population.push_back(strip);
    populationSum.push_back(stripSum);

}

void mutation(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &allScores, vector<int> &mates,
              vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j;
    vector<int> stripSumX;
    vector<vector<int> > stripX;
    vector<int> stripSumY;
    vector<vector<int> > stripY;

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

    localSearch(numScores, maxBoxWidth, maxStripWidth, allScores, mates, adjMatrix, boxWidths, stripSum, strip, stripSumX, stripX, stripSumY, stripY);





}

void localSearch(int numScores, int maxBoxWidth, int maxStripWidth, vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix,
                 vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY){

    int a, b, c, d, i, j, k, l, pairSizeX, pairSizeY;
    int swapType, moveType, feasible;

    //region PairPair
    PairPair:
    swapType = 0;
    moveType = 0;
    /*SWAPPING A PAIR OF BOXES FROM EACH SET*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
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
                                        swapType = 1;
                                        cout << "i: " << i << " j: " << j << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
                                        if(stripX[i].size() == 4){ //If stripX[i] only contains 2 boxes
                                            if(stripY[j].size() == 4){ //If stripY[j] only contains 2 boxes
                                                //Do a straight swap, no need for MBAHRA
                                                stripX[i].swap(stripY[j]);
                                                swap(stripSumX[i], stripSumY[j]);
                                                feasible = 1;
                                            }
                                            else if (d == c+2){ //If stripY[j] contains more than 2 boxes & the two chosen boxes in stripY[j] are adjacent
                                                //Only perform MBAHRA on stripY[j]
                                                moveType = 11;
                                                MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else{ //IIf stripY[j] contains more than 2 boxes & boxes c and d are not adjacent
                                                moveType = 0;
                                                MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else if (stripY[j].size() == 4){ //If stripY[j] only contains 2 boxes but stripX[i] contains > 2 boxes
                                            if(b == a+2){ //If the two chosen boxes in stripX[i] are adjacent to one another
                                                //Only perform MBAHRA on stripX[i]
                                                moveType = 12;
                                                MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else{ //If boxes a and b are not adjacent
                                                moveType = 0;
                                                MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else{ //If stripX[i].size() > 4 && stripY[j[.size() > 4
                                            moveType = 0;
                                            MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                        }
                                        if(feasible == 1){
                                            //++count;
                                            //cout << count << ": PairPair\n";
                                            goto PairSin;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //cout << "NO PAIR PAIR SWAP PERFORMED.\n";
    //endregion

    //region PairSin
    PairSin:
    swapType = 0;
    moveType = 0;
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
                                swapType = 2;
                                if(stripX[i].size() == 4){ //If stripX[i] only contains 2 boxes
                                    if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box
                                        //straight swap
                                        stripX[i].swap(stripY[j]);
                                        swap(stripSumX[i], stripSumY[j]);
                                        feasible = 1;
                                    }
                                    else{ //If stripY[j] contains more than 1 box
                                        //Only perform MBAHRA on stripY[j]
                                        moveType = 21;
                                        MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box, but stripX[i] contains > 2 boxes
                                    if(b == a+2){ //If the two boxes chosen from stripX[i] are adjacent to one another
                                        //Only perform MBAHRA on stripX[i]
                                        moveType = 22;
                                        MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                    }
                                    else{ //If the two boxes chosen from stripX[i] are not adjacent to one another
                                        moveType = 0;
                                        MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else { //If stripX[i].size() > 4 && stripY[j].size() > 2
                                    moveType = 0;
                                    MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                                }

                                if(feasible == 1){
                                    //++count;
                                    //cout << count << ": PairSin\n";
                                    goto SinSin;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //cout << "NO PAIR SIN SWAP PERFORMED.\n";
    //endregion

    //region SinSin
    SinSin:
    swapType = 0;
    moveType = 0;
    /*SWAPPING ONE BOX FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        for(a = 0; a < stripX[i].size()-1; a+=2){ // Starting from the first score on the first box until the first score on the last box
            for(j = 0; j < stripY.size(); ++j){ // For each strip in the set stripY
                for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
                    //Check if boxwidth[a] < boxWidth[c] and that box can fit on strip
                    if(boxWidths[stripX[i][a]][stripX[i][a+1]] < boxWidths[stripY[j][c]][stripY[j][c+1]]
                       && stripSumX[i] - boxWidths[stripX[i][a]][stripX[i][a+1]] + boxWidths[stripY[j][c]][stripY[j][c+1]] <= maxStripWidth){
                        swapType = 3;
                        if(stripX[i].size() == 2){ //If stripX[i] only contains 1 box
                            if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box
                                //straight swap
                                stripX[i].swap(stripY[j]);
                                swap(stripSumX[i], stripSumY[j]);
                                feasible = 1;
                            }
                            else{ //If stripY[j] contains more than 1 box
                                //Only peform MBAHRA on stripY[j]
                                moveType = 31;
                                MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                            }
                        }
                        else if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box but stripX[i].size() > 2
                            //Only perform MBAHRA on stripX[i]
                            moveType = 32;
                            MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                        }
                        else{ //stripX[i].size() > 2 && stripY[j].size() > 2
                            moveType = 0;
                            MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                        }
                        if(feasible == 1){
                            //++count;
                            //cout << count << ": SinSin\n";
                            goto MoveSin;
                        }
                    }
                }
            }
        }
    }
    //cout << "NO SINGLE SINGLE SWAP PERFORMED.\n";
    //endregion

    //region MoveSin
    MoveSin:
    swapType = 0;
    moveType = 0;
    /*MOVING ONE BOX FROM SET STRIPY TO SET STRIPX*/
    for(j = 0; j < stripY.size(); ++j){ //For each strip in the set stripY
        for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
            for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
                if(stripSumX[i] + boxWidths[stripY[j][c]][stripY[j][c+1]] <= maxStripWidth){
                    swapType = 4;
                    if(stripY[j].size() == 2){ //If stripY[j] only contains one box
                        moveType = 41;
                        MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                    }
                    else {
                        MBAHRA(swapType, moveType, feasible, i, a, b, j, c, d, allScores, boxWidths, stripSumX, stripX, stripSumY, stripY);
                    }
                    if(feasible == 1){
                        //++count;
                        //cout << count << ": MoveSin\n";
                        goto PairPair;
                    }
                }
            }
        }
    }
    cout << "NO SINGLE MOVE PERFORMED.\n\n";
    //endregion

    //End:
    cout << "Local Search complete\n-------------------\n";
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

    //Do FFD on stripY

    vector<int> partialBoxes;
    for(i = 0; i < stripY.size(); ++i){
        for(j = 0; j < stripY[i].size(); ++j){
            partialBoxes.push_back(stripY[i][j]);
        }
    }
    sort(partialBoxes.begin(), partialBoxes.end());

    stripY.clear();
    stripY.resize(partialBoxes.size()/2);
    stripSumY.clear();
    for(i = 0; i < partialBoxes.size()/2; ++i){
        stripSumY.push_back(0);
    }

    partialFFD(numScores, maxBoxWidth, maxStripWidth, mates, adjMatrix, boxWidths, partialBoxes, stripSumY, stripY);

    //join sets stripX and stripY together back into vector<vector<int> > strip

    while(!strip.empty()){
        strip.pop_back();
    }

    while(!stripSum.empty()){
        stripSum.pop_back();
    }

    for(i = 0; i < stripX.size(); ++i){
        strip.push_back(stripX[i]);
        stripSum.push_back(stripSumX[i]);
    }
    for(i = 0; i < stripY.size(); ++i){
        strip.push_back(stripY[i]);
        stripSum.push_back(stripSumY[i]);
    }

    cout << "Local search: " << strip.size() << " strips\n\n";
    cout << "Strips (X and Y combined):\n";
    for(i = 0; i < strip.size(); ++i){
        cout << "Strip " << i << ": " ;
        for(j = 0; j < strip[i].size(); ++j){
            cout << strip[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "StripSum:\n";
    for(i = 0; i < stripSum.size(); ++i){
        cout << stripSum[i] << " ";
    }
    cout << endl << endl;

}

void MBAHRA(int swapType, int moveType, int &feasible, int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores,
            vector<vector<int> > &boxWidths, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY){

    feasible = 0;
    int i, j, k;
    int threshold = 70;
    int vacant = 999;
    vector<int> scoresX;
    vector<int> orderX;
    vector<int> originalX;
    vector<int> finalX;
    vector<int> scoresY;
    vector<int> orderY;
    vector<int> originalY;
    vector<int> finalY;

    //region Creating scoresX vector
    /**FOR NEW STRIPX PAIRPAIR**/
    if(swapType == 1) {
        if(moveType == 11){
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
            originalX.push_back(stripY[j1][d1]);
            originalX.push_back(stripY[j1][d1+1]);
            goto Part2;
        }
        else {
            for (k = 0; k < stripX[i1].size(); ++k) {
                if (k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);

            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(allScores[stripY[j1][d1]]);
            originalX.push_back(stripY[j1][d1]);
            scoresX.push_back(allScores[stripY[j1][d1 + 1]]);
            originalX.push_back(stripY[j1][d1 + 1]);
            scoresX.push_back(70);
            scoresX.push_back(70);
        }
    }
        /**FOR NEW STRIPX PAIRSIN**/
    else if(swapType == 2){
        if(moveType == 21){ //stripX[i] contains 2 boxes, stripY[j] contains more than 1 box, MBAHRA on stripY[j] only
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
            goto Part2;
        }
        else {
            for (k = 0; k < stripX[i1].size(); ++k) {
                if (k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(70);
            scoresX.push_back(70);
        }
    }
        /**FOR NEW STRIPX SINSIN**/
    else if(swapType == 3){
        if(moveType == 31){ //stripX[i] contains 1 box, stripY[j] contains more than 1 box, MBAHRA on stripY[j] only
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
            goto Part2;
        }
        else {
            for (k = 0; k < stripX[i1].size(); ++k) {
                if (k == a1 || k == a1 + 1) {
                    continue;
                }
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(70);
            scoresX.push_back(70);
        }
    }
        /**FOR NEW STRIPX MOVE**/
    else if(swapType == 4){
        for(k = 0; k < stripX[i1].size(); ++k){
            scoresX.push_back(allScores[stripX[i1][k]]);
            originalX.push_back(stripX[i1][k]);
        }
        scoresX.push_back(allScores[stripY[j1][c1]]);
        originalX.push_back(stripY[j1][c1]);
        scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
        originalX.push_back(stripY[j1][c1 + 1]);
        scoresX.push_back(70);
        scoresX.push_back(70);
    }
    //endregion

    if(swapType == 1 || swapType == 2 || swapType == 3 || swapType == 4) {
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
            goto Part2;

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
            goto Part2;
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

                goto Part2;
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

    /*********************************************************************************************************/

    Part2:
    if(moveType == 41){
        feasible = 1;
        goto End;
    }

    //region Creating scoresY vector
    /**FOR NEW STRIPY PAIRPAIR**/
    if(swapType == 1) {
        if(moveType == 12){
            originalY.push_back(stripX[i1][a1]);
            originalY.push_back(stripX[i1][a1+1]);
            originalY.push_back(stripX[i1][b1]);
            originalY.push_back(stripX[i1][b1+1]);
            feasible = 1;
            goto End;
        }
        else {
            for (k = 0; k < stripY[j1].size(); ++k) {
                if (k == c1 || k == c1 + 1 || k == d1 || k == d1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);

            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            scoresY.push_back(allScores[stripX[i1][b1]]);
            originalY.push_back(stripX[i1][b1]);
            scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
            originalY.push_back(stripX[i1][b1 + 1]);
            scoresY.push_back(70);
            scoresY.push_back(70);
        }
    }
        /**FOR NEW STRIPY PAIRSIN**/
    else if(swapType == 2){
        if(moveType == 22){
            originalY.push_back(stripX[i1][a1]);
            originalY.push_back(stripX[i1][a1+1]);
            originalY.push_back(stripX[i1][b1]);
            originalY.push_back(stripX[i1][b1+1]);
            feasible = 1;
            goto End;
        }
        else {
            for (k = 0; k < stripY[j1].size(); ++k) {
                if (k == c1 || k == c1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);

            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            scoresY.push_back(allScores[stripX[i1][b1]]);
            originalY.push_back(stripX[i1][b1]);
            scoresY.push_back(allScores[stripX[i1][b1 + 1]]);
            originalY.push_back(stripX[i1][b1 + 1]);
            scoresY.push_back(70);
            scoresY.push_back(70);
        }
    }
        /**FOR NEW STRIPY SINSIN**/
    else if(swapType == 3){
        if(moveType == 32){
            originalY.push_back(stripX[i1][a1]);
            originalY.push_back(stripX[i1][a1+1]);
            feasible = 1;
            goto End;
        }
        else {
            for (k = 0; k < stripY[j1].size(); ++k) {
                if (k == c1 || k == c1 + 1) {
                    continue;
                }
                scoresY.push_back(allScores[stripY[j1][k]]);
                originalY.push_back(stripY[j1][k]);

            }
            scoresY.push_back(allScores[stripX[i1][a1]]);
            originalY.push_back(stripX[i1][a1]);
            scoresY.push_back(allScores[stripX[i1][a1 + 1]]);
            originalY.push_back(stripX[i1][a1 + 1]);
            scoresY.push_back(70);
            scoresY.push_back(70);
        }
    }
        /**FOR NEW STRIPY MOVE**/
    else if(swapType == 4){
        for(k = 0; k < stripY[j1].size(); ++k){
            if(k == c1 || k == c1 + 1){
                continue;
            }
            scoresY.push_back(allScores[stripY[j1][k]]);
            originalY.push_back(stripY[j1][k]);
        }
        scoresY.push_back(70);
        scoresY.push_back(70);
    }
    //endregion

    if(swapType == 1 || swapType == 2 || swapType == 3 || swapType == 4) {
        //region Variables
        int nScoresY = scoresY.size();
        int nBoxY = scoresY.size() / 2;
        int nCompY = (nBoxY + (nBoxY % 2)) / 2;
        vector<int> invOrderY(nScoresY);
        vector<vector<int> > adjMatY(nScoresY, vector<int>(nScoresY, 0));
        vector<int> matesY(nScoresY, 0);

        //MTGMA
        int lastMatchY = vacant;
        int mateMatchY = vacant;
        int vacantFlagY = 0;
        int matchSizeY = 0;
        vector<int> cycleVertexY(nScoresY, 1);
        vector<int> matchListY(nScoresY, vacant);

        //MIS
        int numCyclesY = 0;
        int smallestVertexY;
        int currentVertexY;
        vector<vector<int> > mateInducedY;
        vector<int> lengthMateInducedY;
        vector<int> tempMISY;
        vector<int> checkedY(nScoresY, 0);

        //FCA
        int qstarY;
        int numEdgesY;
        vector<vector<int> > SY(nCompY, vector<int>(nCompY, 0));
        vector<vector<int> > TY;
        vector<int> edgeY;
        vector<int> tY;

        //PatchGraph
        int q, u, v, saveY, SSumY, SqIntSY;
        int fullY = vacant;
        vector<int> SSetY;
        vector<int> tempPGY;
        vector<int> patchVertexY(nScoresY, vacant);
        vector<vector<int> > TpatchY;

        //MakePath
        vector<int> completePathY;
        vector<int> fullCycleY;

        //endregion

        //region Creating Instance
        for (k = 0; k < nScoresY; ++k) {
            orderY.push_back(k);
        }

        /*cout << "Scores:\n";
        for (k = 0; k < nScoresY; ++k) {
            cout << scoresY[k] << " ";
        }
        cout << endl << endl;*/

        /*cout << "Order:\n";
        for (k = 0; k < nScoresY; ++k) {
            cout << orderY[k] << " ";
        }
        cout << endl << endl;*/

        for (i = 1; i < nScoresY; ++i) {
            for (j = i - 1; j >= 0; --j) {
                if (scoresY[i] < scoresY[orderY[j]]) {
                    orderY[j + 1] = orderY[j];
                    orderY[j] = i;
                }
            }
        }

        /*cout << "Order:\n";
        for (k = 0; k < nScoresY; ++k) {
            cout << orderY[k] << " ";
        }
        cout << endl << endl;*/


        for (k = 0; k < nScoresY; ++k) {
            invOrderY[orderY[k]] = k;
        }

        /*cout << "Inverse Order:\n";
        for (k = 0; k < nScoresY; ++k) {
            cout << invOrderY[k] << " ";
        }
        cout << endl << endl;*/

        for (i = 0; i < nScoresY - 1; i += 2) {
            adjMatY[invOrderY[i]][invOrderY[i + 1]] = 2;
            adjMatY[invOrderY[i + 1]][invOrderY[i]] = 2;
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

        /*cout << "Adjacency Matrix\n";
        for (i = 0; i < nScoresY; ++i) {
            for (j = 0; j < nScoresY; ++j) {
                cout << adjMatY[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl << endl;*/


        for (i = 0; i < nScoresY; ++i) {
            for (j = 0; j < nScoresY; ++j) {
                if (adjMatY[i][j] == 2) {
                    matesY[i] = j;
                    break;
                }
            }
        }
        /*cout << "Mates Vector:\n";
        for (i = 0; i < matesY.size(); ++i) {
            cout << matesY[i] << " ";
        }
        cout << endl << endl;*/
        //endregion

        //region MTGMA
        /**MTGMA**/
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
                        && (mateMatchY >
                            i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
                        && (scoresY[lastMatchY] + scoresY[mateMatchY] >=
                            threshold)) { //can mate be matched with last matched vertex?
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
        if (matchSizeY < nBoxY) {
            feasible = 0;
            goto End;
        }

        //region MIS
        /**MIS**/
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
        if (lengthMateInducedY[0] == nScoresY) {
            for (j = 0; j < mateInducedY[0].size(); ++j) {
                fullCycleY.push_back(mateInducedY[0][j]);
            }
            //region MAKEPATH
            /**MAKEPATH FUNCTION**/
            for (i = 0; i < fullCycleY.size() - 1; ++i) {
                if ((fullCycleY[i] == nScoresY - 1 && fullCycleY[i + 1] == nScoresY - 2) ||
                    (fullCycleY[i] == nScoresY - 2 && fullCycleY[i + 1] == nScoresY - 1)) {
                    if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                        for (j = 2; j < fullCycleY.size(); ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }

                    else if (i == fullCycleY.size() -
                                  2) { //if the dominating vertices are at the end of the fullCycle vector
                        for (j = 0; j < fullCycleY.size() - 2; ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }
                    else { //if the dominating vertices are in the middle of the fullCycle vector
                        for (j = i + 2; j < fullCycleY.size(); ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;

                    }

                }
            }
            for (i = 0; i < completePathY.size(); ++i) {
                finalY.push_back(originalY[orderY[completePathY[i]]]);
            }

            feasible = 1;
            goto End;
            //END MAKE PATH FUNCTION
            //endregion

        }

        //region FCA
        /**FCA**/
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
            while (k < numEdgesY - 2 && (adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] != 1 ||
                                         cycleVertexY[edgeY[k]] == cycleVertexY[edgeY[k + 1]])) {
                ++k;
            }
            if (adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] == 1 &&
                cycleVertexY[edgeY[k]] != cycleVertexY[edgeY[k + 1]]) {
                ++qstarY;
                tY.push_back(edgeY[k]);
                SY[qstarY][cycleVertexY[edgeY[k]]] = 1;
                while (k < numEdgesY - 1 && adjMatY[edgeY[k]][matchListY[edgeY[k + 1]]] == 1 &&
                       SY[qstarY][cycleVertexY[edgeY[k + 1]]] == 0) { //add more edges to current T-cycle
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
        if (qstarY == -1) {
            //INFEASIBLE, EXIT
            feasible = 0;
            goto End;
        }

        //region PATCHGRAPH
        /**PATCHGRAPH**/
        vector<int> QSetY(qstarY, 0);
        vector<int> patchCycleY(qstarY, vacant);
        for (i = 0; i < TY.size(); ++i) {
            if (TY[i].size() == numCyclesY) {
                fullY = i;
                break;
            }
        }

        if (fullY != vacant) {
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
                if (mateInducedY[cycleVertexY[TY[fullY][v]]][saveY + 1] == TY[fullY][v]) {
                    for (i = saveY + 1; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                    for (i = mateInducedY[cycleVertexY[TY[fullY][v]]].size();
                         i-- > saveY + 1;) { //from end of cycle to element at position save+1
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                }

                    //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
                    //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
                else if (mateInducedY[cycleVertexY[TY[fullY][v]]][saveY - 1] == TY[fullY][v]) {
                    for (i = saveY; i < mateInducedY[cycleVertexY[TY[fullY][v]]].size(); ++i) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                    for (i = 0; i < saveY; ++i) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                }

                    //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
                    //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
                else if (saveY == 0 &&
                         mateInducedY[cycleVertexY[TY[fullY][v]]][mateInducedY[cycleVertexY[TY[fullY][v]]].size() -
                                                                  1] == TY[fullY][v]) {
                    for (i = 0; i < mateInducedY[cycleVertexY[TY[fullY][v]]].size(); ++i) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                }

                    //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
                    //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
                else if (saveY == mateInducedY[cycleVertexY[TY[fullY][v]]].size() - 1 &&
                         mateInducedY[cycleVertexY[TY[fullY][v]]][0] == TY[fullY][v]) {
                    for (i = mateInducedY[cycleVertexY[TY[fullY][v]]].size(); i-- > 0;) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TY[fullY][v]]][i]);
                    }
                }
                //END ONETCYCLEPATCH FUNCTION
                //endregion
            }

            //region MAKEPATH
            /**MAKEPATH FUNCTION**/
            for (i = 0; i < fullCycleY.size() - 1; ++i) {
                if ((fullCycleY[i] == nScoresY - 1 && fullCycleY[i + 1] == nScoresY - 2) ||
                    (fullCycleY[i] == nScoresY - 2 && fullCycleY[i + 1] == nScoresY - 1)) {
                    if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                        for (j = 2; j < fullCycleY.size(); ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }

                    else if (i == fullCycleY.size() -
                                  2) { //if the dominating vertices are at the end of the fullCycle vector
                        for (j = 0; j < fullCycleY.size() - 2; ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;
                    }
                    else { //if the dominating vertices are in the middle of the fullCycle vector
                        for (j = i + 2; j < fullCycleY.size(); ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePathY.push_back(fullCycleY[j]);
                        }
                        break;

                    }

                }
            }
            for (i = 0; i < completePathY.size(); ++i) {
                finalY.push_back(originalY[orderY[completePathY[i]]]);
            }

            feasible = 1;
            goto End;
            //END MAKE PATH FUNCTION
            //endregion

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
                if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == TpatchY[u][v]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = 0; i < saveY; ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == matchListY[TpatchY[u][v]]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = 0; i < saveY; ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == matchListY[TpatchY[u][v]]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = saveY; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size();
                             i-- > saveY + 1;) { //from end of cycle to element at position save+1
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == TpatchY[u][v]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = saveY; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            x = 1;
                            break;
                        }
                    }
                    if (x == 0) {
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size();
                             i-- > saveY + 1;) { //from end of cycle to element at position save+1
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }
                }

                    //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                    //the last element in the cycle is Tpatch[u][v]
                else if (saveY == 0 &&
                         mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size() -
                                                                   1] == TpatchY[u][v]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                    //the first element in the cycle is Tpatch[u][v]
                else if (saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1 &&
                         mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == TpatchY[u][v]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1; i-- > 0;) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                    //and the last element in the cycle is matchList[Tpatch[u][v]]
                else if (saveY == 0 &&
                         mateInducedY[cycleVertexY[TpatchY[u][v]]][mateInducedY[cycleVertexY[TpatchY[u][v]]].size() -
                                                                   1] == matchListY[TpatchY[u][v]]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                            break;
                        }
                    }
                }

                    //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                    //the first element in the cycle is matchList[Tpatch[u][v]]
                else if (saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1 &&
                         mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == matchListY[TpatchY[u][v]]) {
                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                    for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1; i-- > 0;) {
                        if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                            fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                        }
                        else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
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
                            }
                            else {
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
                            }
                            else {
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
                    if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == TpatchY[u][v]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = 0; i < saveY; ++i) {
                                if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                }
                                else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                    u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
                    else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY - 1] == matchListY[TpatchY[u][v]]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = saveY + 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = 0; i < saveY; ++i) {
                                if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                }
                                else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                    u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
                    else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == matchListY[TpatchY[u][v]]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = saveY; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size();
                                 i-- > saveY + 1;) { //from end of cycle to element at position save+1
                                if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                }
                                else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                    u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
                    else if (mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY + 1] == TpatchY[u][v]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = saveY; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                x = 1;
                                break;
                            }
                        }
                        if (x == 0) {
                            for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size();
                                 i-- > saveY + 1;) { //from end of cycle to element at position save+1
                                if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                }
                                else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                    fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                    u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                    break;
                                }
                            }
                        }
                    }

                        //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
                        //the last element in the cycle is Tpatch[u][v]
                    else if (saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][
                                                   mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1] ==
                                           TpatchY[u][v]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
                        //the first element in the cycle is Tpatch[u][v]
                    else if (saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1 &&
                             mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == TpatchY[u][v]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1; i-- > 0;) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
                        //and the last element in the cycle is matchList[Tpatch[u][v]]
                    else if (saveY == 0 && mateInducedY[cycleVertexY[TpatchY[u][v]]][
                                                   mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1] ==
                                           matchListY[TpatchY[u][v]]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = 1; i < mateInducedY[cycleVertexY[TpatchY[u][v]]].size(); ++i) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                                u = patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]];
                                break;
                            }
                        }
                    }

                        //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
                        //the first element in the cycle is matchList[Tpatch[u][v]]
                    else if (saveY == mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1 &&
                             mateInducedY[cycleVertexY[TpatchY[u][v]]][0] == matchListY[TpatchY[u][v]]) {
                        fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][saveY]);
                        for (i = mateInducedY[cycleVertexY[TpatchY[u][v]]].size() - 1; i-- > 0;) {
                            if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] == vacant) {
                                fullCycleY.push_back(mateInducedY[cycleVertexY[TpatchY[u][v]]][i]);
                            }
                            else if (patchVertexY[mateInducedY[cycleVertexY[TpatchY[u][v]]][i]] != vacant) {
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
                /**MAKEPATH FUNCTION**/
                for (i = 0; i < fullCycleY.size() - 1; ++i) {
                    if ((fullCycleY[i] == nScoresY - 1 && fullCycleY[i + 1] == nScoresY - 2) ||
                        (fullCycleY[i] == nScoresY - 2 && fullCycleY[i + 1] == nScoresY - 1)) {
                        if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                            for (j = 2; j < fullCycleY.size(); ++j) {
                                completePathY.push_back(fullCycleY[j]);
                            }
                            break;
                        }

                        else if (i == fullCycleY.size() -
                                      2) { //if the dominating vertices are at the end of the fullCycle vector
                            for (j = 0; j < fullCycleY.size() - 2; ++j) {
                                completePathY.push_back(fullCycleY[j]);
                            }
                            break;
                        }
                        else { //if the dominating vertices are in the middle of the fullCycle vector
                            for (j = i + 2; j < fullCycleY.size(); ++j) {
                                completePathY.push_back(fullCycleY[j]);
                            }
                            for (j = 0; j < i; ++j) {
                                completePathY.push_back(fullCycleY[j]);
                            }
                            break;

                        }

                    }
                }
                for (i = 0; i < completePathY.size(); ++i) {
                    finalY.push_back(originalY[orderY[completePathY[i]]]);
                }

                feasible = 1;
                goto End;
                //END MAKE PATH FUNCTION
                //endregion

            }
            else if (SSumY < numCyclesY) {
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

    //END MBAHRA STRIPY

    //region End
    End:
    if(feasible == 1){
        if(swapType == 1){ //PairPair
            stripSumX[i1] = stripSumX[i1] - (boxWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + boxWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                            + (boxWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + boxWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
            stripSumY[j1] = stripSumY[j1] - (boxWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + boxWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                            + (boxWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + boxWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
            if(moveType == 11){
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
            else if(moveType == 12){
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
            else {
                stripX[i1].swap(finalX);
                stripY[j1].swap(finalY);
            }
        }

        else if(swapType == 2){ //PairSin
            stripSumX[i1] = stripSumX[i1] - (boxWidths[stripX[i1][a1]][stripX[i1][a1+1]] + boxWidths[stripX[i1][b1]][stripX[i1][b1+1]])
                            + boxWidths[stripY[j1][c1]][stripY[j1][c1+1]];
            stripSumY[j1] = stripSumY[j1] - boxWidths[stripY[j1][c1]][stripY[j1][c1+1]]
                            + (boxWidths[stripX[i1][a1]][stripX[i1][a1+1]] + boxWidths[stripX[i1][b1]][stripX[i1][b1+1]]);
            if(moveType == 21){
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
            else if(moveType == 22){
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
            else {
                stripX[i1].swap(finalX);
                stripY[j1].swap(finalY);
            }
        }

        else if(swapType == 3){ //SinSin
            stripSumX[i1] = stripSumX[i1] - boxWidths[stripX[i1][a1]][stripX[i1][a1+1]] + boxWidths[stripY[j1][c1]][stripY[j1][c1+1]];
            stripSumY[j1] = stripSumY[j1] - boxWidths[stripY[j1][c1]][stripY[j1][c1+1]] + boxWidths[stripX[i1][a1]][stripX[i1][a1+1]];
            if(moveType == 31){
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
            else if(moveType == 32){
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
            else {
                stripX[i1].swap(finalX);
                stripY[j1].swap(finalY);
            }
        }

        else if(swapType == 4){
            if(moveType == 41){
                stripSumX[i1] += boxWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripX[i1].swap(finalX);
                stripY.erase(stripY.begin() + j1);
                stripSumY.erase(stripSumY.begin() + j1);
            }
            else {
                stripSumX[i1] += boxWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripSumY[j1] -= boxWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                stripX[i1].swap(finalX);
                stripY[j1].swap(finalY);
            }
        }
    }
    //endregion


}






























































































