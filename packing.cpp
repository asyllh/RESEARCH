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

    int i, j, k, nScores, nBox;
    int threshold = 70;
    vector<int> scores;
    vector<int> order;

    for(k = 0; k < stripX[i1].size(); ++k){
        if(k == a1 || k == a1 + 1 || k == b1 || k == b1 + 1){
            continue;
        }
        scores.push_back(allScores[stripX[i1][k]]);
    }
    scores.push_back(allScores[stripY[j1][c1]]);
    scores.push_back(allScores[stripY[j1][c1 + 1]]);
    scores.push_back(allScores[stripY[j1][d1]]);
    scores.push_back(allScores[stripY[j1][d1 + 1]]);
    scores.push_back(70);
    scores.push_back(70);

    nScores = scores.size();
    nBox = scores.size() /2;
    vector<int> invOrder(nScores);
    vector<vector<int> > adjMat(nScores, vector<int>(nScores, 0));

    for(k = 0; k < nScores; ++k){
        order.push_back(k);
    }

    cout << "Scores:\n";
    for(k = 0; k < nScores; ++k){
        cout << scores[k] << " ";
    }
    cout << endl << endl;

    cout << "Order:\n";
    for(k = 0; k < nScores; ++k){
        cout << order[k] << " ";
    }
    cout << endl << endl;

    for(i = 1; i < nScores; ++i){
        for(j = i-1; j >= 0; --j){
            if(scores[i] < scores[order[j]]){
                order[j+1] = order[j];
                order[j] = i;
            }
        }
    }

    cout << "Order:\n";
    for(k = 0; k < nScores; ++k){
        cout << order[k] << " ";
    }
    cout << endl << endl;


    for(k = 0; k < nScores; ++k){
        invOrder[order[k]] = k;
    }

    cout << "Inverse Order:\n";
    for(k = 0; k < nScores; ++k){
        cout << invOrder[k] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < nScores; i+=2){
        adjMat[invOrder[i]][invOrder[i+1]] = 2;
        adjMat[invOrder[i+1]][invOrder[i]] = 2;
    }

    cout << "Adjacency Matrix\n";
    for(i = 0; i < nScores; ++i){
        for(j = 0; j < nScores; ++j){
            cout << adjMat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;



















}



























































































