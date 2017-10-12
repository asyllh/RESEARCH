/*--------------/
ALH
packing.cpp
12/10/17
/--------------*/
#include <algorithm>
#include <cmath>
#include "packing.h"
using namespace std;

int lowerBound(double totalBoxWidth, int maxStripWidth){
    int lBound;

    lBound = ceil(totalBoxWidth/maxStripWidth);
    return lBound;
}

void packStripsFFD(int numBox, int maxBoxWidth, int maxStripWidth, double totalBoxWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> boxDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    while(boxDecrease.size() < numBox) { //numBox -1, doesn't include dominating vertices
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

    cout << "Box decrease:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;

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
                    else{
                        continue;
                    }
                }
                else{
                    continue;
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


    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }


    cout << "FFD: " << numStrips << " strips\n-----------------------\n";

    cout << "Lower Bound: " << lowerBound(totalBoxWidth, maxStripWidth) << " strips\n---------------\n";

    cout << "Strips FFD (scores):\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            cout << "Strip " << i << ": ";
            for(j = 0; j < strip[i].size(); ++j){
                cout << strip[i][j] << " ";
            }
            cout << endl;
        }
    }
    cout << endl;

    cout << "Strips FFD (boxWidths):\n";
    for(i = 0; i < stripWidth.size(); ++i){
        if(!stripWidth[i].empty()){
            cout << "Strip " << i << ": ";
            for(j = 0; j < stripWidth[i].size(); ++j){
                cout << stripWidth[i][j] << " ";
            }
            cout << endl;
        }
    }
    cout << endl;

    cout << "Strip Widths FFD:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" <<  stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;



}

