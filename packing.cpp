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
    cout << "Initial Cost (cost from FFD): " << costFFD << endl;
    cout << "-------------------------------\n\n";

}

int initCost(int &totalCost, int maxStripWidth, vector<int> &stripSum){

    int i;
    totalCost = 0;

    for(i = 0; i < stripSum.size(); ++i){
        totalCost += pow((maxStripWidth - stripSum[i]), 2);
    }

    return totalCost;

}

void costSwap(int i, int j, int k, int l, int &totalCost, int maxStripWidth, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    int u, v;
    totalCost = totalCost - pow((maxStripWidth - stripSum[i]), 2) - pow((maxStripWidth - stripSum[k]), 2);
    stripSum[i] = stripSum[i] - boxWidths[strip[i][j]][strip[i][j+1]] + boxWidths[strip[k][l]][strip[k][l+1]];
    stripSum[k] = stripSum[k] - boxWidths[strip[k][l]][strip[k][l+1]] + boxWidths[strip[i][j]][strip[i][j+1]];
    totalCost = totalCost + pow((maxStripWidth - stripSum[i]), 2) + pow((maxStripWidth - stripSum[k]), 2);

    /*cout << "Swap\n";
    for(u = 0; u < strip.size(); ++u){
        cout << "Strip " << u << ": ";
        for(v = 0; v < strip[u].size(); ++v){
            cout << strip[u][v] << " ";
        }
        cout << endl;
    }
    cout << endl;*/
    //cout << "Current cost: " << totalCost << endl;
}

void swapBox(int i, int j, int k, int l, int &totalCost, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    if(strip[i].size() == 2){
        if(l == 0){ //CASE 9
            if(boxWidths[strip[i][j]][strip[i][j+1]] > boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[i][1]][strip[k][2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][0]);
                    swap(strip[i][1], strip[k][1]);
                }
                else if(adjMatrix[strip[i][0]][strip[k][2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][0]);
                    swap(strip[i][1], strip[k][1]);
                    swap(strip[k][0], strip[k][1]);
                }
            }
        }
        else if (l == strip[k].size() - 2){ //CASE 10
            if(boxWidths[strip[i][j]][strip[i][j+1]] > boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[i][0]][strip[k][l-1]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][l]);
                    swap(strip[i][1], strip[k][l+1]);
                }
                else if(adjMatrix[strip[i][1]][strip[k][l-1]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][l]);
                    swap(strip[i][1], strip[k][l+1]);
                    swap(strip[k][l], strip[k][l+1]);
                }
            }
        }
        else { //CASE 11: l middle vector
            if(boxWidths[strip[i][j]][strip[i][j+1]] > boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[i][0]][strip[k][l-1]] == 1 && adjMatrix[strip[i][1]][strip[k][l+2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][l]);
                    swap(strip[i][1], strip[k][l+1]);
                }
                else if(adjMatrix[strip[i][1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][0]][strip[k][l+2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][l]);
                    swap(strip[i][1], strip[k][l+1]);
                    swap(strip[k][l], strip[k][l+1]);
                }
            }
        }
    }

    else if(strip[k].size() == 2){
        if(j == 0){ //CASE 6
            if(boxWidths[strip[i][j]][strip[i][j+1]] < boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[k][1]][strip[i][2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][0]);
                    swap(strip[i][1], strip[k][1]);
                }
                else if(adjMatrix[strip[k][0]][strip[i][2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[i][0], strip[k][0]);
                    swap(strip[i][1], strip[k][1]);
                    swap(strip[i][0], strip[i][1]);
                }
            }
        }
        else if(j == strip[i].size() -2){ //CASE 7
            if(boxWidths[strip[i][j]][strip[i][j+1]] < boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[k][0]][strip[i][j-1]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[k][0], strip[i][j]);
                    swap(strip[k][1], strip[i][j+1]);
                }
                else if(adjMatrix[strip[k][1]][strip[i][j-1]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[k][0], strip[i][j]);
                    swap(strip[k][1], strip[i][j+1]);
                    swap(strip[i][j], strip[i][j+1]);
                }
            }
        }
        else{ //CASE 8: j middle vector
            if(boxWidths[strip[i][j]][strip[i][j+1]] < boxWidths[strip[k][l]][strip[k][l+1]]){
                if(adjMatrix[strip[k][0]][strip[i][j-1]] == 1 && adjMatrix[strip[k][1]][strip[i][j+2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[k][0], strip[i][j]);
                    swap(strip[k][1], strip[i][j+1]);
                }
                else if(adjMatrix[strip[k][1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][0]][strip[i][j+2]] == 1){
                    costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                    swap(strip[k][0], strip[i][j]);
                    swap(strip[k][1], strip[i][j+1]);
                    swap(strip[i][j], strip[i][j+1]);
                }
            }
        }
    }

    else if(j == 0){
        if(l == 0){ //CASE 1
            if(adjMatrix[strip[i][1]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
            }
            else if(adjMatrix[strip[i][0]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
                swap(strip[k][0], strip[k][1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
                swap(strip[i][0], strip[i][1]);
            }
            else if(adjMatrix[strip[i][0]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][1]);
                swap(strip[i][1], strip[k][0]);
            }
        }
        else if (l == strip[k].size()-2){ //CASE 2A
            if(adjMatrix[strip[i][0]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][strip[k].size()-1]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][strip[k].size()-1]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[k][l], strip[k][l+1]);

            }
            else if(adjMatrix[strip[i][0]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][l]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[i][0], strip[i][1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][l]][strip[i][2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l+1]);
                swap(strip[i][1], strip[k][l]);
            }

        }
        else { //CASE 4a: l is in the middle of the strip/vector
            if(adjMatrix[strip[k][l+1]][strip[i][2]] == 1 && adjMatrix[strip[i][0]][strip[k][l-1]] == 1 && adjMatrix[strip[i][1]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][2]] == 1 && adjMatrix[strip[i][1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][0]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[k][l], strip[k][l+1]);
            }
            else if(adjMatrix[strip[k][l]][strip[i][2]] == 1 && adjMatrix[strip[i][0]][strip[k][l-1]] == 1 && adjMatrix[strip[i][1]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[i][0], strip[i][1]);
            }
            else if(adjMatrix[strip[k][l]][strip[i][2]] == 1 && adjMatrix[strip[i][1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][0]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][0], strip[k][l+1]);
                swap(strip[i][1], strip[k][l]);
            }
        }
    }

    else if(j == strip[i].size() - 2){
        if(l == 0){ //CASE 2b
            if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
            }
            else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
                swap(strip[k][0], strip[k][1]);
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
                swap(strip[i][j], strip[i][j+1]);
            }
            else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][1]);
                swap(strip[i][j+1], strip[k][0]);
            }
        }
        else if (l == strip[k].size() - 2){ //CASE 3
            if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
                swap(strip[k][l], strip[k][l+1]);
            }
            else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
                swap(strip[i][j], strip[i][j+1]);
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l+1]);
                swap(strip[i][j+1], strip[k][l]);
            }
        }
        else { //CASE 4b: l is in the middle of the strip/vector
            if(adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
            }
            else if(adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
                swap(strip[k][l], strip[k][l+1]);
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l]);
                swap(strip[i][j+1], strip[k][l+1]);
                swap(strip[i][j], strip[i][j+1]);
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1){
                costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
                swap(strip[i][j], strip[k][l+1]);
                swap(strip[i][j+1], strip[k][l]);
            }
        }

    }

    else if (l == 0){ //CASE 4c: j is in the middle of the strip/vector
        if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1 && adjMatrix[strip[k][1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][0]);
            swap(strip[i][j+1], strip[k][1]);
        }
        else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1 && adjMatrix[strip[k][1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][0]);
            swap(strip[i][j+1], strip[k][1]);
            swap(strip[k][0], strip[k][1]);
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][0]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][0]);
            swap(strip[i][j+1], strip[k][1]);
            swap(strip[i][j], strip[i][j+1]);
        }
        else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][0]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][1]);
            swap(strip[i][j+1], strip[k][0]);
        }

    }

    else if (l == strip[k].size() - 2){ //CASE 4d: j is in the middle of the strip/vector
        if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
            swap(strip[k][l], strip[k][l+1]);
        }
        else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
            swap(strip[i][j], strip[i][j+1]);
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l+1]);
            swap(strip[i][j+1], strip[k][l]);
        }
    }

    else{ //CASE 5: both j and l are in the middle of their respective strips/vectors
        if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1
           && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
            swap(strip[k][l], strip[k][l+1]);
        }
        else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l]);
            swap(strip[i][j+1], strip[k][l+1]);
            swap(strip[i][j], strip[i][j+1]);
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            costSwap(i, j, k, l, totalCost, maxStripWidth, boxWidths, strip, stripSum);
            swap(strip[i][j], strip[k][l+1]);
            swap(strip[i][j+1], strip[k][l]);
        }
    }



}

void checkSwap(int &totalCost, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    int i, j, k, l;

    for(i = 0; i < strip.size() -1; ++i){ // ******************
        for(j = 0; j < strip[i].size() -1; j += 2){ //**************
            for(k = i+1; k < strip.size(); ++k){
                for(l = 0; l < strip[k].size()-1; l += 2){
                    //cout << "i: " << i << "\tj: " << j << "\tk: " << k << "\tl: " << l << "\tCurrent cost: " << totalCost << endl;
                    if(strip[i].size() == 2 && strip[k].size() == 2){
                        continue;
                    }
                    else if(stripSum[i] - boxWidths[strip[i][j]][strip[i][j+1]] + boxWidths[strip[k][l]][strip[k][l+1]] <= maxStripWidth
                       && stripSum[k] - boxWidths[strip[k][l]][strip[k][l+1]] + boxWidths[strip[i][j]][strip[i][j+1]] <= maxStripWidth){
                        swapBox(i, j, k, l, totalCost, maxStripWidth, adjMatrix, boxWidths, strip, stripSum);
                        //cout << "Current cost: " << totalCost << endl;

                    }//end if

                    //else move onto next pair in strip[k]

                }//end for l
                //break;
            }//end for k
            //break;
        }//end for j
        //break;
    }//end for i



    cout << "After Swapping Boxes:\n";
    cout << "Number of strips: " << strip.size() << "\n\n";

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
            cout << i << setw(11) << stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;

    cout << "Cost after swapping boxes:" << totalCost << endl;


    cout << "-------------------------\n\n";


}

void moveBox(int i, int j, int k, int l, int &moved, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    if(j == 0){
        if(l == 0 || l == strip[k].size() - 2){
            if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1) {
                strip[i].insert(strip[i].begin(), strip[k][l + 1]);
                strip[i].insert(strip[i].begin(), strip[k][l]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
            else if (adjMatrix[strip[k][l]][strip[i][j]] == 1) {
                strip[i].insert(strip[i].begin(), strip[k][l]);
                strip[i].insert(strip[i].begin(), strip[k][l + 1]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
        }
        else{
            if(adjMatrix[strip[k][l-1]][strip[k][l+2]] == 1) {
                if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1) {
                    strip[i].insert(strip[i].begin(), strip[k][l + 1]);
                    strip[i].insert(strip[i].begin(), strip[k][l]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
                else if (adjMatrix[strip[k][l]][strip[i][j]] == 1) {
                    strip[i].insert(strip[i].begin(), strip[k][l]);
                    strip[i].insert(strip[i].begin(), strip[k][l + 1]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
            }

        }


    }

    else if(j == strip[i].size() -2){
        if(l == 0 || l == strip[k].size() - 2){
            if (adjMatrix[strip[k][l]][strip[i][j+1]] == 1) {
                strip[i].push_back(strip[k][l]);
                strip[i].push_back(strip[k][l+1]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
            else if (adjMatrix[strip[k][l+1]][strip[i][j+1]] == 1) {
                strip[i].push_back(strip[k][l+1]);
                strip[i].push_back(strip[k][l]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
            else if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1) {
                strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                strip[i].insert(strip[i].begin()+j, strip[k][l]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
            else if (adjMatrix[strip[k][l]][strip[i][j]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1) {
                strip[i].insert(strip[i].begin()+j, strip[k][l]);
                strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
        }
        else{ //l in middle
            if(adjMatrix[strip[k][l-1]][strip[k][l+2]] == 1) {
                if (adjMatrix[strip[k][l]][strip[i][j+1]] == 1) {
                    strip[i].push_back(strip[k][l]);
                    strip[i].push_back(strip[k][l+1]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
                else if (adjMatrix[strip[k][l+1]][strip[i][j+1]] == 1) {
                    strip[i].push_back(strip[k][l+1]);
                    strip[i].push_back(strip[k][l]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
                else if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1) {
                    strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                    strip[i].insert(strip[i].begin()+j, strip[k][l]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
                else if (adjMatrix[strip[k][l]][strip[i][j]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1) {
                    strip[i].insert(strip[i].begin()+j, strip[k][l]);
                    strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
            }

        }

    }

    else { //j in middle
        if(l == 0 || l == strip[k].size() - 2){
            if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1) {
                strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                strip[i].insert(strip[i].begin()+j, strip[k][l]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
            else if (adjMatrix[strip[k][l]][strip[i][j]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1) {
                strip[i].insert(strip[i].begin()+j, strip[k][l]);
                strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                moved = 1;
            }
        }
        else{ //l in middle
            if(adjMatrix[strip[k][l-1]][strip[k][l+2]] == 1) {
                if (adjMatrix[strip[k][l + 1]][strip[i][j]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1) {
                    strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                    strip[i].insert(strip[i].begin()+j, strip[k][l]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
                else if (adjMatrix[strip[k][l]][strip[i][j]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1) {
                    strip[i].insert(strip[i].begin()+j, strip[k][l]);
                    strip[i].insert(strip[i].begin()+j, strip[k][l + 1]);
                    strip[k].erase(strip[k].begin() + l, strip[k].begin() + (l + 2));
                    stripSum[i] += boxWidths[strip[k][l]][strip[k][l + 1]];
                    stripSum[k] -= boxWidths[strip[k][l]][strip[k][l + 1]];
                    moved = 1;
                }
            }

        }


    }





}

void checkMove(int moved, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    int i, j, k, l, u, v;
    int maxSSi, minSSi, minu;
    int maxSSk, minSSk, maxv;
    moved = 0;

    cout << "After Moving Boxes:\n";

    //find strip k with smallest stripSum
    k = min_element(stripSum.begin(), stripSum.end()) - stripSum.begin();

    //find strip i with largest stripSum
    i = max_element(stripSum.begin(), stripSum.end()) - stripSum.begin();

    //cout << "Strip with smallest stripSum is strip " << k << ".\n";
    //cout << "Strip with largest stripSum is strip " << i << ".\n\n";
    l = 0;

    //3-> YES: CAN BOX l FIT ONTO STRIP i?
    Beginning:
    if(stripSum[i] + boxWidths[strip[k][l]][strip[k][l+1]] <= maxStripWidth){
        //cout << "point1\n";

        //4: DOES BOX l MEET MSSC ON STRIP i?
        //try moveBox for all cases of j, when moved = 1, break
        for(j = 0; j < strip[i].size() - 1; j+=2) { //go through all boxes on strip i and see if box l can fit next to/between them
            if(moved == 0) { //if box has not yet been moved, attempt to move the box
                moveBox(i, j, k, l, moved, adjMatrix, boxWidths, strip, stripSum);
                //cout << "point2\n";
            }
            else if (moved == 1){ //else if the box has finally been moved, exit for loop, no need to try to fit anywhere else
                //cout << "point3\n";
                break;
            }
        }

        //5: BOX HAS BEEN MOVED AND STRIPSUM UPDATED
        if(moved == 1){ //i.e. if box meets mssc on strip i, and a box has been moved from strip k to strip i
            moved = 0;
            //6: IS STRIP k EMPTY?
            if(strip[k].empty()){ //7: REMOVE STRIP AND STRIPSUM (is strip k empty now that a box has been removed from it?)
                strip.erase(strip.begin() + k); //delete strip k
                stripSum.erase(stripSum.begin() + k); //delete stripSum
                k = min_element(stripSum.begin(), stripSum.end()) - stripSum.begin();
                i = max_element(stripSum.begin(), stripSum.end()) - stripSum.begin();
                goto Beginning;
            }
            else{ //strip[k] not empty
                //l stays the same, k stays the same
                i = max_element(stripSum.begin(), stripSum.end()) - stripSum.begin();
                goto Beginning;
            }
        }

        else if(moved == 0){ // if box doesn't meet mssc on strip i (so box has not been moved)
            //8: is strip i the last strip to be checked that is larger than k (in terms of stripSum)?
            checkStripi:
            maxSSi = stripSum[i];
            minSSi = stripSum[k];
            minu = i;
            for(u = 0; u < stripSum.size(); ++u){
                if(u == i || u == k){
                    continue;
                }
                else if (stripSum[u] > minSSi && stripSum[u] < maxSSi){
                    minSSi = stripSum[u];
                    minu = u;
                }
            }
            if(minu == i){ // i.e. current strip is last strip checked, no more strips available
                i = max_element(stripSum.begin(), stripSum.end()) - stripSum.begin();
                if(l == strip[k].size() - 2){ //if we are on the last box in strip k
                    maxSSk = stripSum[i];
                    minSSk = stripSum[k];
                    maxv = k;
                    for(v = 0; v < stripSum.size(); ++v){
                        if(v == i || v == k){
                            continue;
                        }
                        else if(stripSum[v] > minSSk && stripSum[v] < maxSSk){
                            maxSSk = stripSum[v];
                            maxv = v;
                        }
                    }
                    if(maxv == k){
                        goto End;
                    }
                    else{ //next smallest strip k has been found
                        k = maxv;
                        l = 0;
                        goto Beginning;
                    }
                }
                else{ // if there are more boxes on strip k that have not yet been checked
                    l += 2; // move to the next box on strip k
                    goto Beginning;
                }
            }
            else { //if we have found the next largest strip i, i.e. there are more strips to be checked
                i = minu; //strip i is now the next largest strip
                goto Beginning;
            }

        }

    }

    else{ //3 -> NO: box cannot fit onto strip
        //8: is strip i the last strip to be checked that is larger than k (in terms of stripSum)?
        goto checkStripi;




    }

    End:
    //cout << "THE END\n\n";

    cout << "Number of strips: " << strip.size() << "\n\n";

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
            cout << i << setw(11) << stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;


/* We find the strip with the largest residual/smallest stripSum, this will be strip k
     * We then use the first box (box at the beginning) of this strip, and attempt to place it onto all the other strips (*)
     * We check all other strips i in order from largest stripSum to smallest stripSum
     * This is because the cost will be higher if we move a box from a strip k to a strip with a large stripSum, than to a
     * strip with a smaller stripSum
     * If the box can't fit onto a strip i, dont bother testing the mssc, move onto the next strip i+1
     * If a box can fit onto a strip i, check to see if we can move the box from strip k onto the strip i, at either
     * the beginning, middle or end of the strip i
     *
     * We then need to change the stripSum values for the two strips that have gained/lost a box
     * Then we go back and check which of the strips now has the largest residual/smallest total strip width
     *    - SHOULD BE THE SAME STRIP!
     * Continue until the strip is empty, or until none of the boxes can be moved
     *
     * Check for empty strips - if we end up with a strip that is empty, i.e. no elements/stripSum == maxStripWidth, we remove this strip entirely
     *
     * We also have to remember to remove the spaces from strip k where the box used to be before it was moved.
     *
     * (*) Is it better to start with the box at the end of strip k, or beginning? - try end first
     * What if box cannot be moved onto any strip, try next box on strip k (i.e. second box from beginning/end)?
     * What if strip k only contains one box (or multiple boxes), and that box (or none of the boxes) does/do not fit onto any of the other strips?
     *    - Need to find strip S with the second largest residual/second smallest stripSum
     * Then, if a move is performed, and the algorithm moves a box from strip S to another strip, when we find the strip with the
     * largest residual/smallest stripSum, it might find the other strip which has a box/boxes than cannot be moved
     * We need to find a way of blocking that strip from being chosen as strip k (unless a box is moved onto it?)
     *
     * USE VECTOR OF POINTERS
     */



}






