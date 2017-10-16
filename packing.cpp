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

void packStripsFFD(int numBox, int maxBoxWidth, int maxStripWidth, double totalBoxWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<int> &stripNumBoxes, vector<vector<int> > &strip, vector<vector<int> > &stripWidth){

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

    cout << "Strips FFD (boxWidths):\n";
    for(i = 0; i < stripWidth.size(); ++i){
        cout << "Strip " << i << ": ";
        for(j = 0; j < stripWidth[i].size(); ++j){
            cout << stripWidth[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Strip" << setw(8) << "#Boxes" << setw(8) << "Width" << setw(12) << "Residual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << setw(9) << stripNumBoxes[i] << setw(10) <<  stripSum[i] << setw(9) << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;

}

void checkSwap(int i, int j, int k, int l, vector<vector<int> > &adjMatrix, vector<vector<int> > &strip){

    int u, v;

    if(j == 0){
        if(l == 0){ //CASE 1
            if(adjMatrix[strip[i][1]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);

                /*for(u = 0; u < strip.size(); ++u){
                    cout << "Strip " << u << ": ";
                    for(v = 0; v < strip[u].size(); ++v){
                        cout << strip[u][v] << " ";
                    }
                    cout << endl;
                }
                cout << endl;*/
            }
            else if(adjMatrix[strip[i][0]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
                swap(strip[k][0], strip[k][1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
                swap(strip[i][0], strip[i][1]);
            }
            else if(adjMatrix[strip[i][0]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][0]);
                swap(strip[i][1], strip[k][1]);
                swap(strip[i][0], strip[i][1]);
                swap(strip[k][0], strip[k][1]);
            }
        }
        else if (l == strip[k].size()-2){ //CASE 2A
            if(adjMatrix[strip[i][0]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][strip[k].size()-1]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][strip[k].size()-1]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[k][l], strip[k][l+1]);

            }
            else if(adjMatrix[strip[i][0]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][l]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[i][0], strip[i][1]);
            }
            else if(adjMatrix[strip[i][1]][strip[k][strip[k].size()-3]] == 1 && adjMatrix[strip[k][l]][strip[i][2]] == 1){
                swap(strip[i][0], strip[k][l]);
                swap(strip[i][1], strip[k][l+1]);
                swap(strip[i][0], strip[i][1]);
                swap(strip[k][l], strip[k][l+1]);
            }

        }
        else { //CASE 4a: l is in the middle of the strip/vector
            if(adjMatrix[strip[k][l+1]][strip[i][2]] == 1 && adjMatrix[strip[i][0]][strip[k][l-1]] == 1 && adjMatrix[strip[i][1]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][2]] == 1 && adjMatrix[strip[i][1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][0]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l]][strip[i][2]] == 1 && adjMatrix[strip[i][0]][strip[k][l-1]] == 1 && adjMatrix[strip[i][1]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l]][strip[i][2]] == 1 && adjMatrix[strip[i][1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][0]][strip[k][l+2]] == 1){
                //perform swap
            }
        }
    }

    else if(j == strip[i].size() - 2){
        if(l == 0){ //CASE 2b
            if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1){
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
            }
            else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1){
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
                swap(strip[k][0], strip[k][1]);
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
                swap(strip[i][j], strip[i][j+1]);
            }
            else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                swap(strip[i][j], strip[k][0]);
                swap(strip[i][j+1], strip[k][1]);
                swap(strip[k][0], strip[k][1]);
                swap(strip[i][j], strip[i][j+1]);
            }



        }
        else if (l == strip[k].size() - 2){ //CASE 3
            if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1){
                //perform swap
            }


        }
        else { //CASE 4b: l is in the middle of the strip/vector
            if(adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1){
                //perform swap
            }
            else if(adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1){
                //perform swap
            }
        }

    }

    else if (l == 0){ //CASE 4c: j is in the middle of the strip/vector
        if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1 && adjMatrix[strip[k][1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][0]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][0]][strip[i][j-1]] == 1 && adjMatrix[strip[k][1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j]][strip[k][2]] == 1 && adjMatrix[strip[k][1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][0]][strip[i][j+2]] == 1){
            //perform swap
        }

    }

    else if (l == strip[k].size() - 2){ //CASE 4d: j is in the middle of the strip/vector
        if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            //perform swap
        }
    }

    else{ //CASE 5: both j and l are in the middle of their respective strips/vectors
        if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1
           && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l+1]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j+1]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            //perform swap
        }
        else if(adjMatrix[strip[i][j+1]][strip[k][l-1]] == 1 && adjMatrix[strip[i][j]][strip[k][l+2]] == 1
                && adjMatrix[strip[k][l+1]][strip[i][j-1]] == 1 && adjMatrix[strip[k][l]][strip[i][j+2]] == 1){
            //perform swap
        }
    }



}

void swapBox(int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<vector<int> > &strip, vector<int> &stripSum){

    int i, j, k, l;

    for(i = 0; i < strip.size() - 1; ++i){
        for(j = 0; j < strip[i].size() -1; j += 2){
            for(k = i+1; k < 3; ++k){
                for(l = 0; l < strip[k].size()-1; l += 2){
                    if(stripSum[i] - boxWidths[strip[i][j]][strip[i][j+1]] + boxWidths[strip[k][l]][strip[k][l+1]] <= maxStripWidth
                       && stripSum[k] - boxWidths[strip[k][l]][strip[k][l+1]] + boxWidths[strip[i][j]][strip[i][j+1]] <= maxStripWidth){
                        checkSwap(i, j, k, l, adjMatrix, strip);
                        //break;





                    }//end if

                    //else move onto next pair in strip[k]

                }//end for l
                //break;
            }//end for k
            //break;
        }//end for j
        //break;
    }//end for i

    cout << "END\n";

}



