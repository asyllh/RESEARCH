/*--------------------/
ALH
packing.cpp
18/08/2017
/--------------------*/
#include <algorithm>
#include "packing.h"
using namespace std;

void packStripsSmallest(int numScores, int numBox, int maxStripWidth, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths){
//REDO THIS ALGORITHM
    // USE WITH NUMSCORES - 2
    int i, j, x, k;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    int numStrips = 0;
    vector<int> checked(numScores - 2, 0);

    strip[0].push_back(0);
    strip[0].push_back(mates[0]);
    stripSum[0] += boxWidths[0][mates[0]];
    /*for(k = 0; k < adjMatrix.size(); ++k){
        adjMatrix[k][0] = 0;
        adjMatrix[k][mates[0]] = 0;
    }*/
    checked[0] = 1;
    checked[mates[0]] = 1;
    x = 0;

    //cout << "NS: " << numScores << endl;

    for(j = 0; j < numScores - 2; ++j){
        if (checked[j] == 1) {
            continue;
        }
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(adjMatrix[strip[i].back()][j] == 1){
                    if(stripSum[i] + boxWidths[j][mates[j]] <= maxStripWidth){
                        strip[i].push_back(j);
                        strip[i].push_back(mates[j]);
                        stripSum[i] += boxWidths[j][mates[j]];
                        checked[j] = 1;
                        checked[mates[j]] = 1;
                        x = 1;
                        break;
                    }
                }
                else {
                    break;
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(j);
                strip[i].push_back(mates[j]);
                stripSum[i] += boxWidths[j][mates[j]];
                checked[j] = 1;
                checked[mates[j]] = 1;
                x = 1;
                break;
            }

        }
        if(x == 1) {
            x = 0;
            j = -1;
        }
    }

    cout << "Strips:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            for(j = 0; j < strip[i].size(); ++j){
                cout << strip[i][j] << " ";
            }
            cout << endl;
            ++numStrips;
        }
    }
    cout << endl;

    cout << "Total number of strips required: " << numStrips << endl << endl;

    cout << "Number of boxes per strip:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            cout << "Strip " << i << ": " << strip[i].size() << endl;

        }
    }
    cout << endl;

    cout << "Strip Widths(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << stripSum[i] << endl;
        }
    }
    cout << endl;

    cout << "Strip Waste(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;



}

void packStripsMIS(int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<vector<int> > &mateInduced, vector<vector<int> > &boxWidths){

    int i, j, k;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);

    i = 0;
    for(j = 0; j < mateInduced[0].size()-1; j+=2){
        if(stripSum[i] + boxWidths[mateInduced[0][j]][mateInduced[0][j+1]] <= maxStripWidth){
            strip[i].push_back(mateInduced[0][j]);
            strip[i].push_back(mateInduced[0][j+1]);
            stripSum[i] += boxWidths[mateInduced[0][j]][mateInduced[0][j+1]];
        }
        else{
            ++i;
            strip[i].push_back(mateInduced[0][j]);
            strip[i].push_back(mateInduced[0][j+1]);
            stripSum[i] += boxWidths[mateInduced[0][j]][mateInduced[0][j+1]];
        }
    }

    if(mateInduced.size() > 1){
        for(k = 1; k < mateInduced.size(); ++k){
            for(j = 0; j < mateInduced[k].size()-1; j+=2){
                for(i = 0; i < strip.size(); ++i){
                    if(!strip[i].empty()){
                        if(adjMatrix[strip[i].back()][mateInduced[k][j]] == 1){
                            if(stripSum[i] + boxWidths[mateInduced[k][j]][mateInduced[k][j+1]] <= maxStripWidth){
                                strip[i].push_back(mateInduced[k][j]);
                                strip[i].push_back(mateInduced[k][j+1]);
                                stripSum[i] += boxWidths[mateInduced[k][j]][mateInduced[k][j+1]];
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
                    else if(strip[i].empty()){
                        strip[i].push_back(mateInduced[k][j]);
                        strip[i].push_back(mateInduced[k][j+1]);
                        stripSum[i] += boxWidths[mateInduced[k][j]][mateInduced[k][j+1]];
                        break;
                    }
                }
            }
        }
    }



    cout << "Strips:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            for(j = 0; j < strip[i].size(); ++j){
                cout << strip[i][j] << " ";
            }
            cout << endl;
            ++numStrips;
        }
    }
    cout << endl;

    cout << "Total number of strips required: " << numStrips << endl << endl;

    cout << "Number of boxes per strip:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            cout << "Strip " << i << ": " << strip[i].size()/2 << endl;

        }
    }
    cout << endl;

    cout << "Strip Widths(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << stripSum[i] << endl;
        }
    }
    cout << endl;

    cout << "Strip Waste(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;

}

void weakMatchPath(int vacant, int numScores, vector<int> &matchList, vector<int> &mates){

    int i, j;
    int startPath = vacant;
    int endPath = vacant;
    int currentVertex;
    int smallestVertex;
    vector<int> checked(numScores, 0);
    vector<int> tempWMP;
    vector<vector<int> > weakPaths;


    for(i = 0; i < matchList.size(); ++i){
        if(matchList[i] == vacant){
            if(startPath == vacant){
                startPath = i;
                continue;
            }
            else if(startPath != vacant){
                endPath = i;
                break;
            }
        }
    }

    cout << "startPath: " << startPath << endl;
    cout << "endPath: " << endPath << endl;

    currentVertex = startPath;

    do {
        tempWMP.push_back(currentVertex);
        checked[currentVertex] = 1;
        tempWMP.push_back(mates[currentVertex]);
        checked[mates[currentVertex]] = 1;
        currentVertex = matchList[mates[currentVertex]];
    } while(currentVertex != mates[endPath]);

    tempWMP.push_back(currentVertex);
    checked[currentVertex] = 1;
    tempWMP.push_back(mates[currentVertex]);
    checked[mates[currentVertex]] = 1;

    weakPaths.push_back(tempWMP);
    tempWMP.clear();

    for(i = 0; i < numScores; ++i){
        if(checked[i] == 0){
            smallestVertex = i;
            break;
        }
    }

    do {
        currentVertex = smallestVertex;
        do {
            tempWMP.push_back(currentVertex);
            checked[currentVertex] = 1;
            tempWMP.push_back(mates[currentVertex]);
            checked[mates[currentVertex]] = 1;
            currentVertex = matchList[mates[currentVertex]];
        } while (currentVertex != smallestVertex);

        weakPaths.push_back(tempWMP);
        tempWMP.clear();

        for (i = 0; i < numScores; ++i) {
            if (checked[i] == 0) {
                smallestVertex = i;
                break;
            }
        }


    } while (smallestVertex != currentVertex);


    for(i = 0; i < weakPaths.size(); ++i){
        for(j = 0; j < weakPaths[i].size(); ++j){
            cout << weakPaths[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;


}

void packStripsBFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, k, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripResidual(numBox, maxStripWidth);
    vector<vector<int> > strip(numBox);
    vector<int> boxDecrease; //contains the score number of the smallest score of the boxes in decreasing order
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);
    /* i.e if mates 0-4 have boxwidth 536, mates 1-3 have BW 494, mates 2-5 have BW 940
     * then boxDecrease will contain the numbers: 2 0 1 (in that specific order)
     */

    while(boxDecrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
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

    cout << "\nBox Decrease:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;

    int minRes;
    int minResj;
    int minResk;
    int maxRes;
    int x;

    strip[0].push_back(boxDecrease[0]);
    strip[0].push_back(mates[boxDecrease[0]]);
    stripResidual[0] -= boxWidths[boxDecrease[0]][mates[boxDecrease[0]]];
    stripWidth[0].push_back(boxWidths[boxDecrease[0]][mates[boxDecrease[0]]]);


    for(i = 1; i < boxDecrease.size(); ++i){
        minRes = maxStripWidth + 1;
        for(j = 0; j < stripResidual.size(); ++j){
            if(stripResidual[j] < minRes){
                minRes = stripResidual[j];
                minResj = j;
            }
        }

        x = 0;
        do {
            if(!strip[minResj].empty()) {
                if (stripResidual[minResj] - boxWidths[boxDecrease[i]][mates[boxDecrease[i]]] >= 0) {
                    if (adjMatrix[strip[minResj].back()][boxDecrease[i]] == 1) {
                        strip[minResj].push_back(boxDecrease[i]);
                        strip[minResj].push_back(mates[boxDecrease[i]]);
                        stripResidual[minResj] -= boxWidths[boxDecrease[i]][mates[boxDecrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[boxDecrease[i]][mates[boxDecrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else if (adjMatrix[strip[minResj].back()][mates[boxDecrease[i]]] == 1) {
                        strip[minResj].push_back(mates[boxDecrease[i]]);
                        strip[minResj].push_back(boxDecrease[i]);
                        stripResidual[minResj] -= boxWidths[boxDecrease[i]][mates[boxDecrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[boxDecrease[i]][mates[boxDecrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else { //if neither score is adjacent, try next strip
                        maxRes = maxStripWidth + 1;
                        for (k = 0; k < stripResidual.size(); ++k) {
                            if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                                maxRes = stripResidual[k];
                                minResk = k;
                            }

                            else if (stripResidual[k] == minRes) {
                                if (minResj != k) {
                                    maxRes = stripResidual[k];
                                    minResk = k;
                                    break;
                                }
                            }
                        }
                        minRes = maxRes;
                        minResj = minResk;
                        continue;
                    }
                }

                else { //if box doesn't fit in strip, try next strip
                    maxRes = maxStripWidth + 1;
                    for (k = 0; k < stripResidual.size(); ++k) {
                        if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                            maxRes = stripResidual[k];
                            minResk = k;
                        }
                        else if (stripResidual[k] == minRes) {
                            if (minResj != k) {
                                maxRes = stripResidual[k];
                                minResk = k;
                                break;
                            }
                        }
                    }
                    minRes = maxRes;
                    minResj = minResk;
                    continue;
                }
            }

            else if(strip[minResj].empty()){
                strip[minResj].push_back(boxDecrease[i]);
                strip[minResj].push_back(mates[boxDecrease[i]]);
                stripResidual[minResj] -= boxWidths[boxDecrease[i]][mates[boxDecrease[i]]];
                stripWidth[minResj].push_back(boxWidths[boxDecrease[i]][mates[boxDecrease[i]]]);
                x = 1;
            }

        } while (x == 0); //end while loop


    } //end for loop

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "BFD: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips BFD (Scores):\n";
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

    cout << "Strips BFD (boxWidths):\n";
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

    cout << "Strip Widths BFD:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripResidual.size(); ++i){
        if(stripResidual[i] != maxStripWidth) {
            cout << i << "\t" << stripNumBoxes[i] <<"\t" << maxStripWidth - stripResidual[i] << "\t" << stripResidual[i] << endl;
        }
    }
    cout << endl << endl;

}

void packStripsFFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> boxDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    while(boxDecrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
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
    cout << endl;*/

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

void packStripsNFD(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> boxDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    while(boxDecrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
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

    strip[0].push_back(boxDecrease[0]);
    strip[0].push_back(mates[boxDecrease[0]]);
    stripSum[0] += boxWidths[boxDecrease[0]][mates[boxDecrease[0]]];
    stripWidth[0].push_back(boxWidths[boxDecrease[0]][mates[boxDecrease[0]]]);

    i = 0;
    for(j = 1; j < boxDecrease.size(); ++j){
        if(stripSum[i] + boxWidths[boxDecrease[j]][mates[boxDecrease[j]]] <= maxStripWidth){
            if(adjMatrix[strip[i].back()][boxDecrease[j]] == 1){
                strip[i].push_back(boxDecrease[j]);
                strip[i].push_back(mates[boxDecrease[j]]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                continue;
            }
            else if (adjMatrix[strip[i].back()][mates[boxDecrease[j]]] == 1){
                strip[i].push_back(mates[boxDecrease[j]]);
                strip[i].push_back(boxDecrease[j]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                continue;
            }
            else{
                ++i;
                strip[i].push_back(boxDecrease[j]);
                strip[i].push_back(mates[boxDecrease[j]]);
                stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
                stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
                continue;

            }
        }
        else{
            ++i;
            strip[i].push_back(boxDecrease[j]);
            strip[i].push_back(mates[boxDecrease[j]]);
            stripSum[i] += boxWidths[boxDecrease[j]][mates[boxDecrease[j]]];
            stripWidth[i].push_back(boxWidths[boxDecrease[j]][mates[boxDecrease[j]]]);
            continue;
        }
    }

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "NFD: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips NFD (scores):\n";
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

    cout << "Strips NFD (boxWidths):\n";
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

    cout << "Strip Widths NFD:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" <<stripNumBoxes[i] << "\t" << stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;





}

void packStripsBFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){
    int i, j, k, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripResidual(numBox, maxStripWidth);
    vector<vector<int> > strip(numBox);
    vector<int> boxIncrease; //contains the score number of the smallest score of the boxes in decreasing order
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);
    /* i.e if mates 0-4 have boxwidth 536, mates 1-3 have BW 494, mates 2-5 have BW 940
     * then boxDecrease will contain the numbers: 2 0 1 (in that specific order)
     */

    while(boxIncrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
        for (i = 0; i < mates.size(); ++i) {
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] < max) { //what happens if two boxes have the same width?
                max = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxIncrease.push_back(mini);
        min = max;
        max = maxBoxWidth;
    }

    cout << "\nBox Increase:\n";
    for(i = 0; i < boxIncrease.size(); ++i){
        cout << boxIncrease[i] << " ";
    }
    cout << endl << endl;

    int minRes;
    int minResj;
    int minResk;
    int maxRes;
    int x;

    strip[0].push_back(boxIncrease[0]);
    strip[0].push_back(mates[boxIncrease[0]]);
    stripResidual[0] -= boxWidths[boxIncrease[0]][mates[boxIncrease[0]]];
    stripWidth[0].push_back(boxWidths[boxIncrease[0]][mates[boxIncrease[0]]]);


    for(i = 1; i < boxIncrease.size(); ++i){
        minRes = maxStripWidth + 1;
        for(j = 0; j < stripResidual.size(); ++j){
            if(stripResidual[j] < minRes){
                minRes = stripResidual[j];
                minResj = j;
            }
        }

        x = 0;
        do {
            if(!strip[minResj].empty()) {
                if (stripResidual[minResj] - boxWidths[boxIncrease[i]][mates[boxIncrease[i]]] >= 0) {
                    if (adjMatrix[strip[minResj].back()][boxIncrease[i]] == 1) {
                        strip[minResj].push_back(boxIncrease[i]);
                        strip[minResj].push_back(mates[boxIncrease[i]]);
                        stripResidual[minResj] -= boxWidths[boxIncrease[i]][mates[boxIncrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[boxIncrease[i]][mates[boxIncrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else if (adjMatrix[strip[minResj].back()][mates[boxIncrease[i]]] == 1) {
                        strip[minResj].push_back(mates[boxIncrease[i]]);
                        strip[minResj].push_back(boxIncrease[i]);
                        stripResidual[minResj] -= boxWidths[boxIncrease[i]][mates[boxIncrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[boxIncrease[i]][mates[boxIncrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else { //if neither score is adjacent, try next strip
                        maxRes = maxStripWidth + 1;
                        for (k = 0; k < stripResidual.size(); ++k) {
                            if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                                maxRes = stripResidual[k];
                                minResk = k;
                            }

                            else if (stripResidual[k] == minRes) {
                                if (minResj != k) {
                                    maxRes = stripResidual[k];
                                    minResk = k;
                                    break;
                                }
                            }
                        }
                        minRes = maxRes;
                        minResj = minResk;
                        continue;
                    }
                }

                else { //if box doesn't fit in strip, try next strip
                    maxRes = maxStripWidth + 1;
                    for (k = 0; k < stripResidual.size(); ++k) {
                        if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                            maxRes = stripResidual[k];
                            minResk = k;
                        }
                        else if (stripResidual[k] == minRes) {
                            if (minResj != k) {
                                maxRes = stripResidual[k];
                                minResk = k;
                                break;
                            }
                        }
                    }
                    minRes = maxRes;
                    minResj = minResk;
                    continue;
                }
            }

            else if(strip[minResj].empty()){
                strip[minResj].push_back(boxIncrease[i]);
                strip[minResj].push_back(mates[boxIncrease[i]]);
                stripResidual[minResj] -= boxWidths[boxIncrease[i]][mates[boxIncrease[i]]];
                stripWidth[minResj].push_back(boxWidths[boxIncrease[i]][mates[boxIncrease[i]]]);
                x = 1;
            }

        } while (x == 0); //end while loop


    } //end for loop

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "BFI: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips BFI (Scores):\n";
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

    cout << "Strips BFI (boxWidths):\n";
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

    cout << "Strip Widths BFI:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripResidual.size(); ++i){
        if(stripResidual[i] != maxStripWidth) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" << maxStripWidth - stripResidual[i] << "\t" << stripResidual[i] << endl;
        }
    }
    cout << endl << endl;

}

void packStripsFFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> boxIncrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    while(boxIncrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
        for (i = 0; i < mates.size(); ++i) {
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] < max) { //what happens if two boxes have the same width?
                max = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxIncrease.push_back(mini);
        min = max;
        max = maxBoxWidth;
    }

    /*cout << "\nBox Increase:\n";
    for(i = 0; i < boxDecrease.size(); ++i){
        cout << boxDecrease[i] << " ";
    }
    cout << endl << endl;*/

    strip[0].push_back(boxIncrease[0]);
    strip[0].push_back(mates[boxIncrease[0]]);
    stripSum[0] += boxWidths[boxIncrease[0]][mates[boxIncrease[0]]];
    stripWidth[0].push_back(boxWidths[boxIncrease[0]][mates[boxIncrease[0]]]);


    for(j = 1; j < boxIncrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[boxIncrease[j]][mates[boxIncrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][boxIncrease[j]] == 1){
                        strip[i].push_back(boxIncrease[j]);
                        strip[i].push_back(mates[boxIncrease[j]]);
                        stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                        stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][mates[boxIncrease[j]]] == 1){
                        strip[i].push_back(mates[boxIncrease[j]]);
                        strip[i].push_back(boxIncrease[j]);
                        stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                        stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
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
                strip[i].push_back(boxIncrease[j]);
                strip[i].push_back(mates[boxIncrease[j]]);
                stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
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

    cout << "FFI: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips FFI (scores):\n";
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

    cout << "Strips FFI (boxWidths):\n";
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

    cout << "Strip Widths FFI:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" << stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;



}

void packStripsNFI(int numBox, int maxBoxWidth, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, mini;
    int min = 0;
    int max = maxBoxWidth;
    int numStrips = 0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> boxIncrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    while(boxIncrease.size() < numBox -1) { //numBox -1, doesn't include dominating vertices
        for (i = 0; i < mates.size(); ++i) {
            if (boxWidths[i][mates[i]] > min && boxWidths[i][mates[i]] < max) { //what happens if two boxes have the same width?
                max = boxWidths[i][mates[i]];
                mini = i;
            }
        }
        boxIncrease.push_back(mini);
        min = max;
        max = maxBoxWidth;
    }

    strip[0].push_back(boxIncrease[0]);
    strip[0].push_back(mates[boxIncrease[0]]);
    stripSum[0] += boxWidths[boxIncrease[0]][mates[boxIncrease[0]]];
    stripWidth[0].push_back(boxWidths[boxIncrease[0]][mates[boxIncrease[0]]]);

    i = 0;
    for(j = 1; j < boxIncrease.size(); ++j){
        if(stripSum[i] + boxWidths[boxIncrease[j]][mates[boxIncrease[j]]] <= maxStripWidth){
            if(adjMatrix[strip[i].back()][boxIncrease[j]] == 1){
                strip[i].push_back(boxIncrease[j]);
                strip[i].push_back(mates[boxIncrease[j]]);
                stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
                continue;
            }
            else if (adjMatrix[strip[i].back()][mates[boxIncrease[j]]] == 1){
                strip[i].push_back(mates[boxIncrease[j]]);
                strip[i].push_back(boxIncrease[j]);
                stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
                continue;
            }
            else{
                ++i;
                strip[i].push_back(boxIncrease[j]);
                strip[i].push_back(mates[boxIncrease[j]]);
                stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
                stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
                continue;

            }
        }
        else{
            ++i;
            strip[i].push_back(boxIncrease[j]);
            strip[i].push_back(mates[boxIncrease[j]]);
            stripSum[i] += boxWidths[boxIncrease[j]][mates[boxIncrease[j]]];
            stripWidth[i].push_back(boxWidths[boxIncrease[j]][mates[boxIncrease[j]]]);
            continue;
        }
    }

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "NFI: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips NFI (scores):\n";
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

    cout << "Strips NFI (boxWidths):\n";
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

    cout << "Strip Widths NFI:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" << stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;





}

void packStripsFFDScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> scoreDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = dummyMates.size(); i-- > 0;){
        if(dummyMates[i] != vacant){
            scoreDecrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores decreasing:\n";
    for(i = 0; i < scoreDecrease.size(); ++i){
        cout << scoreDecrease[i] << " ";
    }
    cout << endl << endl;

    strip[0].push_back(mates[scoreDecrease[0]]);
    strip[0].push_back(scoreDecrease[0]);
    stripSum[0] += boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]]);

    for(j = 1; j < scoreDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][mates[scoreDecrease[j]]] == 1){
                        strip[i].push_back(mates[scoreDecrease[j]]);
                        strip[i].push_back(scoreDecrease[j]);
                        stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                        stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][scoreDecrease[j]] == 1){
                        strip[i].push_back(scoreDecrease[j]);
                        strip[i].push_back(mates[scoreDecrease[j]]);
                        stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                        stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
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
            else if(strip[i].empty()){
                strip[i].push_back(mates[scoreDecrease[j]]);
                strip[i].push_back(scoreDecrease[j]);
                stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
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


    cout << "FFDScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips FFDScores (scores):\n";
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

    cout << "Strips FFDScores (boxWidths):\n";
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

    cout << "Strip Widths FFDScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" <<  stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;





}

void packStripsFFIScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> scoreIncrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = 0; i < dummyMates.size(); ++i){
        if(dummyMates[i] != vacant){
            scoreIncrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores Increasing:\n";
    for(i = 0; i < scoreIncrease.size(); ++i){
        cout << scoreIncrease[i] << " ";
    }
    cout << endl << endl;

    strip[0].push_back(scoreIncrease[0]);
    strip[0].push_back(mates[scoreIncrease[0]]);
    stripSum[0] += boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]]);

    for(j = 1; j < scoreIncrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]] <= maxStripWidth){
                    if(adjMatrix[strip[i].back()][scoreIncrease[j]] == 1){
                        strip[i].push_back(scoreIncrease[j]);
                        strip[i].push_back(mates[scoreIncrease[j]]);
                        stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                        stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][mates[scoreIncrease[j]]] == 1){
                        strip[i].push_back(mates[scoreIncrease[j]]);
                        strip[i].push_back(scoreIncrease[j]);
                        stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                        stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
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
            else if(strip[i].empty()){
                strip[i].push_back(scoreIncrease[j]);
                strip[i].push_back(mates[scoreIncrease[j]]);
                stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
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


    cout << "FFIScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips FFIScores (scores):\n";
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

    cout << "Strips FFIScores (boxWidths):\n";
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

    cout << "Strip Widths FFIScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" <<  stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;








}

void packStripsBFDScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, k;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripResidual(numBox, maxStripWidth);
    vector<vector<int> > strip(numBox);
    vector<int> scoreDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = dummyMates.size(); i-- > 0;){
        if(dummyMates[i] != vacant){
            scoreDecrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores decreasing:\n";
    for(i = 0; i < scoreDecrease.size(); ++i){
        cout << scoreDecrease[i] << " ";
    }
    cout << endl << endl;

    int minRes;
    int minResj;
    int minResk;
    int maxRes;
    int x;

    strip[0].push_back(mates[scoreDecrease[0]]);
    strip[0].push_back(scoreDecrease[0]);
    stripResidual[0] -= boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]]);


    for(i = 1; i < scoreDecrease.size(); ++i){
        minRes = maxStripWidth + 1;
        for(j = 0; j < stripResidual.size(); ++j){ //this loop will give us the strip with the smallest residual
            if(stripResidual[j] < minRes){
                minRes = stripResidual[j];
                minResj = j;
            }
        }

        x = 0;
        do {
            if(!strip[minResj].empty()) {
                if (stripResidual[minResj] - boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]] >= 0) {
                    if (adjMatrix[strip[minResj].back()][mates[scoreDecrease[i]]] == 1) {
                        strip[minResj].push_back(mates[scoreDecrease[i]]);
                        strip[minResj].push_back(scoreDecrease[i]);
                        stripResidual[minResj] -= boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else if (adjMatrix[strip[minResj].back()][scoreDecrease[i]] == 1) {
                        strip[minResj].push_back(scoreDecrease[i]);
                        strip[minResj].push_back(mates[scoreDecrease[i]]);
                        stripResidual[minResj] -= boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else { //if neither score is adjacent, try next strip (i.e the strip with the second smallest residual)
                        maxRes = maxStripWidth + 1;
                        for (k = 0; k < stripResidual.size(); ++k) {
                            if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                                maxRes = stripResidual[k];
                                minResk = k;
                            }

                            else if (stripResidual[k] == minRes) {
                                if (minResj != k) {
                                    maxRes = stripResidual[k];
                                    minResk = k;
                                    break;
                                }
                            }
                        }
                        minRes = maxRes;
                        minResj = minResk;
                        continue;
                    }
                }

                else { //if box doesn't fit in strip, try next strip
                    maxRes = maxStripWidth + 1;
                    for (k = 0; k < stripResidual.size(); ++k) {
                        if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                            maxRes = stripResidual[k];
                            minResk = k;
                        }
                        else if (stripResidual[k] == minRes) {
                            if (minResj != k) {
                                maxRes = stripResidual[k];
                                minResk = k;
                                break;
                            }
                        }
                    }
                    minRes = maxRes;
                    minResj = minResk;
                    continue;
                }
            }

            else if(strip[minResj].empty()){
                strip[minResj].push_back(mates[scoreDecrease[i]]);
                strip[minResj].push_back(scoreDecrease[i]);
                stripResidual[minResj] -= boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]];
                stripWidth[minResj].push_back(boxWidths[scoreDecrease[i]][mates[scoreDecrease[i]]]);
                x = 1;
            }

        } while (x == 0); //end while loop


    } //end for loop


    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "BFDScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips BFDScores (Scores):\n";
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

    cout << "Strips BFDScores (boxWidths):\n";
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

    cout << "Strip Widths BFDScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripResidual.size(); ++i){
        if(stripResidual[i] != maxStripWidth) {
            cout << i << "\t\t" << stripNumBoxes[i] <<"\t\t" << maxStripWidth - stripResidual[i] << "\t\t" << stripResidual[i] << endl;
        }
    }
    cout << endl << endl;


}

void packStripsBFIScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j, k;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripResidual(numBox, maxStripWidth);
    vector<vector<int> > strip(numBox);
    vector<int> scoreIncrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = 0; i < dummyMates.size(); ++i){
        if(dummyMates[i] != vacant){
            scoreIncrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores Increasing:\n";
    for(i = 0; i < scoreIncrease.size(); ++i){
        cout << scoreIncrease[i] << " ";
    }
    cout << endl << endl;

    int minRes;
    int minResj;
    int minResk;
    int maxRes;
    int x;

    strip[0].push_back(scoreIncrease[0]);
    strip[0].push_back(mates[scoreIncrease[0]]);
    stripResidual[0] -= boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]]);


    for(i = 1; i < scoreIncrease.size(); ++i){
        minRes = maxStripWidth + 1;
        for(j = 0; j < stripResidual.size(); ++j){ //this loop will give us the strip with the smallest residual
            if(stripResidual[j] < minRes){
                minRes = stripResidual[j];
                minResj = j;
            }
        }

        x = 0;
        do {
            if(!strip[minResj].empty()) {
                if (stripResidual[minResj] - boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]] >= 0) {
                    if (adjMatrix[strip[minResj].back()][scoreIncrease[i]] == 1) {
                        strip[minResj].push_back(scoreIncrease[i]);
                        strip[minResj].push_back(mates[scoreIncrease[i]]);
                        stripResidual[minResj] -= boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else if (adjMatrix[strip[minResj].back()][mates[scoreIncrease[i]]] == 1) {
                        strip[minResj].push_back(mates[scoreIncrease[i]]);
                        strip[minResj].push_back(scoreIncrease[i]);
                        stripResidual[minResj] -= boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]];
                        stripWidth[minResj].push_back(boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]]);
                        x = 1;
                        continue;
                    }

                    else { //if neither score is adjacent, try next strip (i.e the strip with the second smallest residual)
                        maxRes = maxStripWidth + 1;
                        for (k = 0; k < stripResidual.size(); ++k) {
                            if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                                maxRes = stripResidual[k];
                                minResk = k;
                            }

                            else if (stripResidual[k] == minRes) {
                                if (minResj != k) {
                                    maxRes = stripResidual[k];
                                    minResk = k;
                                    break;
                                }
                            }
                        }
                        minRes = maxRes;
                        minResj = minResk;
                        continue;
                    }
                }

                else { //if box doesn't fit in strip, try next strip
                    maxRes = maxStripWidth + 1;
                    for (k = 0; k < stripResidual.size(); ++k) {
                        if (stripResidual[k] > minRes && stripResidual[k] < maxRes) {
                            maxRes = stripResidual[k];
                            minResk = k;
                        }
                        else if (stripResidual[k] == minRes) {
                            if (minResj != k) {
                                maxRes = stripResidual[k];
                                minResk = k;
                                break;
                            }
                        }
                    }
                    minRes = maxRes;
                    minResj = minResk;
                    continue;
                }
            }

            else if(strip[minResj].empty()){
                strip[minResj].push_back(scoreIncrease[i]);
                strip[minResj].push_back(mates[scoreIncrease[i]]);
                stripResidual[minResj] -= boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]];
                stripWidth[minResj].push_back(boxWidths[scoreIncrease[i]][mates[scoreIncrease[i]]]);
                x = 1;
            }

        } while (x == 0); //end while loop


    } //end for loop


    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "BFIScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips BFIScores (Scores):\n";
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

    cout << "Strips BFIScores (boxWidths):\n";
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

    cout << "Strip Widths BFIScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripResidual.size(); ++i){
        if(stripResidual[i] != maxStripWidth) {
            cout << i << "\t\t" << stripNumBoxes[i] <<"\t\t" << maxStripWidth - stripResidual[i] << "\t\t" << stripResidual[i] << endl;
        }
    }
    cout << endl << endl;


}

void packStripsNFDScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> scoreDecrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = dummyMates.size(); i-- > 0;){
        if(dummyMates[i] != vacant){
            scoreDecrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores decreasing:\n";
    for(i = 0; i < scoreDecrease.size(); ++i){
        cout << scoreDecrease[i] << " ";
    }
    cout << endl << endl;

    strip[0].push_back(mates[scoreDecrease[0]]);
    strip[0].push_back(scoreDecrease[0]);
    stripSum[0] += boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreDecrease[0]][mates[scoreDecrease[0]]]);


    i = 0;
    for(j = 1; j < scoreDecrease.size(); ++j){
        if(stripSum[i] + boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]] <= maxStripWidth){
            if(adjMatrix[strip[i].back()][mates[scoreDecrease[j]]] == 1){
                strip[i].push_back(mates[scoreDecrease[j]]);
                strip[i].push_back(scoreDecrease[j]);
                stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
                continue;
            }
            else if (adjMatrix[strip[i].back()][scoreDecrease[j]] == 1){
                strip[i].push_back(scoreDecrease[j]);
                strip[i].push_back(mates[scoreDecrease[j]]);
                stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
                continue;
            }
            else{
                ++i;
                strip[i].push_back(mates[scoreDecrease[j]]);
                strip[i].push_back(scoreDecrease[j]);
                stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
                continue;

            }
        }
        else{
            ++i;
            strip[i].push_back(mates[scoreDecrease[j]]);
            strip[i].push_back(scoreDecrease[j]);
            stripSum[i] += boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]];
            stripWidth[i].push_back(boxWidths[scoreDecrease[j]][mates[scoreDecrease[j]]]);
            continue;
        }
    }

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "NFDScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips NFDScores (scores):\n";
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

    cout << "Strips NFDScores (boxWidths):\n";
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

    cout << "Strip Widths NFDScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" <<stripNumBoxes[i] << "\t" << stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;













}

void packStripsNFIScores(int vacant, int numBox, int maxStripWidth, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths){

    int i, j;
    int numStrips = 0;
    vector<int> dummyMates;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<int> scoreIncrease;
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes(numBox, 0);

    for(i = 0; i < mates.size() -2; ++i){
        dummyMates.push_back(mates[i]);
    }

    /*cout << "Dummy Mates:\n";
    for(i = 0; i < dummyMates.size(); ++i){
        cout << dummyMates[i] << " ";
    }
    cout << endl;*/

    for(i = 0; i < dummyMates.size(); ++i){
        if(dummyMates[i] != vacant){
            scoreIncrease.push_back(i);
            dummyMates[dummyMates[i]] = vacant;
        }
    }

    cout << "scores Increasing:\n";
    for(i = 0; i < scoreIncrease.size(); ++i){
        cout << scoreIncrease[i] << " ";
    }
    cout << endl << endl;

    strip[0].push_back(scoreIncrease[0]);
    strip[0].push_back(mates[scoreIncrease[0]]);
    stripSum[0] += boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]];
    stripWidth[0].push_back(boxWidths[scoreIncrease[0]][mates[scoreIncrease[0]]]);

    i = 0;
    for(j = 1; j < scoreIncrease.size(); ++j){
        if(stripSum[i] + boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]] <= maxStripWidth){
            if(adjMatrix[strip[i].back()][scoreIncrease[j]] == 1){
                strip[i].push_back(scoreIncrease[j]);
                strip[i].push_back(mates[scoreIncrease[j]]);
                stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
                continue;
            }
            else if (adjMatrix[strip[i].back()][mates[scoreIncrease[j]]] == 1){
                strip[i].push_back(mates[scoreIncrease[j]]);
                strip[i].push_back(scoreIncrease[j]);
                stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
                continue;
            }
            else{
                ++i;
                strip[i].push_back(scoreIncrease[j]);
                strip[i].push_back(mates[scoreIncrease[j]]);
                stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
                stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
                continue;

            }
        }
        else{
            ++i;
            strip[i].push_back(scoreIncrease[j]);
            strip[i].push_back(mates[scoreIncrease[j]]);
            stripSum[i] += boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]];
            stripWidth[i].push_back(boxWidths[scoreIncrease[j]][mates[scoreIncrease[j]]]);
            continue;
        }
    }

    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            stripNumBoxes[i] = strip[i].size() / 2;
            ++numStrips;
        }
    }

    cout << "NFIScores: " << numStrips << " strips\n-----------------------\n";

    cout << "Strips NFIScores (scores):\n";
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

    cout << "Strips NFIScores (boxWidths):\n";
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

    cout << "Strip Widths NFIScores:\n";
    cout << "Strip\t#Boxes\tWidth\tResidual\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] !=0) {
            cout << i << "\t" << stripNumBoxes[i] << "\t" << stripSum[i] << "\t" << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl << endl;


}

void packStripsSmallestSearch(int numScores, int numBox, int maxStripWidth, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths){

    // USE WITH NUMSCORES - 2
    int i, j, x, k, u, v;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    int numStrips = 0;
    vector<int> checked(numScores - 2, 0);

    strip[0].push_back(0);
    strip[0].push_back(mates[0]);
    stripSum[0] += boxWidths[0][mates[0]];
    /*for(k = 0; k < adjMatrix.size(); ++k){
        adjMatrix[k][0] = 0;
        adjMatrix[k][mates[0]] = 0;
    }*/
    checked[0] = 1;
    checked[mates[0]] = 1;
    x = 0;
    i = 0;

    do {
        for (j = 0; j < numScores - 2; ++j) {
            if (checked[j] == 1) {
                continue;
            }
            if (!strip[i].empty()) {
                if (adjMatrix[strip[i].back()][j] == 1) {
                    if (stripSum[i] + boxWidths[j][mates[j]] <= maxStripWidth) {
                        strip[i].push_back(j);
                        strip[i].push_back(mates[j]);
                        stripSum[i] += boxWidths[j][mates[j]];
                        checked[j] = 1;
                        checked[mates[j]] = 1;
                        /*for (k = 0; k < adjMatrix.size(); ++k) {
                            adjMatrix[k][j] = 0;
                            adjMatrix[k][mates[j]] = 0;
                        }*/
                        x = 1;
                    }

                }
            }
            else if (strip[i].empty()) {
                strip[i].push_back(j);
                strip[i].push_back(mates[j]);
                stripSum[i] += boxWidths[j][mates[j]];
                checked[j] = 1;
                checked[mates[j]] = 1;
                /*for (k = 0; k < adjMatrix.size(); ++k) {
                    adjMatrix[k][j] = 0;
                    adjMatrix[k][mates[j]] = 0;
                }*/
                x = 1;

            }

            if(x == 1){
                x = 0;
                j = -1;
            }

        }

        v = -1;
        for (u = 0; u < checked.size(); ++u) {
            if (checked[u] == 0) {
                v = u;
                break;
            }
        }
        ++i;

    }while(v != -1);


    cout << "Strips:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            for(j = 0; j < strip[i].size(); ++j){
                cout << strip[i][j] << " ";
            }
            cout << endl;
            ++numStrips;
        }
    }
    cout << endl;

    cout << "Total number of strips required: " << numStrips << endl << endl;

    cout << "Number of boxes per strip:\n";
    for(i = 0; i < strip.size(); ++i){
        if(!strip[i].empty()){
            cout << "Strip " << i << ": " << strip[i].size() << endl;

        }
    }
    cout << endl;

    cout << "Strip Widths(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << stripSum[i] << endl;
        }
    }
    cout << endl;

    cout << "Strip Waste(mm):\n";
    for(i = 0; i < stripSum.size(); ++i){
        if(stripSum[i] != 0){
            cout << "Strip " << i << ": " << maxStripWidth - stripSum[i] << endl;
        }
    }
    cout << endl;



}























