/*--------------/
ALH
packing.cpp
Combined Program with Heuristics and EA
17/03/18
/--------------*/

#include <algorithm>
#include <iomanip>
#include "packing.h"
using namespace std;

void Swap(int &a, int &b){
    int temp = a;
    a = b;
    b = temp;
}

int LowerBound(double totalItemWidth, int stripWidth){
    int lBound;

    lBound = ceil(totalItemWidth/stripWidth);
    return lBound;
}

double Fitness(int stripWidth, vector<int> &stripSum, vector<vector<int> > &strip){

    int i;
    double total = 0.0;
    double final;

    for(i = 0; i < strip.size(); ++i){
        double a = static_cast<double>(stripSum[i]) / static_cast<double>(stripWidth);
        total += pow(a, 2);
    }

    final = total / static_cast<double>(strip.size());

    return final;

}

void Optimality(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int stripSize, int LB){

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
} //End Optimality

// FFD checking vicinal sum constraint for both sides of each item.
void MFFD(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem,
          int maxItemWidth, int stripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners,
          vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){
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
                if(stripSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= stripWidth){
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

    while(stripSum.back() == 0){
        stripSum.pop_back();
        strip.pop_back();
    }

    //int stripSize = strip.size();
    //int LB = LowerBound(totalItemWidth, stripWidth);

    //double avg = static_cast<double>(numItem) / static_cast<double>(stripSize);

    //cout << "Lower Bound: " << LB << endl;
    //cout << "# strips MFFD: " << stripSize << endl;
    //cout << "Avg # items per strip: " << avg << endl << endl;
    //cout << endl;

    //cout << stripSize << endl;
    //cout << avg << endl;

    //Optimality(opt, opt90, opt80, opt70, opt60, opt50, optLow, stripSize, LB);

} //End MFFD


// Packing each strip in turn, choosing smallest score width that meets vicinal sum constraint.
void PairSmallest(int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores, int numItem, int stripWidth,
                  double totalItemWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                  vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){
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
                if(adjMatrix[strip[j].back()][scoreIncrease[i]] == 1 && stripSum[j] + itemWidths[scoreIncrease[i]][partners[scoreIncrease[i]]] <= stripWidth){
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

    while(stripSum.back() == 0){
        stripSum.pop_back();
        strip.pop_back();
    }

    //int stripSize = strip.size();
    //int LB = LowerBound(totalItemWidth, stripWidth);

    //double avg = static_cast<double>(numItem) / static_cast<double>(stripSize);

    //cout << "Lower Bound: " << LB << endl;
    //cout << "# strips PairSmallest: " << stripSize << endl;
    //cout << "Avg # items per strip: " << avg << endl << endl;
    //cout << endl;

    //cout << stripSize << endl;
    //cout << avg << endl;

    //Optimality(opt, opt90, opt80, opt70, opt60, opt50, optLow, stripSize, LB);

} //End PairSmallest


// FFD including AHCA, instead of attempting to place item on end of strip, run AHCA to find feasible solution.
void MFFDPlus(int tau, int &opt, int &opt90, int &opt80, int &opt70, int &opt60, int &opt50, int &optLow, int numScores,
              int numItem, int maxItemWidth, int stripWidth, double totalItemWidth, vector<int> &allScores, vector<int> &partners,
              vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j, mini;
    int min = 0;
    int max = maxItemWidth;
    int feasible = 0;
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
                if(stripSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= stripWidth){
                    feasible = 0;
                    AHCAH(tau, i, j, feasible, allScores, partners, adjMatrix, itemWidths, itemDecrease, stripSum,
                          strip);
                    if(feasible == 1){
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

    /*
    for(j = 1; j < itemDecrease.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= stripWidth){
                    if(strip[i].size() == 2){ //If the strip only contains one item, don't run AHCA, just do checks instead
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
                        else if(adjMatrix[strip[i].front()][itemDecrease[j]] == 1){
                            strip[i].insert(strip[i].begin(), itemDecrease[j]);
                            strip[i].insert(strip[i].begin(), partners[itemDecrease[j]]);
                            stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        }
                        else if(adjMatrix[strip[i].front()][partners[itemDecrease[j]]] == 1){
                            strip[i].insert(strip[i].begin(), partners[itemDecrease[j]]);
                            strip[i].insert(strip[i].begin(), itemDecrease[j]);
                            stripSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        }
                    }
                    else {
                        feasible = 0;
                        AHCAH(tau, i, j, feasible, allScores, partners, adjMatrix, itemWidths, itemDecrease, stripSum,
                              strip);
                        if (feasible == 1) {
                            break;
                        }
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
    }*/

    while(stripSum.back() == 0){
        stripSum.pop_back();
        strip.pop_back();
    }

    //int stripSize = strip.size();
    //int LB = LowerBound(totalItemWidth, stripWidth);

    //double avg = static_cast<double>(numItem) / static_cast<double>(stripSize);

    //cout << "Lower Bound: " << LB << endl;
    //cout << "# strips MFFDPlus: " << stripSize << endl;
    //cout << "Avg # items per strip: " << avg << endl << endl;
    //cout << endl;

    //cout << stripSize << endl;
    //cout << avg << endl;

    //Optimality(opt, opt90, opt80, opt70, opt60, opt50, optLow, stripSize, LB);

} //End MFFDPlus


void FFD(int numScores, int numItem, int maxItemWidth, vector<int> &partners, vector<vector<int> > &itemWidths, vector<int> &itemOrder){

    int i, mini;
    int min = 0;
    int max = maxItemWidth;
    vector<int> checked(numScores, 0);

    while(itemOrder.size() < numItem) {
        for (i = 0; i < numScores; ++i) {
            if(checked[i] == 1){
                continue;
            }
            if (itemWidths[i][partners[i]] > min && itemWidths[i][partners[i]] <= max) {
                min = itemWidths[i][partners[i]];
                mini = i;
            }
        }
        itemOrder.push_back(mini);
        checked[mini] = 1;
        checked[partners[mini]] = 1;
        max = min;
        min = 0;
    }


}

void FFR(int numScores, int numItem, vector<int> &partners, vector<vector<int> > &itemWidths, vector<int> &itemOrder){

    int r;
    vector<int> checked(numScores, 0);

    while(itemOrder.size() < numItem){
        r = rand() % numScores;
        while(checked[r] == 1){
            r = rand() % numScores;
        }
        if(r < partners[r]){
            itemOrder.push_back(r);
            checked[r] = 1;
            checked[partners[r]] = 1;
        }
        else{
            itemOrder.push_back(partners[r]);
            checked[r] = 1;
            checked[partners[r]] = 1;
        }
    }



}

void FFShell(int numScores, int numItem, int maxItemWidth, int stripWidth, vector<int> &partners, vector<vector<int> > &adjMatrix,
             vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip, bool decrease){

    int i, j;
    vector<int> itemOrder;

    if(decrease){
        FFD(numScores, numItem, maxItemWidth, partners, itemWidths, itemOrder);
    }
    else{
        FFR(numScores, numItem, partners, itemWidths, itemOrder);
    }

    strip[0].push_back(itemOrder[0]);
    strip[0].push_back(partners[itemOrder[0]]);
    stripSum[0] += itemWidths[itemOrder[0]][partners[itemOrder[0]]];


    for(j = 1; j < itemOrder.size(); ++j){
        for(i = 0; i < strip.size(); ++i){
            if(!strip[i].empty()){
                if(stripSum[i] + itemWidths[itemOrder[j]][partners[itemOrder[j]]] <= stripWidth){
                    if(adjMatrix[strip[i].back()][itemOrder[j]] == 1){
                        strip[i].push_back(itemOrder[j]);
                        strip[i].push_back(partners[itemOrder[j]]);
                        stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                        break;
                    }
                    else if (adjMatrix[strip[i].back()][partners[itemOrder[j]]] == 1){
                        strip[i].push_back(partners[itemOrder[j]]);
                        strip[i].push_back(itemOrder[j]);
                        stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                        break;
                    }
                }
            }
            else if (strip[i].empty()){
                strip[i].push_back(itemOrder[j]);
                strip[i].push_back(partners[itemOrder[j]]);
                stripSum[i] += itemWidths[itemOrder[j]][partners[itemOrder[j]]];
                break;
            }
        }
    }

    while(stripSum.back() == 0){
        stripSum.pop_back();
        strip.pop_back();
    }


}

void PartialFFD(int numScores, int maxItemWidth, int stripWidth, vector<int> &partners, vector<vector<int> > &adjMatrix,
                vector<vector<int> > &itemWidths, vector<int> &partialItem, vector<int> &partialSum, vector<vector<int> > &partialSol){

    int i, j, mini;
    int min = 0;
    int max = maxItemWidth;
    vector<int> itemDecrease;
    vector<int> checked(numScores, 0);


    while(itemDecrease.size() < partialItem.size()/2) {
        for (i = 0; i < partialItem.size(); ++i) {
            if(checked[partialItem[i]] == 1){
                continue;
            }
            if (itemWidths[partialItem[i]][partners[partialItem[i]]] > min && itemWidths[partialItem[i]][partners[partialItem[i]]] <= max) {
                min = itemWidths[partialItem[i]][partners[partialItem[i]]];
                mini = partialItem[i];
            }
        }
        itemDecrease.push_back(mini);
        checked[mini] = 1;
        checked[partners[mini]] = 1;
        max = min;
        min = 0;
    }


    partialSol[0].push_back(itemDecrease[0]);
    partialSol[0].push_back(partners[itemDecrease[0]]);
    partialSum[0] += itemWidths[itemDecrease[0]][partners[itemDecrease[0]]];


    for(j = 1; j < itemDecrease.size(); ++j){
        for(i = 0; i < partialSol.size(); ++i){
            if(!partialSol[i].empty()){
                if(partialSum[i] + itemWidths[itemDecrease[j]][partners[itemDecrease[j]]] <= stripWidth){
                    if(adjMatrix[partialSol[i].back()][itemDecrease[j]] == 1){
                        partialSol[i].push_back(itemDecrease[j]);
                        partialSol[i].push_back(partners[itemDecrease[j]]);
                        partialSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        break;
                    }
                    else if (adjMatrix[partialSol[i].back()][partners[itemDecrease[j]]] == 1){
                        partialSol[i].push_back(partners[itemDecrease[j]]);
                        partialSol[i].push_back(itemDecrease[j]);
                        partialSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                        break;
                    }
                }
            }
            else if (partialSol[i].empty()){
                partialSol[i].push_back(itemDecrease[j]);
                partialSol[i].push_back(partners[itemDecrease[j]]);
                partialSum[i] += itemWidths[itemDecrease[j]][partners[itemDecrease[j]]];
                break;
            }
        }
    }

    while(partialSum.back() == 0){
        partialSum.pop_back();
        partialSol.pop_back();
    }

}

void CreateInitPop(int tau, int numPop, int numScores, int numItem, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners,
                   vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths,
                   vector<vector<int> > &populationSum, vector<vector<vector<int> > > &population){

    int i;
    vector<vector<int> > strip(numItem);
    vector<int> stripSum(numItem, 0);

    FFShell(numScores, numItem, maxItemWidth, stripWidth, partners, adjMatrix, itemWidths, stripSum, strip, true);

    Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);

    population.push_back(strip);
    populationSum.push_back(stripSum);

    for(i = 0; i < numPop - 1; ++i){
        strip.clear();
        strip.resize(numItem);
        stripSum.clear();
        stripSum.resize(numItem, 0);

        FFShell(numScores, numItem, maxItemWidth, stripWidth, partners, adjMatrix, itemWidths, stripSum, strip, false);

        Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);

        population.push_back(strip);
        populationSum.push_back(stripSum);

    }

}

void Mutation(int tau, int numScores, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
              vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j;
    vector<int> stripSumX;
    vector<vector<int> > stripX;
    vector<int> stripSumY;
    vector<vector<int> > stripY;
    vector<int> randOrder;


    for(i = 0; i < strip.size(); ++i){
        randOrder.push_back(i);
    }

    random_shuffle(randOrder.begin(), randOrder.end());

    int r = rand() % (strip.size() -1) + 1;

    for(i = 0; i < r; ++i){
        stripX.push_back(strip[randOrder[i]]);
        stripSumX.push_back(stripSum[randOrder[i]]);
    }

    for(i = r; i < randOrder.size(); ++i){
        stripY.push_back(strip[randOrder[i]]);
        stripSumY.push_back(stripSum[randOrder[i]]);
    }


    LocalSearch(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip, stripSumX, stripX, stripSumY, stripY);

}

void LocalSearch(int tau, int numScores, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                 vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip, vector<int> &stripSumX, vector<vector<int> > &stripX,
                 vector<int> &stripSumY, vector<vector<int> > &stripY){

    int a, b, c, d, i, j, k, l, pairSizeX, pairSizeY;
    int swapType, moveType, feasible;

    //region PairPair
    PairPair:
    /*SWAPPING A PAIR OF BOXES FROM EACH SET*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        if(stripX[i].size() >= 4){ //If there are at least 2 boxes on stripX[i] (note that each element represents a score, so 4 elements = 2 boxes)
            //Go through each pair of boxes on stripX[i]
            for(a = 0; a < stripX[i].size()-3; a+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                for(b = a+2; b < stripX[i].size()-1; b+=2){ //Starting from the first score on the second box until the first score on the last box
                    pairSizeX = itemWidths[stripX[i][a]][stripX[i][a+1]] + itemWidths[stripX[i][b]][stripX[i][b+1]]; //Sum box widths
                    //Check if there exists a pair of boxes on a strip in set stripY that have a combined width larger than pairSizeX
                    for(j = 0; j < stripY.size(); ++j){ //For each strip in the set stripY
                        if(stripY[j].size() >= 4){ //If there are at least 2 boxes on stripY[j]
                            //Go through each pair of boxes on stripY[j]
                            for(c = 0; c < stripY[j].size()-3; c+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                                for(d = c+2; d < stripY[j].size()-1; d+=2){ //Starting from the first score on the second box until the first score on the last box
                                    pairSizeY = itemWidths[stripY[j][c]][stripY[j][c+1]] + itemWidths[stripY[j][d]][stripY[j][d+1]]; //Sum box widths
                                    //Check if pairSizeX < pairSizeY and that boxes can fit onto strip
                                    if(pairSizeX < pairSizeY && stripSumX[i] - pairSizeX + pairSizeY <= stripWidth){
                                        swapType = 1;
                                        //cout << "i: " << i << " j: " << j << " a: " << a << " b: " << b << " c: " << c << " d: " << d << endl;
                                        if(stripX[i].size() == 4){ //If stripX[i] only contains 2 boxes
                                            if(stripY[j].size() == 4){ //If stripY[j] only contains 2 boxes
                                                //Do a straight Swap, no need for AHCA
                                                stripX[i].swap(stripY[j]);
                                                Swap(stripSumX[i], stripSumY[j]);
                                                feasible = 1;
                                            }
                                            else if (d == c+2){ //If stripY[j] contains more than 2 boxes & the two chosen boxes in stripY[j] are adjacent
                                                //Only perform AHCA on stripY[j]
                                                moveType = 11;
                                                InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                         itemWidths,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else{ //IIf stripY[j] contains more than 2 boxes & boxes c and d are not adjacent
                                                moveType = 0;
                                                InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                         itemWidths,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else if (stripY[j].size() == 4){ //If stripY[j] only contains 2 boxes but stripX[i] contains > 2 boxes
                                            if(b == a+2){ //If the two chosen boxes in stripX[i] are adjacent to one another
                                                //Only perform AHCA on stripX[i]
                                                moveType = 12;
                                                InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                         itemWidths,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                            else{ //If boxes a and b are not adjacent
                                                moveType = 0;
                                                InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                         itemWidths,
                                                         stripSumX, stripX, stripSumY, stripY);
                                            }
                                        }
                                        else{ //If stripX[i].size() > 4 && stripY[j[.size() > 4
                                            moveType = 0;
                                            InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                     itemWidths,
                                                     stripSumX, stripX, stripSumY, stripY);
                                        }
                                        if(feasible == 1){
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
    //endregion

    //region PairSin
    PairSin:
    /*SWAPPING A PAIR OF BOXES FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        if(stripX[i].size() >= 4){ //If there are at least 2 boxes on stripX[i]
            for(a = 0; a < stripX[i].size()-3; a+=2){ //Starting from the first score on the first box until the first score on the penultimate box
                for(b = a+2; b < stripX[i].size()-1; b+=2){ //Starting from the first score on the second box until the first score on the last box
                    pairSizeX = itemWidths[stripX[i][a]][stripX[i][a+1]] + itemWidths[stripX[i][b]][stripX[i][b+1]]; //Sum box widths
                    //Check if there exists a box on a strip in set stripY whose width is larger than pairSizeX
                    for(j = 0; j < stripY.size(); ++j){
                        //Go through each box on stripY[j]
                        for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box unil the first score on the last box
                            //Check if pairSizeX < width of box in stripY, and that box can fit onto strip
                            if(pairSizeX <= itemWidths[stripY[j][c]][stripY[j][c+1]] && stripSumX[i] - pairSizeX + itemWidths[stripY[j][c]][stripY[j][c+1]] <= stripWidth){
                                swapType = 2;
                                if(stripX[i].size() == 4){ //If stripX[i] only contains 2 boxes
                                    if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box
                                        //straight Swap
                                        stripX[i].swap(stripY[j]);
                                        Swap(stripSumX[i], stripSumY[j]);
                                        feasible = 1;
                                    }
                                    else{ //If stripY[j] contains more than 1 box
                                        //Only perform AHCA on stripY[j]
                                        moveType = 21;
                                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                 itemWidths,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box, but stripX[i] contains > 2 boxes
                                    if(b == a+2){ //If the two boxes chosen from stripX[i] are adjacent to one another
                                        //Only perform AHCA on stripX[i]
                                        moveType = 22;
                                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                 itemWidths,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                    else{ //If the two boxes chosen from stripX[i] are not adjacent to one another
                                        moveType = 0;
                                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores,
                                                 itemWidths,
                                                 stripSumX, stripX, stripSumY, stripY);
                                    }
                                }
                                else { //If stripX[i].size() > 4 && stripY[j].size() > 2
                                    moveType = 0;
                                    InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                             stripSumX, stripX, stripSumY, stripY);
                                }
                                if(feasible == 1){
                                    goto SinSin;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //endregion

    //region SinSin
    SinSin:
    /*SWAPPING ONE BOX FROM SET STRIPX WITH ONE BOX FROM SET STRIPY*/
    for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
        for(a = 0; a < stripX[i].size()-1; a+=2){ // Starting from the first score on the first box until the first score on the last box
            for(j = 0; j < stripY.size(); ++j){ // For each strip in the set stripY
                for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
                    //Check if boxwidth[a] < boxWidth[c] and that box can fit on strip
                    if(itemWidths[stripX[i][a]][stripX[i][a+1]] < itemWidths[stripY[j][c]][stripY[j][c+1]]
                       && stripSumX[i] - itemWidths[stripX[i][a]][stripX[i][a+1]] + itemWidths[stripY[j][c]][stripY[j][c+1]] <= stripWidth){
                        swapType = 3;
                        if(stripX[i].size() == 2){ //If stripX[i] only contains 1 box
                            if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box
                                //straight Swap
                                stripX[i].swap(stripY[j]);
                                Swap(stripSumX[i], stripSumY[j]);
                                feasible = 1;
                            }
                            else{ //If stripY[j] contains more than 1 box
                                //Only peform AHCA on stripY[j]
                                moveType = 31;
                                InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                         stripSumX, stripX, stripSumY, stripY);
                            }
                        }
                        else if(stripY[j].size() == 2){ //If stripY[j] only contains 1 box but stripX[i].size() > 2
                            //Only perform AHCA on stripX[i]
                            moveType = 32;
                            InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                     stripSumX, stripX, stripSumY, stripY);
                        }
                        else{ //stripX[i].size() > 2 && stripY[j].size() > 2
                            moveType = 0;
                            InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                     stripSumX, stripX, stripSumY, stripY);
                        }
                        if(feasible == 1){
                            goto MoveSin;
                        }
                    }
                }
            }
        }
    }
    //endregion

    //region MoveSin
    MoveSin:
    /*MOVING ONE BOX FROM SET STRIPY TO SET STRIPX*/
    for(j = 0; j < stripY.size(); ++j){ //For each strip in the set stripY
        for(c = 0; c < stripY[j].size()-1; c+=2){ //Starting from the first score on the first box until the first score on the last box
            for(i = 0; i < stripX.size(); ++i){ //For each strip in the set stripX
                if(stripSumX[i] + itemWidths[stripY[j][c]][stripY[j][c+1]] <= stripWidth){
                    swapType = 4;
                    if(stripY[j].size() == 2){ //If stripY[j] only contains one box
                        moveType = 41;
                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    else if(stripY[j].size() == 4){
                        moveType = 42;
                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    else {
                        moveType = 0;
                        InitAHCA(tau, swapType, moveType, feasible, i, a, b, j, c, d, allScores, itemWidths,
                                 stripSumX, stripX, stripSumY, stripY);
                    }
                    if(feasible == 2){
                        goto End;
                    }
                    else if(feasible == 1){
                        goto PairPair;
                    }
                }
            }
        }
    }
    //endregion



    End:
    if(feasible == 2){
        strip.clear();
        stripSum.clear();

        for(i = 0; i < stripX.size(); ++i){
            strip.push_back(stripX[i]);
            stripSum.push_back(stripSumX[i]);
        }
    }
    else{
        //Do FFD on stripY
        vector<int> partialItem;
        for(i = 0; i < stripY.size(); ++i){
            for(j = 0; j < stripY[i].size(); ++j){
                partialItem.push_back(stripY[i][j]);
            }
        }
        sort(partialItem.begin(), partialItem.end());

        stripY.clear();
        stripY.resize(partialItem.size()/2);
        stripSumY.clear();
        stripSumY.resize(partialItem.size()/2, 0);

        PartialFFD(numScores, maxItemWidth, stripWidth, partners, adjMatrix, itemWidths, partialItem, stripSumY, stripY);

        //join sets stripX and stripY together back into vector<vector<int> > strip

        strip.clear();
        stripSum.clear();

        for(i = 0; i < stripX.size(); ++i){
            strip.push_back(stripX[i]);
            stripSum.push_back(stripSumX[i]);
        }
        for(i = 0; i < stripY.size(); ++i){
            strip.push_back(stripY[i]);
            stripSum.push_back(stripSumY[i]);
        }
    }

}

void InitAHCA(int tau, int swapType, int moveType, int &feasible, int i1, int a1, int b1, int j1, int c1, int d1, vector<int> &allScores,
              vector<vector<int> > &itemWidths, vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY){

    int k;
    feasible = 0;
    vector<int> scoresX;
    vector<int> scoresY;
    vector<int> originalX;
    vector<int> originalY;
    vector<int> finalX;
    vector<int> finalY;

    //region swapType == 1
    /**PAIRPAIR**/
    if(swapType == 1){
        if(moveType == 11){ //X[i] = 2 items, Y[j] > 2 items, c and d adjacent, AHCA on Y[j] only.
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
            originalX.push_back(stripY[j1][d1]);
            originalX.push_back(stripY[j1][d1+1]);
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
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHCA on scoresY, originalY, finalY
            AHCAEA(tau, feasible, scoresY, originalY, finalY);
            if(feasible == 1){
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
        } //End moveType = 11

        else if(moveType == 12){ //Y[j] = 2 items, X[i] > 2 items, a and b adjacent, AHCA on X[i] only.
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                originalY.push_back(stripX[i1][a1]);
                originalY.push_back(stripX[i1][a1+1]);
                originalY.push_back(stripX[i1][b1]);
                originalY.push_back(stripX[i1][b1+1]);
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
        } //End moveType = 12

        else{ //moveType == 0, AHCA on both X and Y.
            /* X[i] = 2 items, Y[j] > 2 items, c and d not adjacent
             * Y[i] = 2 items, X[i] > 2 items, a and b not adjacent
             * X[i] > 2 items, Y[j] > 2 items */
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                feasible = 0;
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
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHCA on scoresY, originalY, finalY
                AHCAEA(tau, feasible, scoresY, originalY, finalY);
                if(feasible == 1){
                    stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]])
                                    + (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]]);
                    stripSumY[j1] = stripSumY[j1] - (itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]] + itemWidths[stripY[j1][d1]][stripY[j1][d1 + 1]])
                                    + (itemWidths[stripX[i1][a1]][stripX[i1][a1 + 1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1 + 1]]);
                    stripX[i1].swap(finalX);
                    stripY[j1].swap(finalY);
                }
            }
        } //End moveType = 0
    } //End swapType = 1
    //endregion

    //region swapType == 2
    /**PAIRSIN**/
    if(swapType == 2){
        if(moveType == 21){ //X[i] = 2 items, Y[j] > 1 item, AHCA on Y[j] only.
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
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
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHCA on scoresY, originalY, finalY
            AHCAEA(tau, feasible, scoresY, originalY, finalY);
            if(feasible == 1){
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]])
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]]
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]]);
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
        } //End moveType = 21
        else if(moveType == 22){ //Y[j] = 1 item, X[i] > 2 items, a and b adjacent, AHCA on X[i] only.
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                originalY.push_back(stripX[i1][a1]);
                originalY.push_back(stripX[i1][a1+1]);
                originalY.push_back(stripX[i1][b1]);
                originalY.push_back(stripX[i1][b1+1]);
                stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]])
                                + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]]
                                + (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]]);
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
        } // End moveType = 22

        else{ //moveType == 0, AHCA on both X and Y.
            /* Y[j] = 1 item, X[i] > 2 items, a and b not adjacent
             * X[i] > 2 items, Y[j] > 1 item */
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                feasible = 0;
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
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHCA on scoresY, originalY, finalY
                AHCAEA(tau, feasible, scoresY, originalY, finalY);
                if(feasible == 1){
                    stripSumX[i1] = stripSumX[i1] - (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]])
                                    + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                    stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]]
                                    + (itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripX[i1][b1]][stripX[i1][b1+1]]);
                    stripX[i1].swap(finalX);
                    stripY[j1].swap(finalY);
                }
            }
        } //End moveType = 0
    } //End swapType = 2
    //endregion

    //region swapType == 3
    /**SINSIN**/
    if(swapType == 3){
        if(moveType == 31){ //X[i] = 1 item, Y[j] > 1 item, AHCA on Y[j] only.
            originalX.push_back(stripY[j1][c1]);
            originalX.push_back(stripY[j1][c1+1]);
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
            scoresY.push_back(tau);
            scoresY.push_back(tau);
            //Run AHCA on scoresY, originalY, finalY
            AHCAEA(tau, feasible, scoresY, originalY, finalY);
            if(feasible == 1){
                stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]] + itemWidths[stripX[i1][a1]][stripX[i1][a1+1]];
                stripX[i1].swap(originalX);
                stripY[j1].swap(finalY);
            }
        } //End moveType = 31

        else if(moveType == 32){ //Y[j] = 1 item, X[i] > 1 item, AHCA on X[i] only.
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                originalY.push_back(stripX[i1][a1]);
                originalY.push_back(stripX[i1][a1+1]);
                stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]] + itemWidths[stripX[i1][a1]][stripX[i1][a1+1]];
                stripX[i1].swap(finalX);
                stripY[j1].swap(originalY);
            }
        } //End moveType = 32

        else{ //moveType == 0, X[i] > 1 item, Y[j] > 1 item, AHCA on both X and Y.
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
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                feasible = 0;
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
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHCA on scoresY, originalY, finalY
                AHCAEA(tau, feasible, scoresY, originalY, finalY);
                if(feasible == 1){
                    stripSumX[i1] = stripSumX[i1] - itemWidths[stripX[i1][a1]][stripX[i1][a1+1]] + itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                    stripSumY[j1] = stripSumY[j1] - itemWidths[stripY[j1][c1]][stripY[j1][c1+1]] + itemWidths[stripX[i1][a1]][stripX[i1][a1+1]];
                    stripX[i1].swap(finalX);
                    stripY[j1].swap(finalY);
                }
            }
        } //End moveType = 0
    } //End swapType = 3
    //endregion

    //region swapType == 4
    /**MOVESIN**/
    if(swapType == 4){
        if(moveType == 41){ //Y[j] = 1 item, AHCA on X[i] only (if feasible, row of Y[j] becomes empty and is deleted).
            for(k = 0; k < stripX[i1].size(); ++k){
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1+1]];
                stripX[i1].swap(finalX);
                stripY.erase(stripY.begin() + j1);
                stripSumY.erase(stripSumY.begin() + j1);
                if(stripY.empty()){
                    feasible = 2;
                }
            }
        } //End moveType = 41

        else if(moveType == 42){ //Y[j] = 2 items, AHCA on X[i] only (if feasible, Y[j] will only have 1 item left).
            for(k = 0; k < stripX[i1].size(); ++k){
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                if(c1 == 0){
                    stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] -= itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripY[j1].erase(stripY[j1].begin(), stripY[j1].begin() + 2);
                    stripX[i1].swap(finalX);
                }
                else if(c1 == 2){
                    stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] -= itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripY[j1].pop_back();
                    stripY[j1].pop_back();
                    stripX[i1].swap(finalX);
                }
                else{
                    cout << "[ERROR]: c1 in stripY[j1] is neither 0 nor 2, check that stripY[j1].size() == 4\n";
                    exit(1);
                }

            }
        } //End moveType = 42

        else{ //moveType = 0, Y[j] > 2 items, AHCA on both X and Y.
            for(k = 0; k < stripX[i1].size(); ++k){
                scoresX.push_back(allScores[stripX[i1][k]]);
                originalX.push_back(stripX[i1][k]);
            }
            scoresX.push_back(allScores[stripY[j1][c1]]);
            originalX.push_back(stripY[j1][c1]);
            scoresX.push_back(allScores[stripY[j1][c1 + 1]]);
            originalX.push_back(stripY[j1][c1 + 1]);
            scoresX.push_back(tau);
            scoresX.push_back(tau);
            //Run AHCA on scoresX, originalX, finalX
            AHCAEA(tau, feasible, scoresX, originalX, finalX);
            if(feasible == 1){
                feasible = 0;
                for(k = 0; k < stripY[j1].size(); ++k){
                    if(k == c1 || k == c1 + 1){
                        continue;
                    }
                    scoresY.push_back(allScores[stripY[j1][k]]);
                    originalY.push_back(stripY[j1][k]);
                }
                scoresY.push_back(tau);
                scoresY.push_back(tau);
                //Run AHCA on scoresY, originalY, finalY
                AHCAEA(tau, feasible, scoresY, originalY, finalY);
                if(feasible == 1){
                    stripSumX[i1] += itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripSumY[j1] -= itemWidths[stripY[j1][c1]][stripY[j1][c1 + 1]];
                    stripX[i1].swap(finalX);
                    stripY[j1].swap(finalY);
                }
            }
        } //End moveType = 0
    } //End swapType = 4
    //endregion

} //End void initAHCA

void EA(int tau, int recomb, int numScores, int maxItemWidth, int stripWidth, int &bestEnd, double &bestFitness, vector<int> &allScores, vector<int> &partners,
        vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths, vector<vector<int> > &populationSum,
        vector<vector<vector<int> > > &population){

    int i, j, k, l;
    vector<vector<int> > stripX;
    vector<vector<int> > stripY;
    vector<vector<int> > offspring;
    vector<int> stripSumX;
    vector<int> stripSumY;
    vector<int> offspringSum;

    /**Choose two solutions from population at random**/
    k = rand() % population.size();
    l = rand() % population.size();
    while (k == l){
        l = rand() % population.size();
    }


    for(i = 0; i < population[k].size(); ++i){
        stripX.push_back(population[k][i]);
        stripSumX.push_back(populationSum[k][i]);
    }

    for(i = 0; i < population[l].size(); ++i){
        stripY.push_back(population[l][i]);
        stripSumY.push_back(populationSum[l][i]);
    }


    double parent1Cost = Fitness(stripWidth, stripSumX, stripX);
    double parent2Cost = Fitness(stripWidth, stripSumY, stripY);

    //If GGA operator is chosen
    if(recomb == 1){
        GGA(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);
    }
        //If GPX' operator chosen
    else if(recomb == 2) {
        GPX(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);
    }

    double offspringCost = Fitness(stripWidth, offspringSum, offspring);

    if(parent1Cost < parent2Cost){ //Fitness of parent solution put in stripX worse than Fitness of parent solution put in stripY, replace population[k]
        offspring.swap(population[k]);
        offspringSum.swap(populationSum[k]);
        if(offspringCost > bestFitness){
            bestFitness = offspringCost;
            bestEnd = k;
        }
    }
    else if(parent1Cost > parent2Cost){ //Fitness of parent solution put in stripY worse than Fitness of parent solution put in stripX, replace population[l]
        offspring.swap(population[l]);
        offspringSum.swap(populationSum[l]);
        if(offspringCost > bestFitness){
            bestFitness = offspringCost;
            bestEnd = l;
        }
    }
    else if(parent1Cost == parent2Cost){
        double r = double (rand()) / double(RAND_MAX);
        if(r < 0.5){
            offspring.swap(population[k]);
            offspringSum.swap(populationSum[k]);
            if(offspringCost > bestFitness){
                bestFitness = offspringCost;
                bestEnd = k;
            }
        }
        else if (r >= 0.5){
            offspring.swap(population[l]);
            offspringSum.swap(populationSum[l]);
            if(offspringCost > bestFitness){
                bestFitness = offspringCost;
                bestEnd = l;
            }
        }
    }


}


void GGA(int tau, int numScores, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
         vector<vector<int> > &itemWidths, vector<int> &offspringSum, vector<vector<int> > &offspring,
         vector<int> &stripSumX, vector<vector<int> > &stripX, vector<int> &stripSumY, vector<vector<int> > &stripY){

    int i, j, k, l;
    vector<int> checked(numScores, 0);
    vector<int> absentItems;

    /**choose these randomly**/
    k = rand() % stripY.size();
    l = rand() % stripY.size();
    while(k >= l || (k == 0 && l == stripY.size() - 1)){
        k = rand() % stripY.size();
        l = rand() % stripY.size();
    }

    for(i = k; i <= l; ++i){
        for(j = 0; j < stripY[i].size(); ++j){
            checked[stripY[i][j]] = 1;
        }
    }

    for(i = 0; i < stripX.size(); ++i){
        for(j = 0; j < stripX[i].size(); ++j){
            if(checked[stripX[i][j]] == 1){
                stripX.erase(stripX.begin() + i);
                stripSumX.erase(stripSumX.begin() + i);
                --i;
                break;
            }
        }
    }


    for(i = k; i <= l; ++i){
        stripX.push_back(stripY[i]);
        stripSumX.push_back(stripSumY[i]);
    }

    for(i = 0; i < stripX.size(); ++i){
        for(j = 0; j < stripX[i].size(); ++j){
            checked[stripX[i][j]] = 1;
        }
    }

    for(i = 0; i < checked.size(); ++i){
        if(checked[i] == 0){
            absentItems.push_back(i);
        }
    }


    if(absentItems.empty()){
        stripX.swap(offspring);
        stripSumX.swap(offspringSum);
        Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring);
    }
    else if(absentItems.size() % 2 == 0){
        stripY.clear();
        stripY.resize(absentItems.size()/2);
        stripSumY.clear();
        stripSumY.resize(absentItems.size()/2, 0);

        PartialFFD(numScores, maxItemWidth, stripWidth, partners, adjMatrix, itemWidths, absentItems, stripSumY, stripY);

        LocalSearch(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);

        Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring);

    }
    else{
        cout << "[ERROR]: absentItems.size() is odd, not valid.\n";
        exit(1);
    }

}

void GPX(int tau, int numScores, int maxItemWidth, int stripWidth, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
         vector<vector<int> > &itemWidths, vector<int> &offspringSum, vector<vector<int> > &offspring, vector<int> &stripSumX, vector<vector<int> > &stripX,
         vector<int> &stripSumY, vector<vector<int> > &stripY){

    int i, j, k;
    double r = 2.0;
    vector<int> checked(numScores, 0);
    vector<int> absentItems;

    //If stripX and stripY each have a strip whose stripSums are equal and the largest of all other strips, choose between them randomly.
    if(*max_element(stripSumX.begin(), stripSumX.end()) == *max_element(stripSumY.begin(), stripSumY.end())){
        r = double(rand()) / double(RAND_MAX);
    }

    /**Fullest strip is in stripX**/
    if(*max_element(stripSumX.begin(), stripSumX.end()) > *max_element(stripSumY.begin(), stripSumY.end()) || r < 0.5){
        while(!stripX.empty() && !stripY.empty()) {
            k = distance(stripSumX.begin(), max_element(stripSumX.begin(), stripSumX.end()));

            //Mark the items in the chosen strip from stripX in the checked vector
            for (j = 0; j < stripX[k].size(); ++j) {
                checked[stripX[k][j]] = 1;
            }

            //Put the chosen strip from stripX into offspring
            offspring.push_back(stripX[k]);
            offspringSum.push_back(stripSumX[k]);
            stripX.erase(stripX.begin() + k);
            stripSumX.erase(stripSumX.begin() + k);

            //Go through strips in stripY, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for (i = 0; i < stripY.size(); ++i) {
                for (j = 0; j < stripY[i].size(); ++j) {
                    if (checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(stripY.empty()) {
                break;
            }

            //Now go to stripY and find the fullest strip
            k = distance(stripSumY.begin(), max_element(stripSumY.begin(), stripSumY.end()));
            for (j = 0; j < stripY[k].size(); ++j) {
                checked[stripY[k][j]] = 1;
            }
            offspring.push_back(stripY[k]);
            offspringSum.push_back(stripSumY[k]);
            stripY.erase(stripY.begin() + k);
            stripSumY.erase(stripSumY.begin() + k);

            for (i = 0; i < stripX.size(); ++i) {
                for (j = 0; j < stripX[i].size(); ++j) {
                    if (checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }

        } // End while


    } //End If fullest strip is in stripX

        /**Fullest strip is in stripY**/
    else if(*max_element(stripSumX.begin(), stripSumX.end()) < *max_element(stripSumY.begin(), stripSumY.end()) || r >= 0.5){
        while(!stripX.empty() && !stripY.empty()) {
            k = distance(stripSumY.begin(), max_element(stripSumY.begin(), stripSumY.end()));

            //Mark the items in the chosen strip from stripY in the checked vector
            for (j = 0; j < stripY[k].size(); ++j) {
                checked[stripY[k][j]] = 1;
            }

            //Put the chosen strip from stripX into offspring
            offspring.push_back(stripY[k]);
            offspringSum.push_back(stripSumY[k]);
            stripY.erase(stripY.begin() + k);
            stripSumY.erase(stripSumY.begin() + k);

            //Go through strips in stripX, delete the strips that contain any items that have been checked (i.e. that are in offspring)
            for (i = 0; i < stripX.size(); ++i) {
                for (j = 0; j < stripX[i].size(); ++j) {
                    if (checked[stripX[i][j]] == 1) {
                        stripX.erase(stripX.begin() + i);
                        stripSumX.erase(stripSumX.begin() + i);
                        --i;
                        break;
                    }
                }
            }

            if(stripX.empty()) {
                break;
            }

            //Now go to stripX and find the fullest strip
            k = distance(stripSumX.begin(), max_element(stripSumX.begin(), stripSumX.end()));
            for (j = 0; j < stripX[k].size(); ++j) {
                checked[stripX[k][j]] = 1;
            }

            offspring.push_back(stripX[k]);
            offspringSum.push_back(stripSumX[k]);
            stripX.erase(stripX.begin() + k);
            stripSumX.erase(stripSumX.begin() + k);

            for (i = 0; i < stripY.size(); ++i) {
                for (j = 0; j < stripY[i].size(); ++j) {
                    if (checked[stripY[i][j]] == 1) {
                        stripY.erase(stripY.begin() + i);
                        stripSumY.erase(stripSumY.begin() + i);
                        --i;
                        break;
                    }
                }
            }
        } // End while
    } //End If fullest strip is in stripY


    for(i = 0; i < checked.size(); ++i){
        if(checked[i] == 0){
            absentItems.push_back(i);
        }
    }

    if(absentItems.empty()){
        Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring);
    }
    else if(absentItems.size() % 2 == 0){
        stripX.swap(offspring);
        stripSumX.swap(offspringSum);
        offspring.clear();
        offspringSum.clear();
        stripY.clear();
        stripY.resize(absentItems.size()/2);
        stripSumY.clear();
        stripSumY.resize(absentItems.size()/2, 0);

        PartialFFD(numScores, maxItemWidth, stripWidth, partners, adjMatrix, itemWidths, absentItems, stripSumY, stripY);

        LocalSearch(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring, stripSumX, stripX, stripSumY, stripY);

        Mutation(tau, numScores, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, offspringSum, offspring);

    }
    else{
        cout << "[ERROR]: absentItems.size() is odd, not valid.\n";
        exit(1);
    }


}

void AHCAH(int tau, int i1, int j1, int &feasible, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
           vector<vector<int> > &itemWidths, vector<int> &itemDecrease, vector<int> &stripSum, vector<vector<int> > &strip){

    int k;
    vector<int> scores;
    vector<int> original;
    feasible = 0;

    int dom = 150; /****/ //dominating score widths

    //Creating scores vector
    for(k = 0; k < strip[i1].size(); ++k){
        scores.push_back(allScores[strip[i1][k]]);
        original.push_back(strip[i1][k]);
    }
    scores.push_back(allScores[itemDecrease[j1]]);
    original.push_back(itemDecrease[j1]);
    scores.push_back(allScores[partners[itemDecrease[j1]]]);
    original.push_back(partners[itemDecrease[j1]]);
    scores.push_back(dom); /****/
    scores.push_back(dom); /****/

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
        InitInstance(tau, nScores, adjMat, scores, order, partnersX);

        MCM(nScores, matchSize, adjMat, partnersX, matchList, cycleVertex);
        if (matchSize < nItem) {
            feasible = 0;
            break;
        }

        MPS(nScores, nCycles, partnersX, matchList, mpStructure);
        if (mpStructure[0].size() == nScores) {
            for (j = 2; j < mpStructure[0].size(); ++j) {
                altHam.push_back(mpStructure[0][j]);
            }

            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            feasible = 1;
            break;
        }

        BR(qstar, matchSize, adjMat, matchList, cycleVertex, edge, mpStructure, C, S);
        if (qstar == -1) {
            feasible = 0;
            break;
        }

        CP(nScores, nComp, feasible, qstar, nCycles, partnersX, matchList, cycleVertex, edge, adjMat, C, S, altHam);
        if (feasible == 1) {
            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            break;
        }

    } while (once == 0);

    if(feasible == 1){
        stripSum[i1] += itemWidths[itemDecrease[j1]][partners[itemDecrease[j1]]];
        strip[i1].swap(final);
    }

} //End AHCA


void AHCAEA(int tau, int &feasible, vector<int> &scores, vector<int> &original, vector<int> &final){

    feasible = 0;
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
        InitInstance(tau, nScores, adjMat, scores, order, partnersX);

        MCM(nScores, matchSize, adjMat, partnersX, matchList, cycleVertex);
        if (matchSize < nItem) {
            feasible = 0;
            break;
        }

        MPS(nScores, nCycles, partnersX, matchList, mpStructure);
        if (mpStructure[0].size() == nScores) {
            for (j = 2; j < mpStructure[0].size(); ++j) {
                altHam.push_back(mpStructure[0][j]);
            }

            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            feasible = 0;
            break;
        }

        BR(qstar, matchSize, adjMat, matchList, cycleVertex, edge, mpStructure, C, S);
        if (qstar == -1) {
            feasible = 0;
            break;
        }

        CP(nScores, nComp, feasible, qstar, nCycles, partnersX, matchList, cycleVertex, edge, adjMat, C, S, altHam);
        if (feasible == 1) {
            for (i = 0; i < altHam.size(); ++i) {
                final.push_back(original[order[altHam[i]]]);
            }
            break;
        }

    } while (once == 0);


}


void InitInstance(int tau, int nScores, vector<vector<int> > &adjMat, vector<int> &scores, vector<int> &order, vector<int> &partnersX){

    int i, j;
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
        for (j = 0; j < nScores; ++j) {
            if (adjMat[i][j] == 2) {
                partnersX[i] = j;
            }
        }

    }



}


void MCM(int nScores, int &matchSize, vector<vector<int> > &adjMat, vector<int> &partnersX, vector<int> &matchList,
         vector<int> &cycleVertex){

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

}


void MPS(int nScores, int &nCycles, vector<int> &partnersX, vector<int> &matchList, vector<vector<int> > &mpStructure){

    int i, current;
    int smallest = nScores - 2;
    vector<int> temp;
    vector<int> checked(nScores, 0);

    do {
        current = smallest;
        temp.clear();
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

}


void BR(int &qstar, int matchSize, vector<vector<int> > &adjMat, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &edge,
        vector<vector<int> > &mpStructure, vector<vector<int> > &C, vector<vector<int> > &S){

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

}


void CP(int nScores, int nComp, int &feasible, int qstar, int nCycles, vector<int> &partnersX, vector<int> &matchList,
        vector<int> &cycleVertex, vector<int> &edge, vector<vector<int> > &adjMat, vector<vector<int> > &C, vector<vector<int> > &S, vector<int> &altHam){

    int a, i, j, k, l, q, u, v, SSum, SqIntS;
    int v1 = 0;
    int v2 = 0;
    int vacant = 999;
    int full = vacant;
    int current = nScores - 2;
    int maxRowSize = 0;
    int maxRow;
    int type;
    vector<int> temp;
    vector<int> SSet2;
    vector<int> SSet3;
    vector<int> edgeCopy;
    vector<int> connectML;
    vector<int> QSet(nComp, 0);
    vector<int> inCycle(nScores, 0);
    vector<int> connectCycle(nComp, vacant); //was patchCycleX
    vector<vector<int> > Cconnect;
    vector<vector<int> > Cq;


    for (i = 0; i < C.size(); ++i) {
        if (C[i].size() == nCycles) {
            full = i;
            break;
        }
    }

    if (full != vacant) {

        connectML = matchList;

        for (v = 0; v < C[full].size() - 1; ++v) {
            connectML[C[full][v]] = matchList[C[full][v + 1]];
            connectML[matchList[C[full][v + 1]]] = C[full][v];
        }
        connectML[C[full][C[full].size() - 1]] = matchList[C[full][0]];
        connectML[matchList[C[full][0]]] = C[full][C[full].size() - 1];

        do {
            altHam.push_back(current);
            altHam.push_back(partnersX[current]);
            current = connectML[partnersX[current]];
        } while (altHam.size() < nScores);

        altHam.erase(altHam.begin(), altHam.begin() + 2);
        feasible = 1;

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
            edgeCopy = edge;
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
            bool added = false;
            int lastRow = Cconnect.size();

            while (!added) {
                k = 0;
                do {
                    while (k < nEdgesC - 2 && (adjMat[edgeCopy[k]][matchList[edgeCopy[k + 1]]] != 1 || cycleVertex[edgeCopy[k]] == cycleVertex[edgeCopy[k + 1]])) {
                        ++k;
                    }

                    if (adjMat[edgeCopy[k]][matchList[edgeCopy[k + 1]]] == 1 && cycleVertex[edgeCopy[k]] != cycleVertex[edgeCopy[k + 1]]
                        && ((SSet2[cycleVertex[edgeCopy[k]]] == 0 && SSet2[cycleVertex[edgeCopy[k + 1]]] == 1)
                        || (SSet2[cycleVertex[edgeCopy[k]]] == 1 && SSet2[cycleVertex[edgeCopy[k + 1]]] == 0))) {
                        temp.push_back(edgeCopy[k]);
                        temp.push_back(edgeCopy[k + 1]);
                        SSet2[cycleVertex[edgeCopy[k]]] = 1;
                        SSet2[cycleVertex[edgeCopy[k + 1]]] = 1;
                        ++SSum;
                        if (SSum < nCycles) {
                            ++k;
                            while (k < nEdgesC - 1 && SSet2[cycleVertex[edgeCopy[k + 1]]] == 0 && adjMat[edgeCopy[k]][matchList[edgeCopy[k + 1]]] == 1) {
                                ++k;
                                temp.push_back(edgeCopy[k]);
                                SSet2[cycleVertex[edgeCopy[k]]] = 1;
                                ++SSum;
                            }
                        }
                        Cconnect.push_back(temp);
                        added = true;
                        temp.clear();
                    }
                    ++k;
                } while (k < nEdgesC - 1 && SSum < nCycles);

                if (SSum == nCycles) {
                    type = 2;
                    break;
                }
                else if(added == true){ // &&SSum < nCycles
                    for(i = lastRow; i < Cconnect.size(); ++i){
                        for(j = 0; j < Cconnect[i].size(); ++j){
                            for(l = 0; l < edgeCopy.size(); ++l){
                                if(edgeCopy[l] == Cconnect[i][j]){
                                    edgeCopy.erase(edgeCopy.begin()+l);
                                    break;
                                }
                            }
                        }
                    }
                    lastRow = Cconnect.size();
                    added = false;
                    nEdgesC = edgeCopy.size();
                }
                else if(added == false){
                    break;
                }

            }//end while
        }// End find overlaps
        //endregion

        //region TYPE 3
        //TYPE 3: Original
        if(type == 0){
            Cconnect.clear();
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
            }
            else if(type == 3){
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

            connectML = matchList;

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
                altHam.push_back(partnersX[current]);
                current = connectML[partnersX[current]];
            } while (altHam.size() < nScores);

            altHam.erase(altHam.begin(), altHam.begin() + 2);
            feasible = 1;
        }

        else {
            feasible = 0;
        }
    }

}


