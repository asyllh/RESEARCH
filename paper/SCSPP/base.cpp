/*--------------/
ALH
base.cpp
12/10/17
/--------------*/
#include <algorithm>
#include <iomanip>
#include "base.h"
using namespace std;

void resetVectors(int numScores, int numBox, vector<int> &allScores, vector<int> &mates,
                  vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths, vector<int> &stripSum, vector<vector<int> > &strip){

    int i, j;

    for(i = 0; i < numScores; ++i){
        mates[i] = 0;
        for(j = 0; j < numScores; ++j){
            adjMatrix[i][j] = 0;
            boxWidths[i][j] = 0;
        }
    }

    while(!strip.empty()){
        strip.pop_back();
    }

    while(!stripSum.empty()){
        stripSum.pop_back();
    }

    vector<int> temp;
    for(i = 0; i < numBox; ++i){
        strip.push_back(temp);
        stripSum.push_back(0);
    }

    allScores.clear();



}

void output(int opt, int opt90, int opt80, int opt70, int opt60, int opt50, int optLow, int numInstances){

    cout << setprecision(2) << fixed;
    cout << "Number of Instances: " << numInstances << endl;
    cout << "#        opt of LB: " << opt << "\t\tPercentage: "<< static_cast<double>(opt)/numInstances << "\t\t# 100%: " << opt << endl;
    cout << "# 90% - <opt of LB: " << opt90 << "\t\tPercentage: "<< static_cast<double>(opt90)/numInstances << "\t\t# >= 90%: " << opt + opt90 << endl;
    cout << "# 80% - <90% of LB: " << opt80 << "\t\tPercentage: "<< static_cast<double>(opt80)/numInstances << "\t\t# >= 80%: " << opt + opt90 + opt80 << endl;
    cout << "# 70% - <80% of LB: " << opt70 << "\t\tPercentage: "<< static_cast<double>(opt70)/numInstances << "\t\t# >= 70%: " << opt + opt90 + opt80 + opt70 << endl;
    cout << "# 60% - <70% of LB: " << opt60 << "\t\tPercentage: "<< static_cast<double>(opt60)/numInstances << "\t\t# >= 60%: " << opt + opt90 + opt80 + opt70 + opt60 << endl;
    cout << "# 50% - <60% of LB: " << opt50 << "\t\tPercentage: "<< static_cast<double>(opt50)/numInstances << "\t\t# >= 50%: " << opt + opt90 + opt80 + opt70 + opt60 + opt50 << endl;
    cout << "#       <50% of LB: " << optLow << "\t\tPercentage: "<< static_cast<double>(optLow)/numInstances << endl;

}

void createInstance(int numScores, int numBox, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, double &totalBoxWidth,
                    vector<int> &allScores, vector<int> &mates, vector<vector<int> > &adjMatrix, vector<vector<int> > &boxWidths){

    int i, j;
    int threshold = 70;
    int count = 1;
    vector<int> randOrder;
    vector<int> checkBox(numScores, 0);
    totalBoxWidth = 0.0;

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for (i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }


    //Sort all of the scores in the allScores vector in ascending order
    sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    /*cout << "All scores:\n";
    for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl;*/
    //Filling in adjacency matrix - if sum of two scores >= threshold (70), then insert 1 into the matrix, else leave as 0
    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= threshold) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    //Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for (i = 0; i < numScores; ++i) {
        randOrder.push_back(i);
    }

    //Randomly shuffle all values in randOrder vector
    random_shuffle(randOrder.begin(), randOrder.end());

    //Assign mates to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for (i = 0; i < numBox; ++i) {
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    /*cout << "AdjMatrix:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;*/



    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                mates[i] = j;
                break;
            }
        }
    }
    /*cout << "Mates Vector:\n";
    for(i = 0; i < mates.size(); ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;*/

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && boxWidths[i][j] == 0){
                boxWidths[i][j] = rand() % (maxBoxWidth - minBoxWidth + 1) + minBoxWidth;
                boxWidths[j][i] = boxWidths[i][j];
                break;
            }

        }
    }


    //cout << "Box#" << setw(10) << "Mates" << setw(10) << "Scores" << setw(10) << "Width\n";
    for(i = 0; i < numScores; ++i){
        if(checkBox[i] == 1){
            continue;
        }
        //cout << count << setw(10)  << i << "-" << mates[i] << setw(10) << allScores[i] << "-" << allScores[mates[i]] << setw(8) << boxWidths[i][mates[i]] << endl;
        totalBoxWidth += boxWidths[i][mates[i]];
        checkBox[i] = 1;
        checkBox[mates[i]] = 1;
        ++count;

    }

   //cout << "Total Box Widths: " << totalBoxWidth << endl << endl;



    /*cout << "allBoxes:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allBoxes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}

void createInstanceUser(int threshold, int numScores, double &totalBoxWidth, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes){

    int i, j, m1, m2, k;
    int count = 1;
    vector<int> checkBox(numScores, 0);

    for(i = 0; i < userInput.size(); ++i){
        for(j = 0; j < userInput[i].size() - 1; ++j){
            allScores.push_back(userInput[i][j]);
        }
    }

    sort(allScores.begin(), allScores.end());

    cout << "AllScores User:\n";
    for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl << endl;

    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= threshold) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    vector<int>::iterator it1;
    vector<int>::iterator it2;

    for(i = 0; i < userInput.size(); ++i) {
        it1 = find(allScores.begin(), allScores.end(), userInput[i][0]);
        m1 = it1 - allScores.begin();
        it2 = find(allScores.begin(), allScores.end(), userInput[i][1]);
        m2 = it2 - allScores.begin();
        //cout << "position of pairs " << userMates[i][0] << " and " << userMates[i][1] << ": " << m1 << "-" << m2 << endl;
        adjMatrix[m1][m2] = 2;
        adjMatrix[m2][m1] = 2;
        boxWidths[m1][m2] = userInput[i][2];
        boxWidths[m2][m1] = userInput[i][2];
    }

    /*cout << "AdjMatrix:\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

    /*cout << "BoxWidths:\n";
    for(i = 0; i < numScores; ++i){
        if(checkBox[i] == 1){
            continue;
        }
        cout << i << "-" << mates[i] << ": " << boxWidths[i][mates[i]] << endl;
        checkBox[i] = 1;
        checkBox[mates[i]] = 1;
    }*/

    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                mates[i] = j;
                break;
            }
        }
    }

    cout << "Mates:\n";
    for(i = 0; i < mates.size(); ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;

    k = 1;
    for(i = 0; i < numScores; ++i){
        for(j = i+1; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                allBoxes[i][j] = k;
                allBoxes[j][i] = k * -1;
                ++k;
                break;
            }
        }
    }

    cout << "Box#" << setw(10) << "Mates" << setw(10) << "Width\n";
    for(i = 0; i < numScores; ++i){
        if(checkBox[i] == 1){
            continue;
        }
        cout << count << setw(10) << i << "-" << mates[i] << setw(8) << boxWidths[i][mates[i]] << endl;
        totalBoxWidth += boxWidths[i][mates[i]];
        checkBox[i] = 1;
        checkBox[mates[i]] = 1;
        ++count;

    }

    cout << "Total Box Widths: " << totalBoxWidth << endl << endl;

    /*cout << "allBoxes:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allBoxes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}