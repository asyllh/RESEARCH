/*--------------/
ALH
base.cpp
12/10/17
/--------------*/
#include <algorithm>
#include <iomanip>
#include "base.h"
using namespace std;

void resetVectors(int numScores, int numItem, vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix,
                  vector<vector<int> > &itemWidths, vector<int> &stripSum, vector<vector<int> > &strip){
    int i, j;

    for(i = 0; i < numScores; ++i){
        partners[i] = 0;
        for(j = 0; j < numScores; ++j){
            adjMatrix[i][j] = 0;
            itemWidths[i][j] = 0;
        }
    }

    while(!strip.empty()){
        strip.pop_back();
    }

    while(!stripSum.empty()){
        stripSum.pop_back();
    }

    vector<int> temp;
    for(i = 0; i < numItem; ++i){
        strip.push_back(temp);
        stripSum.push_back(0);
    }

    allScores.clear();

} //End resetVectors


void output(int opt, int opt90, int opt80, int opt70, int opt60, int opt50, int optLow, int numInstances){

    cout << setprecision(2) << fixed;
    cout << "Number of Instances: " << numInstances << endl;
    cout << "#        opt of LB: " << opt << "\t\tPercentage: "<< (static_cast<double>(opt)/numInstances) * 100 << "\t\t# 100%: " << opt << endl;
    cout << "# 90% - <opt of LB: " << opt90 << "\t\tPercentage: "<< (static_cast<double>(opt90)/numInstances) * 100 << "\t\t# >= 90%: " << opt + opt90 << endl;
    cout << "# 80% - <90% of LB: " << opt80 << "\t\tPercentage: "<< (static_cast<double>(opt80)/numInstances) * 100 << "\t\t# >= 80%: " << opt + opt90 + opt80 << endl;
    cout << "# 70% - <80% of LB: " << opt70 << "\t\tPercentage: "<< (static_cast<double>(opt70)/numInstances) * 100 << "\t\t# >= 70%: " << opt + opt90 + opt80 + opt70 << endl;
    cout << "# 60% - <70% of LB: " << opt60 << "\t\tPercentage: "<< (static_cast<double>(opt60)/numInstances) * 100 << "\t\t# >= 60%: " << opt + opt90 + opt80 + opt70 + opt60 << endl;
    cout << "# 50% - <60% of LB: " << opt50 << "\t\tPercentage: "<< (static_cast<double>(opt50)/numInstances) * 100 << "\t\t# >= 50%: " << opt + opt90 + opt80 + opt70 + opt60 + opt50 << endl;
    cout << "#       <50% of LB: " << optLow << "\t\tPercentage: "<< (static_cast<double>(optLow)/numInstances) * 100 << endl;

} //End output


void createInstance(int instance, int tau, int numScores, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, double &totalItemWidth,
                    vector<int> &allScores, vector<int> &partners, vector<vector<int> > &adjMatrix, vector<vector<int> > &itemWidths){

    int i, j;
    int count = 1;
    vector<int> randOrder;
    vector<int> checkItem(numScores, 0);
    totalItemWidth = 0.0;

    for (i = 0; i < numScores; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }

    sort(allScores.begin(), allScores.end());

    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= tau) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    for (i = 0; i < numScores; ++i) {
        randOrder.push_back(i);
    }

    random_shuffle(randOrder.begin(), randOrder.end());

    for (i = 0; i < numItem; ++i) {
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                partners[i] = j;
                break;
            }
        }
    }

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && itemWidths[i][j] == 0){
                itemWidths[i][j] = rand() % (maxItemWidth - minItemWidth + 1) + minItemWidth;
                itemWidths[j][i] = itemWidths[i][j];
                break;
            }

        }
    }

    //cout << "Box#" << setw(10) << "Mates" << setw(10) << "Scores" << setw(10) << "Width\n";
    /*for(i = 0; i < numScores; ++i){
        if(checkItem[i] == 1){
            continue;
        }
        //cout << count << setw(10)  << i << "-" << mates[i] << setw(10) << allScores[i] << "-" << allScores[mates[i]] << setw(8) << itemWidths[i][mates[i]] << endl;
        totalItemWidth += itemWidths[i][partners[i]];
        checkItem[i] = 1;
        checkItem[partners[i]] = 1;
        ++count;

    }*/



    if(instance == 0){
        cout << "INSTANCE = " << instance << endl;
        cout << "Box#" << setw(10) << "Mates" << setw(10) << "Scores" << setw(10) << "Width\n";
        for(i = 0; i < numScores; ++i){
            if(checkItem[i] == 1){
                continue;
            }
            cout << count << setw(10)  << i << "-" << partners[i] << setw(10) << allScores[i] << "-" << allScores[partners[i]] << setw(8) << itemWidths[i][partners[i]] << endl;
            totalItemWidth += itemWidths[i][partners[i]];
            checkItem[i] = 1;
            checkItem[partners[i]] = 1;
            ++count;

        }

    }
    else{
        //cout << "Box#" << setw(10) << "Mates" << setw(10) << "Scores" << setw(10) << "Width\n";
        for(i = 0; i < numScores; ++i){
            if(checkItem[i] == 1){
                continue;
            }
            //cout << count << setw(10)  << i << "-" << mates[i] << setw(10) << allScores[i] << "-" << allScores[mates[i]] << setw(8) << itemWidths[i][mates[i]] << endl;
            totalItemWidth += itemWidths[i][partners[i]];
            checkItem[i] = 1;
            checkItem[partners[i]] = 1;
            ++count;

        }

    }



} //End createInstance
