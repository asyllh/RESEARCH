/*--------------/
ALH
main.cpp
Evolutionary Algorithm with Local Search
05/12/2017
06/03/2018
/--------------*/
#include <iostream>
#include <fstream>
#include <ctime>
#include <cstring>
#include <iomanip>
using namespace std;

#include "base.h"
#include "packing.h"

void programInfo(){

    cout << "Evolutionary Algorithm for the SCSPP:\n-------------\n"
         << "PARAMETERS:\n"
         << "       -i <int>    [Number of instances. Default = 1000.]\n"
         << "       -t <int>    [Constraint value. Default = 70.]\n"
         << "       -n <int>    [Number of items. Default = 500.]\n"
         << "       -a <int>    [Minimum score width. Default = 1.]\n"
         << "       -b <int>    [Maximum score width. Default = 70.]\n"
         << "       -w <int>    [Minimum item width. Default = 150.]\n"
         << "       -W <int>    [Maximum item width. Default = 1000.]\n"
         << "       -l <int>    [Length of strips. Default = 5000.]\n"
         << "       -s <int>    [Random seed. Default = 1.]\n"
         << "---------------\n\n";
}

void argumentCheck(int numInstances, int tau, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth,
                   int maxStripWidth, int randomSeed){
    cout << "Evolutionary Algorithm for the SCSPP\n------------------------------\n";
    if(tau == 0){
        cout << "[ERROR]: Constraint value cannot be zero.\n";
        exit(1);
    }
    if(maxStripWidth == 0){
        cout << "[ERROR]: Strip cannot have length zero.\n";
        exit(1);
    }
    if(2*minWidth >= tau){
        cout << "[ERROR]: Constraint value is less than or equal to twice the minimum score width, vicinal sum constraint always valid.\n";
        cout << "         Problem instance is therefore classical strip-packing problem without score constraint (i.e. tau = 0).\n";
        exit(1);
    }
    if(2*maxWidth < tau){
        cout << "[ERROR]: Constraint value is greater than double maximum score width, vicinal sum constraint never valid.\n";
        exit(1);
    }
    if(2*maxWidth >= minItemWidth){
        cout << "[ERROR]: Minimum item width is less than double maximum score width, scores may overlap.\n";
        exit(1);
    }
    if(minWidth > maxWidth){
        cout << "[ERROR]: Minimum score width is greater than maximum score width.\n";
        exit(1);
    }
    if(maxItemWidth > maxStripWidth){
        cout << "[ERROR]: Maximum item width is larger than length of strip.\n";
        exit(1);
    }

    cout << std::left << setw(20) << "Number of instances:" << std::right << setw(10) << numInstances << endl
         << std::left << setw(20) << "Constraint value:" << std::right << setw(10) << tau << endl
         << std::left << setw(20) << "Number of items:" << std::right << setw(10) << numItem << endl
         << std::left << setw(20) << "Minimum score width:" << std::right << setw(10) << minWidth << endl
         << std::left << setw(20) << "Maximum score width:" << std::right << setw(10) << maxWidth << endl
         << std::left << setw(20) << "Minimum item width:" << std::right << setw(10) << minItemWidth << endl
         << std::left << setw(20) << "Maxmimum item width:" << std::right << setw(10) << maxItemWidth << endl
         << std::left << setw(20) << "Length of strips:" << std::right << setw(10) << maxStripWidth << endl
         << std::left << setw(20) << "Random seed:" << std::right << setw(10) << randomSeed << endl;
    cout << "------------------------------\n\n";
}

int main(int argc, char **argv){
    if(argc <= 1){
       programInfo();
        exit(1);
    }

    int x;
    int numInstances = 1000;
    int tau = 70;
    int numItem = 500;
    int minWidth = 1;
    int maxWidth = 70;
    int minItemWidth = 150;
    int maxItemWidth = 1000;
    int maxStripWidth = 5000;
    int randomSeed = 1;

    //region User Arguments
    for(x = 1; x < argc; ++x){
        if(strcmp("-i", argv[x]) == 0){
            numInstances = atoi(argv[++x]);
        }
        else if(strcmp("-t", argv[x]) == 0){
            tau = atoi(argv[++x]);
        }
        else if(strcmp("-n", argv[x]) == 0){
            numItem = atoi(argv[++x]);
        }
        else if(strcmp("-a", argv[x]) == 0){
            minWidth = atoi(argv[++x]);
        }
        else if(strcmp("-b", argv[x]) == 0){
            maxWidth = atoi(argv[++x]);
        }
        else if(strcmp("-w", argv[x]) == 0){
            minItemWidth = atoi(argv[++x]);
        }
        else if(strcmp("-W", argv[x]) == 0){
            maxItemWidth = atoi(argv[++x]);
        }
        else if(strcmp("-l", argv[x]) == 0){
            maxStripWidth = atoi(argv[++x]);
        }
        else if(strcmp("-s", argv[x]) == 0){
            randomSeed = atoi(argv[++x]);
        }
    }
    //endregion

    argumentCheck(numInstances, tau, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, maxStripWidth, randomSeed);

    int i;
    int instance;
    int numScores = numItem * 2;
    double totalItemWidth;
    vector<int> allScores;
    vector<int> partners(numScores, 0);
    vector<vector<int> > itemWidths(numScores, vector<int>(numScores, 0));
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0));
    vector<vector<int> > allItems(numScores, vector<int>(numScores, 0));
    vector<vector<vector<int> > > population;
    vector<vector<int> > populationSum;
    double parent1cost;
    double parent2cost;
    srand(randomSeed);


    time_t startTime, endTime; //start clock
    startTime = clock();

    createInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths, allItems);

    createInitialPopulation(tau, numScores, numItem, maxItemWidth, maxStripWidth, allScores, partners, adjMatrix, itemWidths, populationSum, population);

    EA(tau, numScores, maxItemWidth, maxStripWidth, parent1cost, parent2cost, allScores, partners, adjMatrix, itemWidths, populationSum, population);


    endTime = clock();
    double totalTimeMS = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    double totalTimeS = totalTimeMS / 1000;
    cout << "\nCPU Time = " << totalTimeMS << " milliseconds (" << totalTimeS << " seconds.)\nEND.\n";


}//END INT MAIN


