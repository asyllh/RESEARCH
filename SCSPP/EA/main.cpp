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
#include <cmath>
using namespace std;

#include "base.h"
#include "packing.h"

void ProgramInfo(){

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
         << "       -p <int>    [Number of solutions in population.]\n"
         << "       -r <int>    [Recombination operator. 1: GGA. 2: GPX'.]\n"
         << "       -s <int>    [Random seed. Default = 1.]\n"
         << "---------------\n\n";
}

void ArgumentCheck(int numInstances, int tau, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth,
                   int stripLength, int numPop, int xOver, int randomSeed){

    bool error = false;

    cout << "Evolutionary Algorithm for the SCSPP\n------------------------------\n";
    if(tau == 0){
        cout << "[ERROR]: Constraint value cannot be zero.\n";
        error = true;
    }
    if(stripLength == 0){
        cout << "[ERROR]: Strip cannot have length zero.\n";
        error = true;
    }
    if(2*minWidth >= tau){
        cout << "[ERROR]: Constraint value is less than or equal to twice the minimum score width, vicinal sum constraint always valid.\n";
        cout << "         Problem instance is therefore classical strip-packing problem without score constraint (i.e. tau = 0).\n";
        error = true;
    }
    if(2*maxWidth < tau){
        cout << "[ERROR]: Constraint value is greater than double maximum score width, vicinal sum constraint never valid.\n";
        error = true;
    }
    if(2*maxWidth >= minItemWidth){
        cout << "[ERROR]: Minimum item width is less than double maximum score width, scores may overlap.\n";
        error = true;
    }
    if(minWidth > maxWidth){
        cout << "[ERROR]: Minimum score width is greater than maximum score width.\n";
        error = true;
    }
    if(maxItemWidth > stripLength){
        cout << "[ERROR]: Maximum item width is larger than length of strip.\n";
        error = true;
    }
    if(numPop < 5){
        cout << "[ERROR]: Insufficient number of solutions in population.\n";
        error = true;
    }
    if(xOver != 1 && xOver != 2){
        cout << "[ERROR]: Invalid choice of recombination operator. Please choose either 1: GGA, or 2: GPX'.\n";
        error = true;
    }

    if(error){
        cout << "[EXIT PROGRAM.]\n";
        exit(1);
    }

    cout << std::left << setw(20) << "Number of instances:" << std::right << setw(10) << numInstances << endl
         << std::left << setw(20) << "Constraint value:" << std::right << setw(10) << tau << endl
         << std::left << setw(20) << "Number of items:" << std::right << setw(10) << numItem << endl
         << std::left << setw(20) << "Minimum score width:" << std::right << setw(10) << minWidth << endl
         << std::left << setw(20) << "Maximum score width:" << std::right << setw(10) << maxWidth << endl
         << std::left << setw(20) << "Minimum item width:" << std::right << setw(10) << minItemWidth << endl
         << std::left << setw(20) << "Maxmimum item width:" << std::right << setw(10) << maxItemWidth << endl
         << std::left << setw(20) << "Length of strips:" << std::right << setw(10) << stripLength << endl
         << std::left << setw(20) << "Population size:" << std::right << setw(10) << numPop << endl;
    if(xOver == 1){
        cout << std::left << setw(20) << "Recombination operator: " << std::right << setw(6) << "GGA" << endl;
    }
    else if(xOver == 2){
        cout << std::left << setw(20) << "Recombination operator: " << std::right << setw(6) << "GPX" << endl;
    }
    cout << std::left << setw(20) << "Random seed:" << std::right << setw(10) << randomSeed << endl;
    cout << "------------------------------\n\n";
}

int main(int argc, char **argv){
    if(argc <= 1){
       ProgramInfo();
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
    int stripLength = 5000;
    int numPop = 0;
    int xOver = 1;
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
            stripLength = atoi(argv[++x]);
        }
        else if(strcmp("-p", argv[x]) == 0){
            numPop = atoi(argv[++x]);
        }
        else if(strcmp("-r", argv[x]) == 0){
            xOver = atoi(argv[++x]);
        }
        else if(strcmp("-s", argv[x]) == 0){
            randomSeed = atoi(argv[++x]);
        }
    }
    //endregion

    ArgumentCheck(numInstances, tau, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, stripLength, numPop, xOver, randomSeed);

    int i, j, k, bestStart, bestEnd;
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
    vector<vector<int> > bestSolnStart;
    vector<int> bestSolnStartSum;
    double bestFitness = 0.0;
    double tempFitness;

    srand(randomSeed);

    time_t startTime, endTime;
    startTime = clock();

    CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores,
                   partners, adjMatrix, itemWidths, allItems);


    int LB = LowerBound(stripLength, totalItemWidth);
    cout << "Lower bound: " << LB << " strips." << endl << endl;

    CreateInitPop(tau, numPop, numScores, numItem, maxItemWidth, stripLength, allScores, partners, adjMatrix,
                  itemWidths, populationSum, population);


    //Finding the solution in the population that has the best fitness value
    for(i = 0; i < population.size(); ++i){
        tempFitness = Fitness(stripLength, populationSum[i], population[i]);
        if(tempFitness > bestFitness){
            bestFitness = tempFitness;
            bestStart = i;
        }
    }

    bestSolnStart = population[bestStart];
    bestSolnStartSum = populationSum[bestStart];


    cout << "START - Best solution in the population:\n"
         << "Solution: " << bestStart << "\nFitness: " << bestFitness << "\nSize: " << bestSolnStart.size() << " strips." << endl << endl;

    bestEnd = bestStart;


    for(instance = 0; instance < numInstances; ++instance) {
        EA(tau, xOver, numScores, maxItemWidth, stripLength, bestEnd, bestFitness, allScores, partners, adjMatrix, itemWidths, populationSum, population);
    }


    cout << "END - Best solution in the population:\n"
         << "Solution: " << bestEnd << "\nFitness: " << bestFitness << "\nSize: " << population[bestEnd].size() << " strips." << endl << endl;



    endTime = clock();
    double totalTimeMS = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    double totalTimeS = totalTimeMS / 1000;
    cout << "\nCPU Time = " << totalTimeMS << " milliseconds (" << totalTimeS << " seconds.)\nEND.\n";


}//END INT MAIN


