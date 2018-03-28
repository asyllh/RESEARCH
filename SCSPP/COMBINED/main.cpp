/*--------------/
ALH
main.cpp
Combined Program with Heuristics and EA
17/03/2018
/--------------*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <chrono>
#include <iomanip>
using namespace std;

#include "base.h"
#include "packing.h"

struct Timer{

    std::chrono::high_resolution_clock::time_point startTime, endTime;
    std::chrono::duration<float> totalTime;

    Timer(){
        startTime = std::chrono::high_resolution_clock::now();
    }

    ~Timer(){
        endTime = std::chrono::high_resolution_clock::now();
        totalTime = endTime - startTime;

        float totalTimems = totalTime.count() * 1000.0f;
        //cout << "\nCPU Time: " << totalTimems << "ms (" << totalTime.count() << "s)" << endl;
        cout << totalTime.count() << endl;
        //cout << totalTimems << endl;
    }
};

void ProgramInfo(){

    cout << "SCSPP:\n-------------\n"
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
         << "---------------\n"
         << "ALGORITHM:\n"
         << "       -x <int>    [1: Basic approximate FFD.]\n"
         << "                   [2: Pack strips in turn, choosing smallest feasible score width.]\n"
         << "                   [3: FFD combined with AHCA (exact algorithm).]\n"
         << "                   [4: Evolutionary Algorithm using Local Search.]\n"
         << "       -r <int>    [Recombination operator. 1: GGA. 2: GPX'.]\n"
         << "       -p <int>    [Number of solutions in population.]\n"
         << "---------------\n\n";

}

void ArgumentCheck(int numInstances, int tau, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth,
                   int stripLength, int algType, int numPop, int xOver, int randomSeed){

    bool error = false;

    cout << "SCSPP\n------------------------------\n";
    if(tau == 0){
        //cout << "[ERROR]: Constraint value cannot be zero.\n";
        cout << "[WARNING]: Constraint value is zero, problem instances are equivalent to classical SPP without score constraints.\n";
        //error = true;
    }
    if(stripLength == 0){
        cout << "[ERROR]: Strip cannot have length zero.\n";
        error = true;
    }
    if(2*minWidth >= tau){
        //cout << "[ERROR]: Constraint value is less than or equal to twice the minimum score width, vicinal sum constraint always valid.\n";
        cout << "[WARNING]: Constraint value is less than or equal to twice the minimum score width, vicinal sum constraint always valid.\n";
        cout << "         Problem instance is therefore classical strip-packing problem without score constraint (i.e. tau = 0).\n";
        //error = true;
    }
    if(2*maxWidth < tau){
        //cout << "[ERROR]: Constraint value is greater than double maximum score width, vicinal sum constraint never valid.\n";
        cout << "[WARNING]: Constraint value is greater than double maximum score width, vicinal sum constraint never valid.\n";
        cout << "           Number of strips required = number of items.\n";
        //error = true;
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
    if(algType == 4 && numPop < 5){
        cout << "[ERROR]: Insufficient number of solutions in population.\n";
        error = true;
    }
    if(algType == 0){
        cout << "[ERROR]: No algorithm specified.\n";
        error = true;
    }
    if(algType == 4 && (xOver != 1 && xOver != 2)){
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
         << std::left << setw(20) << "Length of strips:" << std::right << setw(10) << stripLength << endl;
    if(algType == 1){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "BasicFFD\n";
    }
    else if(algType == 2){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "PairSmallest\n";
    }
    else if(algType == 3){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "FFDincAHCA\n";
    }
    else if(algType == 4  && xOver == 1){
        cout << std::left << setw(19) << "Algorithm: " << std::right << setw(6) << "EA with GGA" << endl;
        cout << std::left << setw(20) << "Population size:" << std::right << setw(10) << numPop << endl;
    }
    else if(algType == 4 && xOver == 2){
        cout << std::left << setw(18) << "Algorithm: " << std::right << setw(6) << "EA with GPX'" << endl;
        cout << std::left << setw(20) << "Population size:" << std::right << setw(10) << numPop << endl;
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
    int algType = 0;
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
        else if(strcmp("-x", argv[x]) == 0){
            algType = atoi(argv[++x]);
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

    ArgumentCheck(numInstances, tau, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, stripLength, algType, numPop, xOver, randomSeed);

    int a, i, j, k, instance, bestStart, bestEnd;
    int opt = 0, opt90 = 0, opt80= 0, opt70 = 0, opt60 = 0, opt50 = 0, optLow = 0;
    int numScores = numItem * 2;
    double totalItemWidth = 0.0;
    vector<int> allScores;
    vector<int> stripSum(numItem, 0);
    vector<int> partners(numScores, 0);
    vector<vector<int> > strip(numItem);
    vector<vector<int> > itemWidths(numScores, vector<int>(numScores, 0));
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0));
    vector<vector<vector<int> > > population;
    vector<vector<int> > populationSum;
    vector<vector<int> > bestSolnStart;
    vector<int> bestSolnStartSum;
    double bestFitness = 0.0;
    double tempFitness;
    vector<int> alpha = { 145, 109, 97, 87, 79, 72, 65, 57, 47, 34, 0 };

    srand(randomSeed);

    //Timer timer;

    switch(algType){
        case 1:
            for(a = 0; a < alpha.size(); ++a) {
                tau = alpha[a];
                {
                    Timer timer;
                    for (instance = 0; instance < numInstances; ++instance) {
                        CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                       totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                        BasicFFD(opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, numItem, maxItemWidth,
                                 stripLength, totalItemWidth, allScores, partners, adjMatrix, itemWidths, stripSum,
                                 strip);
                        ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    }
                }
                //Output(opt, opt90, opt80, opt70, opt60, opt50, optLow, numInstances);
            }
            break;

        case 2:
            for(a = 0; a < alpha.size(); ++a) {
                tau = alpha[a];
                {
                    Timer timer;
                    for (instance = 0; instance < numInstances; ++instance) {
                        CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                       totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                        PairSmallest(opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, numItem, stripLength,
                                     totalItemWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                        ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    }
                }
                //Output(opt, opt90, opt80, opt70, opt60, opt50, optLow, numInstances);
            }
            break;

        case 3:
            for(a = 0; a < alpha.size(); ++a) {
                tau = alpha[a];
                {
                    Timer timer;
                    for (instance = 0; instance < numInstances; ++instance) {
                        CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                       totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                        FFDincAHCA(tau, opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, numItem,
                                   maxItemWidth, stripLength, totalItemWidth, allScores, partners, adjMatrix,
                                   itemWidths, stripSum, strip);
                        ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    }
                }
                //Output(opt, opt90, opt80, opt70, opt60, opt50, optLow, numInstances);
            }
            break;

        case 4:
            CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths);
            int LB = LowerBound(stripLength, totalItemWidth);
            cout << "Lower bound: " << LB << " strips." << endl << endl;
            CreateInitPop(tau, numPop, numScores, numItem, maxItemWidth, stripLength, allScores, partners, adjMatrix, itemWidths, populationSum, population);

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
            bestEnd = bestStart;

            cout << "START - Best solution in the population:\n"
                 << "Solution: " << bestStart << "\nFitness: " << bestFitness << "\nSize: " << bestSolnStart.size() << " strips." << endl << endl;

            for(instance = 0; instance < numInstances; ++instance) {
                EA(tau, xOver, numScores, maxItemWidth, stripLength, bestEnd, bestFitness, allScores, partners, adjMatrix, itemWidths, populationSum, population);
            }

            cout << "END - Best solution in the population:\n"
                 << "Solution: " << bestEnd << "\nFitness: " << bestFitness << "\nSize: " << population[bestEnd].size() << " strips." << endl << endl;
            break;

    }



}//END INT MAIN