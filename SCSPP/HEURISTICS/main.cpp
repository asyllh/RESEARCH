/*--------------/
ALH
mainexact.cpp
11/12/2017
/--------------*/
#include <iostream>
#include <fstream>
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
         << "       -m <int>    [Minimum item width. Default = 150.]\n"
         << "       -M <int>    [Maximum item width. Default = 1000.]\n"
         << "       -W <int>    [Width of strips. Default = 5000.]\n"
         << "       -s <int>    [Random seed. Default = 1.]\n"
         << "---------------\n"
         << "ALGORITHM:\n"
         << "       -x <int>    [1: Basic approximate FFD.]\n"
         << "                   [2: Pack strips in turn, choosing smallest feasible score width.]\n"
         << "                   [3: FFD combined with AHCA (exact algorithm).]\n"
         << "---------------\n\n";

}

void ArgumentCheck(int numInstances, int tau, int numItem, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth,
                   int stripWidth, int algType, int randomSeed){

    bool error = false;

    cout << "SCSPP\n------------------------------\n";
    if(tau == 0){
        //cout << "[ERROR]: Constraint value cannot be zero.\n";
        cout << "[WARNING]: Constraint value is zero, problem instances are equivalent to classical SPP without score constraints.\n";
        //error = true;
    }
    if(stripWidth == 0){
        cout << "[ERROR]: Strip cannot have width zero.\n";
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
    if(maxItemWidth > stripWidth){
        cout << "[ERROR]: Maximum item width is larger than width of strip.\n";
        error = true;
    }
    if(algType == 0){
        cout << "[ERROR]: No algorithm specified.\n";
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
         << std::left << setw(20) << "Width of strips:" << std::right << setw(10) << stripWidth << endl;
    if(algType == 1){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "MFFD\n";
    }
    else if(algType == 2){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "PairSmallest\n";
    }
    else if(algType == 3){
        cout << std::left << setw(18) << "Algorithm:" << std::right << setw(13) << "MFFD+\n";
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
    int stripWidth = 5000;
    int algType = 0;
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
        else if(strcmp("-m", argv[x]) == 0){
            minItemWidth = atoi(argv[++x]);
        }
        else if(strcmp("-M", argv[x]) == 0){
            maxItemWidth = atoi(argv[++x]);
        }
        else if(strcmp("-W", argv[x]) == 0){
            stripWidth = atoi(argv[++x]);
        }
        else if(strcmp("-x", argv[x]) == 0){
            algType = atoi(argv[++x]);
        }
        else if(strcmp("-s", argv[x]) == 0){
            randomSeed = atoi(argv[++x]);
        }
    }
    //endregion

    ArgumentCheck(numInstances, tau, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, stripWidth, algType, randomSeed);

    int a, i, j, k, instance;
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
    vector<int> delta = { 145, 109, 97, 87, 79, 72, 65, 57, 47, 34, 0 };
    //double Delta = 0.0;

    Timer timer;

    switch(algType){
        case 1:
            for (a = 0; a < delta.size(); ++a) {
                tau = delta[a];
                srand(randomSeed);
                for (instance = 0; instance < numInstances; ++instance) {
                    CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                   totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                    MFFD(numScores, numItem, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                }
            }
            break;

        case 2:
            for (a = 0; a < delta.size(); ++a) {
                tau = delta[a];
                srand(randomSeed);
                for (instance = 0; instance < numInstances; ++instance) {
                    CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                   totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                    PairSmallest(numScores, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                }
            }
            break;

        case 3:
            for (a = 0; a < delta.size(); ++a) {
                tau = delta[a];
                srand(randomSeed);
                for (instance = 0; instance < numInstances; ++instance) {
                    CreateInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth,
                                   totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                    MFFDPlus(tau, numScores, numItem, maxItemWidth, stripWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                    ResetVar(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                }
            }
            break;

    }

}//END INT MAIN


