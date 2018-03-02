/*--------------/
ALH
mainexact.cpp
11/12/2017
/--------------*/
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
using namespace std;

#include "base.h"
#include "packing.h"

void programInfo(){

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
         << "       -A          [Basic approximate FFD.]\n"
         << "       -B          [Pack strips in turn, choosing smallest feasible score width.]\n"
         << "       -C          [FFD combined with AHCA (exact algorithm). Default.]\n"
         << "       -s <int>    [Random seed. Default = 1.]";
}

void argumentCheck(int tau, int minWidth, int maxWidth, int minItemWidth, int maxItemWidth, int maxStripWidth){
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
}

int main(int argc, char **argv){
    if(argc < 1){
        programInfo();
        exit(1);
    }
    cout << "SCSSP\n-------------\n";

    int x;
    int numInstances = 1000;
    int tau = 70;
    int numItem = 500;
    int minWidth = 1;
    int maxWidth = 70;
    int minItemWidth = 150;
    int maxItemWidth = 1000;
    int maxStripWidth = 5000;
    int packType = 3;
    int randomSeed = 1;

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
        else if(strcmp("-A", argv[x]) == 0){
            packType = 1;
        }
        else if(strcmp("-B", argv[x]) == 0){
            packType = 2;
        }
        else if(strcmp("-C", argv[x]) == 0){
            packType = 3;
        }
        else if(strcmp("-s", argv[x]) == 0){
            randomSeed = atoi(argv[++x]);
        }
    }


    //Variables
    /*int numInstances = atoi(argv[1]);
    int tau = atoi(argv[2]);
    int numItem = atoi(argv[3]);
    int minWidth = atoi(argv[4]);
    int maxWidth = atoi(argv[5]);
    int minItemWidth = atoi(argv[6]);
    int maxItemWidth = atoi(argv[7]);
    int maxStripWidth = atoi(argv[8]);
    int packType = atoi(argv[9]);
    int randomSeed = atoi(argv[10]);*/

    argumentCheck(tau, minWidth, maxWidth, minItemWidth, maxItemWidth, maxStripWidth);

    cout << "Number of instances:\t" << numInstances << endl
         << "Constraint value:\t" << tau << endl
         << "Number of items:\t" << numItem << endl
         << "Minimum score width:\t" << minWidth << endl
         << "Maximum score width:\t" << maxWidth << endl
         << "Minimum item width:\t" << minItemWidth << endl
         << "Maxmimum item width:\t" << maxItemWidth << endl
         << "Length of strips:\t" << maxStripWidth << endl;
    if(packType == 1){
        cout << "Heuristic:\tBasic approximate FFD.\n";
    }
    else if(packType == 2){
        cout << "Heuristic:\tPack strips in turn, choosing smallest feasible score width.\n";
    }
    else if(packType == 3){
        cout << "Heuristic:\tFFD combined with AHCA (exact algorithm).\n";
    }
    cout << "Random seed:\t" << randomSeed << endl;

    int i, j, instance, choice;
    int opt = 0, opt90 = 0, opt80= 0, opt70 = 0, opt60 = 0, opt50 = 0, optLow = 0;
    int numScores = numItem * 2;
    double totalItemWidth = 0.0;
    vector<int> allScores;
    vector<int> stripSum(numItem, 0);
    vector<int> partners(numScores, 0);
    vector<vector<int> > strip(numItem);
    vector<vector<int> > itemWidths(numScores, vector<int>(numScores, 0));
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0));
    srand(randomSeed);
    int cp = 0;
    int na = 0;
    int type0 = 0;
    int type1 = 0;
    int type2 = 0;
    int type3 = 0;
    //srand(unsigned(time(NULL))); //seed


    time_t startTime, endTime;
    startTime = clock();


    switch(packType){
        case 1:
            cout << "FFD Approx:\n";
            for(instance = 0; instance < numInstances; ++instance){
                createInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                basicFFD(opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, numItem, maxItemWidth, maxStripWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                resetVectors(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
            }
            break;

        case 2:
            cout << "FFD Smallest:\n";
            for(instance = 0; instance < numInstances; ++instance){
                createInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                pairSmallest(opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, maxStripWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                resetVectors(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
            }
            break;

        case 3:
            cout << "FFD Exact:\n";
            for(instance = 0; instance < numInstances; ++instance){
                createInstance(tau, numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths);
                FFDincAHCA(cp, na, type0, type1, type2, type3, instance, tau, opt, opt90, opt80, opt70, opt60, opt50, optLow, numScores, numItem, maxItemWidth, maxStripWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
                resetVectors(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
            }
            break;

        default:
            cout << "No packType selected\n\n";
            break;

    }

    resetVectors(numScores, numItem, allScores, partners, adjMatrix, itemWidths, stripSum, strip);
    output(opt, opt90, opt80, opt70, opt60, opt50, optLow, numInstances);
    cout << endl << endl;
    cout << "Number of times CP used " << cp << endl;
    cout << "Number of times noConnect " << na << endl;
    cout << "Number of times type 0 used " << type0 << endl;
    cout << "Number of times type 1 used " << type1 << endl;
    cout << "Number of times type 2 used " << type2 << endl;
    cout << "Number of times type 3 used " << type3 << endl;

    endTime = clock();
    double totalTime = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    cout << "\nCPU Time = " << totalTime << " milliseconds.\nEND.\n";


}//END INT MAIN


