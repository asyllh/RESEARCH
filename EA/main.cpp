/*--------------/
ALH
main.cpp
Evolutionary Algorithm with Local Search
05/12/2017
/--------------*/
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

#include "base.h"
#include "packing.h"

int main(int argc, char **argv){
    //region USAGE - ARGUMENTS REQUIRED
    if(argc < 8){
        cout << "Minimum Score Separation Problem: AHCA.\n";
        cout << "Arguments are the following:\n";
        cout << "- Number of instances (integer)\n";
        cout << "- Number of items (integer)\n";
        cout << "- Minimum width of scores (millimeters, min = 1)\n";
        cout << "- Maximum width of scores (millimeters, max = 70)\n";
        cout << "- Minimum width of boxes (millimeters, min = 140)\n";
        cout << "- Maximum width of boxes (millimeters, max = 1000)\n";
        cout << "- Maximum width of strips (millimeters)\n";
        cout << "- Random Seed (integer)\n";
        cout << "- Name of input file (must have .txt extension)\n";
        exit(1);
    }
    //endregion

    //region VARIABLES
    //int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
    int numItem = atoi(argv[1]); //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    int minWidth = atoi(argv[2]); //minimum width of scores (millimeters)
    int maxWidth = atoi(argv[3]); //maximum width of scores (millimeters)
    int minItemWidth = atoi(argv[4]); //min box width (mm)
    int maxItemWidth = atoi(argv[5]); //max box width (mm)
    int maxStripWidth = atoi(argv[6]);
    int randomSeed = atoi(argv[7]); //random seed

    int i, j, q, n;
    int instance; //counter for instances loop
    int numScores = numItem * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    //int threshold = 70; //adjacency threshold of scores, minimum knife distance
    vector<int> allScores;
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0)); //adjaceny matrix (createInstance, MTGMA, MIS, FCA)
    vector<int> partners(numScores, 0); //contains vertex index for partners, e.g if vertex 0 is partners with vertex 4, then partners[0] = 4 (createInstance, MIS)
    vector<vector<int> > itemWidths(numScores, vector<int>(numScores, 0)); // holds widths in mm for the widths of boxes
    vector<vector<int> > allItems(numScores, vector<int>(numScores, 0)); //contains values from 1 to numItem to denote which box the scores belong to
    double totalItemWidth = 0.0;
    //vector<int> stripSum(numItem, 0);
    //vector<vector<int> > strip(numItem);
    vector<vector<vector<int> > > population;
    vector<vector<int> > populationSum;
    double parent1cost;
    double parent2cost;
    //endregion



    cout << "local search\n-------------\n";

    //srand(unsigned(time(NULL))); //seed
    srand(randomSeed);


    time_t startTime, endTime; //start clock
    startTime = clock();

    createInstance(numScores, numItem, minWidth, maxWidth, minItemWidth, maxItemWidth, totalItemWidth, allScores, partners, adjMatrix, itemWidths, allItems);

    createInitialPopulation(numScores, numItem, maxItemWidth, maxStripWidth, allScores, partners, adjMatrix, itemWidths, populationSum, population);

    EA(numScores, maxItemWidth, maxStripWidth, parent1cost, parent2cost, allScores, partners, adjMatrix, itemWidths, populationSum, population);


    endTime = clock();
    double totalTime = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    cout << "\nCPU Time = " << totalTime << " milliseconds.\nEND.\n";


}//END INT MAIN


