/*--------------/
ALH
main.cpp
12/10/17
/--------------*/
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

#include "base.h"
#include "packing.h"

int main(int argc, char **argv){
    //region USAGE - ARGUMENTS REQUIRED
    if(argc < 9){
        cout << "Minimum Score Separation Problem: MBAHRA.\n";
        cout << "Arguments are the following:\n";
        cout << "- Number of instances (integer)\n";
        cout << "- Number of boxes (integer)\n";
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

    //VARIABLES FROM ARGUMENTS
    //int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
    int numBox = atoi(argv[1]); //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    int minWidth = atoi(argv[2]); //minimum width of scores (millimeters)
    int maxWidth = atoi(argv[3]); //maximum width of scores (millimeters)
    int minBoxWidth = atoi(argv[4]); //min box width (mm)
    int maxBoxWidth = atoi(argv[5]); //max box width (mm)
    int maxStripWidth = atoi(argv[6]);
    int randomSeed = atoi(argv[7]); //random seed

    //VARIABLES
    int i, j, k, q, n;
    int instance; //counter for instances loop
    int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    int threshold = 70; //adjacency threshold of scores, minimum knife distance
    int vacant = 999;
    vector<int> allScores;
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0)); //adjaceny matrix (createInstance, MTGMA, MIS, FCA)
    vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4 (createInstance, MIS)
    vector<vector<int> > boxWidths(numScores, vector<int>(numScores, 0)); // holds widths in mm for the widths of boxes
    vector<vector<int> > allBoxes(numScores, vector<int>(numScores, 0)); //contains values from 1 to numBox to denote which box the scores belong to
    srand(randomSeed); //seed
    vector<vector<int> > userInput;
    vector<int> tempUser(3, 0);
    int choice;
    double totalBoxWidth = 0.0;
    vector<int> stripSum(numBox, 0);
    vector<vector<int> > strip(numBox);
    vector<vector<int> > stripWidth(numBox);
    vector<int> stripNumBoxes;
    int totalCost;
    int moved;


    cout << "MSSP - MBAHRA\n-------------\n";

    //region READ FILE
    ifstream inStream;
    inStream.open(argv[8]);
    if(inStream.fail()){
        cout << "ERROR: file cannot be opened.\n";
        exit(1);
    }
    inStream >> n;

    for(i = 0; i < n; ++i){
        for(j = 0; j < 3; ++j) {
            inStream >> tempUser[j];
        }
        userInput.push_back(tempUser);
    }
    inStream.close();

    //endregion

    time_t startTime, endTime; //start clock
    startTime = clock();

    //region READ FILE CHOICE

        /*cout << "Please choose from the following:\n";
        cout << "1: Create MSSP instance using random values\n";
        cout << "2: Create MSSP instance from file\n";
        cout << "Enter choice:  ";
        cin >> choice;
        while(choice !=1 && choice != 2){
            cout << "Please enter either 1 or 2:";
            cin >> choice;
        }
        switch(choice){
            case 1:
                createInstance(threshold, minWidth, maxWidth, minBoxWidth, maxBoxWidth, numScores, numBox, totalBoxWidth, allScores, adjMatrix, mates, boxWidths, allBoxes);
                break;

            case 2:
                createInstanceUser(threshold, numScores, totalBoxWidth, allScores, userInput, adjMatrix, mates, boxWidths, allBoxes);
                break;

            default:
                cout << "Please enter 1 or 2.\n";
                break;

        }*/

        //endregion

    //createInstance(threshold, minWidth, maxWidth, minBoxWidth, maxBoxWidth, numScores, numBox, totalBoxWidth, allScores, adjMatrix, mates, boxWidths, allBoxes);

    createInstanceUser(threshold, numScores, totalBoxWidth, allScores, userInput, adjMatrix, mates, boxWidths, allBoxes);

    packStripsFFD(totalCost, numBox, maxBoxWidth, maxStripWidth, totalBoxWidth, adjMatrix, mates, boxWidths, stripSum, stripNumBoxes, strip, stripWidth);

    checkSwap(totalCost, maxStripWidth, adjMatrix, boxWidths, strip, stripSum);

    checkMove(moved, maxStripWidth, adjMatrix, boxWidths, strip, stripSum);







    endTime = clock();
    double totalTime = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    cout << "CPU Time = " << totalTime << " milliseconds.\n";

}//END INT MAIN


