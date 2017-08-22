/*--------------------/
ALH
mbahra.h
18/08/2017
/--------------------*/
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

#include "mbahra.h"
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
        exit(1);
    }
    //endregion

    //region VARIABLES
    //VARIABLES FROM ARGUMENTS
    int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
    int numBox = atoi(argv[2]) + 1; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    int minWidth = atoi(argv[3]); //minimum width of scores (millimeters)
    int maxWidth = atoi(argv[4]); //maximum width of scores (millimeters)
    int minBoxWidth = atoi(argv[5]); //min box width (mm)
    int maxBoxWidth = atoi(argv[6]); //max box width (mm)
    int maxStripWidth = atoi(argv[7]);
    int randomSeed = atoi(argv[8]); //random seed


    //VARIABLES
    int i, j, k, q, n;
    int instance; //counter for instances loop
    int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    int numComp = (numBox + (numBox % 2)) / 2;
    int threshold = 70; //adjacency threshold of scores, minimum knife distance
    int vacant = 999; //large empty value
    int feasible = 0; //number of feasible instances
    int infeasible = 0; //number of infeasible instances
    int noMatch = 0; //number of instances with |M| < n (weak match, therefore immediately infeasible)
    int oneCycle = 0; //number of instances where the MIS consists of only one cycle (therefore immediately feasible)
    int noFam = 0; //number of instances with no family of T-cycles (qstar = -1, therefore immediately infeasible)
    int noPatch = 0; //number of instances where T-cycles do not produce a connected patching graph (SSum < numCycles, therefore infeasible)
    int fullT = 0; //number of instances where only one T-cycle is required to connected all cycles in MIS (patching graph is connected using only one T-cycle, therefore feasible)
    int splitT = 0; //number of instances where multiple T-cycles are required to connected all cycles in MIS (patching graph is connected using multiple T-cycles, therefore feasible)
    int problem = 0; //number of problematic instances, SSum > numCycles (ERROR)
    //int startPath = vacant;
    //int endPath = vacant;

    int matchSize; //size (cardinality) of the matching list (matchList.size()) (&MTGMA, FCA)
    int numCycles; //number of cycles in the MIS (mateInduced.size()) (&MIS, patchGraph)
    int qstar; //number of T-cycles (&FCA, patchGraph)
    //vector<int> allScores(numScores, 0); //vector containing all score widths (createInstance, MTGMA)
    vector<int> allScores;
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0)); //adjaceny matrix (createInstance, MTGMA, MIS, FCA)
    vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4 (createInstance, MIS)
    vector<int> matchList(numScores, vacant); //contains vertex index for matching vertices, e.g. if vertex 0 is matched with vertex 9, then matchList[0] = 9 (MTGMA, MIS, FCA, patchGraph)
    vector<int> cycleVertex(numScores, 1); //contains the number of the cycle of the mate-induced structure that the vertex i belongs to (MTGMA, FCA, patchGraph)
    vector<vector<int> > mateInduced; //each row of the matrix corresponds to one cycle, and contains the indices of the score widths (MIS, FCA, patchGraph)
    vector<int> lengthMateInduced; //each elements holds the value of the length of the corresponding cycle in the MIS (MIS)
    vector<vector<int> > S(numComp, vector<int>(numComp, 0)); // == 1 if edge from cycle j is used in T-cycle q (T[q][j]) (FCA, patchGraph)
    vector<vector<int> > T; //each row hold the lower vertex of the edges that make up one T-cycle
    vector<vector<int> > boxWidths(numScores, vector<int>(numScores, 0)); // holds widths in mm for the widths of boxes
    vector<int> fullCycle; //holds initial cycle of final solution
    vector<int> completePath; //holds final path, i.e. fullCycle without dominating vertices
    vector<vector<int> > allBoxes(numScores, vector<int>(numScores, 0)); //contains values from 1 to numBox to denote which box the scores belong to
    //if value is positive, i.e "3", then the 3rd box has the smallest score on the LHS, and the larger score on the RHS
    //if value is negative, i.e "-3", then the 3rd box is rotated, i.e. smallest score on RHS and larger score on LHS.
    srand(randomSeed); //seed
    vector<vector<int> > userInput;
    vector<int> tempUser(3, 0);
    int choice;
    //endregion

    cout << "MSSP - MBAHRA\n-------------\n";

    //region READ FILE
    ifstream inStream;
    inStream.open(argv[9]);
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

    for(instance = 0; instance < numInstances; ++instance) {
        //resetVectors(vacant, numScores, numComp, allScores, adjMatrix, cycleVertex, matchList, mates, S, boxWidths, allBoxes);
        //clearVectors(allScores, mateInduced, lengthMateInduced, T, fullCycle, completePath);

        //region READ FILE CHOICE
        cout << "Please choose from the following:\n";
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
                createInstance(threshold, minWidth, maxWidth, minBoxWidth, maxBoxWidth, numScores, numBox, allScores, adjMatrix, mates, boxWidths, allBoxes);
                break;

            case 2:
                createInstanceUser(threshold, numScores, allScores, userInput, adjMatrix, mates, boxWidths, allBoxes);
                break;

            default:
                cout << "Please enter 1 or 2.\n";
                break;

        }
        //endregion
        //continue;

        //packStripsBFD(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsFFD(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsNFD(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsBFI(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsFFI(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsNFI(numBox, maxBoxWidth, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsFFDScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsFFIScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsBFDScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsBFIScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);
        //packStripsNFDScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);
        packStripsNFIScores(vacant, numBox, maxStripWidth, adjMatrix, mates, boxWidths);

        continue;

        //packStripsSmallest(numScores, numBox, maxStripWidth, mates, adjMatrix, boxWidths);
        //continue; //do not do MTGMA/MIS/FCA/PATCH

        //region MTGMA
        //MTGMA(vacant, threshold, numScores, matchSize, allScores, adjMatrix, cycleVertex, matchList);
        //continue;
        //If the number of matches (i.e. the size of the matching list M) is less than the number of boxes (n), then instance is infeasible ( |M| < n )
        //continue;
        /*if(matchSize == numBox - 1){
            weakMatchPath(vacant, numScores, matchList, mates);
        }*/
        //continue;
        if (matchSize < numBox) {
            cout << instance << ": INFEASIBLE - Not enough matching edges.\n\n";
            ++infeasible;
            ++noMatch;
            continue;
        }
        //endregion

        //region MIS
        MIS(numScores, numCycles, adjMatrix, mates, matchList, mateInduced, lengthMateInduced);
        //packStripsMIS(numBox, maxStripWidth, adjMatrix, mateInduced, boxWidths);
        //continue;
        //If the mate-induced structure only consists of one cycle, then the problem has been solved and is feasible (just remove one matching edge to find feasible path)
        if (lengthMateInduced[0] == numScores) { //if all of the vertices are in the first (and only) cycle of the mate-induced structure
            for(j = 0; j < mateInduced[0].size(); ++j){
                fullCycle.push_back(mateInduced[0][j]);
            }
            makePath(numScores, fullCycle, completePath, boxWidths, allScores, allBoxes);
            ++feasible;
            ++oneCycle;
            continue;
        }
        //endregion

        //region FCA
        FCA(qstar, vacant, matchSize, adjMatrix, cycleVertex, matchList, mateInduced, S, T);
        //If no family of T-cycle found
        if (qstar == -1) {
            cout << instance << ": Infeasible, qstar = -1, no family of T-cycles found.\n\n";
            ++infeasible;
            ++noFam;
            continue;
        }
        //endregion

        //region patchGraph
        //Check if patching graph is connected
        patchGraph(qstar, vacant, instance, numScores, numCycles, feasible, infeasible, fullT, splitT, noPatch, problem, matchList, cycleVertex, mateInduced, S, T, fullCycle, completePath, boxWidths, allScores, allBoxes);
        //endregion

    } //end of for loop instances


    //region OUTPUT
    /*cout << "------------------------------------------------------------------\n";
    cout << "INPUT:\n";
    cout << "# of instances: " << numInstances << endl;
    cout << "# of boxes: " << numBox - 1 << endl;
    cout << "Min score width: " << minWidth << "mm\n";
    cout << "Max score width: " << maxWidth << "mm\n";
    cout << "Min box width: " << minBoxWidth << "mm\n";
    cout << "Max box width: " << maxBoxWidth << "mm\n";
    cout << "Max strip width: " << maxStripWidth << "mm\n";
    cout << "Random seed: " << randomSeed << endl;
    cout << "# of scores: " << numScores << endl;
    cout << "Threshold: " << threshold << "mm\n\n";

    cout << "EVALUATION:\n";
    cout << "# feasible instances: " << feasible << endl;
    cout << "# infeasible instances: " << infeasible << endl;
    cout << "# instances that did not have enough matching edges (|M| < n) (I): " << noMatch << endl;
    cout << "# instances where MIS consisted of one complete cycle (F): " << oneCycle << endl;
    cout << "# instances where no T-cycles were found (I): " << noFam << endl;
    cout << "# instances that required one T-cycle to create solution (F): " << fullT << endl;
    cout << "# instances that required multiple T-cycles to create solution (F): " << splitT << endl;
    cout << "# instances where patching graph was unconnected (I):  " << noPatch << endl;
    cout << "# instances that are problematic: " << problem << endl;
    */
    //endregion
    endTime = clock();
    double totalTime = (((endTime - startTime) / double(CLOCKS_PER_SEC)) * 1000);
    cout << "CPU Time = " << totalTime << " milliseconds.\n";

}//END INT MAIN


/* TO USE PACK STRIPS FUNCTIONS:
 *
 * change variable vector<int> matchList(numScores, vacant) to vector<int> matchList(numScores -2, vacant);
 * void resetVector: for(j = 0; j < numScores -2; ++j) { matchList[i] = vacant }
 * void createInstance: change adjMatrix output from i < numScores, j< numScores to i < numScores -2, j < numScores -2
 * void MIS: change vector<int> checked(numScores, 0) to vector<int> checked(numScores -2, 0)
 * void MTGMA: change for(i = 0; i < numScores; ++i) to for(i =0; i < numScores -2; ++i)
 * void MTGMA: change for(j = numScores-1; j > i; --j) to for(j = numScores - 3; j > i; --j)
 * after MTGMA: change if(matchSize < numBox) to if(matchSize < numBox -1)
 *
 */
