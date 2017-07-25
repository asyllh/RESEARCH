/*--------------------/
ALH
Rewritten Minimum Score Separation Problem
17/07/2017
/--------------------*/

#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <math.h>
#include <float.h>
#include <queue>
#include <limits.h>
#include <algorithm>
#include <iomanip> //header providing parametric manipulators
using namespace std;


int main(){

    //Variables
    int i, j, k, r;
    unsigned int randomSeed = 3;
    int numInstances = 100; //number of instances of mssp, use in main for loop
    unsigned int numBox = 5; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    unsigned int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    int minWidth = 1; //minimum width of score (millimeters)
    int maxWidth = 60; //maximum width of score (millimeters)
    int threshold = 70; //adjacency threshold of scores, minimum knife distance
    int mateMatch; //vertex number for matching algorithm, mate takes the value of the index of the vertex that the current vertex is mates with
    int lastMatch;
    int matchSize; //size/cardinality of the matching list
    int vacant = 999;
    int verticesNotMatched; //counts number of vertices that have not been matched to another vertex in MTGMA
    int smallestVertex;
    int currentVertex;
    int totalCycles; //number of cycles in the mate-induced structure, i.e mateInduced.size() (each row in the mate-induced structure matrix represents a cycle)
    int feasible = 0; //number of feasible instances
    int infeasible = 0; //number of infeasible instances
    int currentEdge; //counter used for list of edges
    int numEdges; //number of (non-empty) edges

    vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4
    vector<int> matchList(numScores, 0); //contains vertex index for matching vertices, e.g. if vertex 0 is matched with vertex 9, then matchList[0] = 9
    vector<int> allScores(numScores, 0); //vector containing all score widths
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores ,0)); //adjacency matrix, 0 if width sum < threshold, 1 if width sum >= threshold, 2 if scores are mates (either side of same box)
    vector<int> checked(numScores, 0); //contains 0 if vector i has not yet been included in the mate-induced structure, 1 if vector i has been placed in MIS
    vector<vector<int> > mateInduced; //size numscore by noComp, i.e. number of rows = numScores, number of columns = noComp
    vector<int> cycle; //used in building the mate-induced structure
    vector<int> lengthMateInduced; //each element holds the value corresponding to the length of the relative cycle in the mate-induced structure
    vector<int> cycleVertex(numScores, 0); //contains the number of the cycle of the mate-induced structure that the vertex i belongs to
    // (e.g. if vertex 4 is in the first cycle of the MIS, then cycleVertex[4] = 0 (0 = first cycle, 1 = second cycle etc))
    vector<int> edge(numBox, 0); //contains the number of lower vertex of each (non-empty) edge

    srand(randomSeed); //seed

    //DON'T FORGET TO CLEAR VECTORS AND VARIABLES FOR NEXT INSTANCE

    cout << "Minimum Score Separation Problem\nMatching-Based Alternating Hamiltonicity Recognition Algorithm\n\n";

    time_t startTime, endTime; //start clock
    startTime = clock();

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for(i = 0; i < numScores -2; ++i){
        allScores[i] = rand() % (maxWidth - minWidth + 1) + minWidth;
    }
    //add two dominating vertices with score widths = 71 (these scores will be either side of same box, mates)
    allScores[numScores - 2] = 71;
    allScores[numScores - 1] = 71;

    //Print out allScores vector
    cout << "All scores:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << " ";
    }
    cout << endl << endl;

    //Sort all of the scores in the allScores vector in ascending order
    sort (allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    //Print out allScores vector (scores now in ascending order)
    cout << "All scores - non-decreasing order:\n";
    for(i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << " ";
    }
    cout << endl << endl;

    //Filling in adjacency matrix - if sum of two scores >= threshold (70), then insert 1 into the matrix, else leave as 0
    for(i = 0; i < allScores.size()-1; ++i){
        for(j = i+1; j< allScores.size(); ++j){
            if(allScores[i] + allScores[j] >= threshold){
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    //Create vector to be used to assign mates
    vector<int> randOrder(numScores, 0);

    //Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for(i = 0; i < numScores-2; ++i){
        randOrder[i] = i;
    }
    randOrder[numScores-2] = numScores - 2;
    randOrder[numScores-1] = numScores - 1;
    //Randomly shuffle all values in randOrder vector EXCEPT the last two values (dominating vertices, must stay as mates)
    random_shuffle(randOrder.begin(), randOrder.begin()+8);

    //Print out randOrder vector
    /*cout << "Random Order:\n";
    for(i = 0; i < randOrder.size(); ++i){
        cout << randOrder[i] << endl;
    }
    cout << endl;*/

    //Assign mates to each score (i.e. pair up scores to define which scores are either side of the same box)
    //In the adjacency matrix, this will be represented by value 2
    //Therefore there will be a value of 2 in every row and every column, non repeating
    for(i = 0; i < numBox; ++i){
        adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
        adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
    }

    //Print out adjacency matrix inc threshold and mates
    cout << "Adjacency Matrix:\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;


    //MATCHING ALGORITHM MTGMA
    //Fill matchingList vector with values 0,..., numScores-1 (i.e. the index of each element)
    for(i = 0; i < numScores; ++i){
        matchList[i] = vacant;
    }
    matchSize = 0;
    lastMatch = vacant;
    verticesNotMatched = 0;

    for(i = 0; i < numScores; ++i){ //check all vertices
        if(matchList[i] == vacant){ //if vertex has not yet been matched
            for(j = numScores -1; j > i; --j){ //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
                if(adjMatrix[i][j] == 1 && matchList[j] == vacant){ //if vertices i and j are adjacent, and if vertex j has not yet been matched
                    matchList[i] = j;
                    matchList[j] = i;
                    lastMatch = i;
                    ++matchSize;
                    break;
                }
                else if(adjMatrix[i][j] == 2 && matchList[j] == vacant){ //if potential match == mate
                    // mark for FCA
                }
            }//end for j
            if(matchList[i] == vacant){ //if vertex has still not been matched
                for(k = 0; k < numScores; ++k){
                    if(adjMatrix[i][k] == 2){ //if vertex i and vertex k are mates
                        mateMatch = k;
                        break;
                    }
                }
                if((allScores[i] + allScores[mateMatch] >= threshold) //match with mate?
                    && (matchList[mateMatch] == vacant) //is mate unmatched?
                    && (lastMatch != vacant) //has the previous vertex been matched?
                    && (mateMatch > i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
                    && (allScores[lastMatch] + allScores[mateMatch] >= threshold)){ //can mate be matched with last matched vertex?
                    // if so, then swap mates
                    matchList[i] = matchList[lastMatch];
                    matchList[lastMatch] = mateMatch;
                    matchList[mateMatch] = lastMatch;
                    matchList[matchList[i]] = i;
                    lastMatch = i;
                    ++matchSize;
                }
                else{
                    //one more unconnected vertex
                }


            }//end if

        }//end if matchList[i] == i
    }//end for i

    //If the number of matches (i.e. the size of the matching list M) is less than the number of boxes (n), then instance is infeasible ( |M| < n )
    if(matchSize < numBox){
        ++infeasible;
        cout << "Instance is infeasible, not enough matching edges available (|M| < n)." << endl;
        goto End;
        //continue;
    }
    cout << "Size of M (matchSize): " << matchSize << endl;


    cout << "Matching List:\n";
    for(i = 0; i < numScores; ++i){
        cout << matchList[i] << " ";
    }
    cout << endl << endl;

    cout << "Matching Vertices Values:\n";
    for(i = 0; i < numScores; ++i){
        if(matchList[i] == vacant){
            cout << i << "\t" << allScores[i] << "\t" << "No Match" << endl;
            ++verticesNotMatched;
        }
        else {
            cout << i << "\t" << allScores[i] << "\t" << allScores[matchList[i]] << "\t" << matchList[i] << endl;
        }

    }

    if(verticesNotMatched == 0){
        cout << "All vertices have been matched.\n\n";
    }
    else {
        cout << "Number of unmatched vertices: " << verticesNotMatched << endl;
    }



    //MATE-INDUCED STRUCTURE
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                mates[i] = j;
                break;
            }
        }
    }
    cout << "Mates Vector:\n";
    for(i = 0; i < numScores; ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < numScores; ++i){
        checked[i] = 0;
    }

    //find the smallest vertex not yet checked for mate-induced structure - start with this vertex
    for(i = 0; i < numScores; ++i){
        if(checked[i] == 0){
            smallestVertex = i;
            break;
        }
    }

    //Building the mate-induced structure
    do{
        currentVertex = smallestVertex;
        do{
            cycle.push_back(currentVertex);
            checked[currentVertex] = 1;
            cycle.push_back(mates[currentVertex]);
            checked[mates[currentVertex]] = 1;
            currentVertex = matchList[mates[currentVertex]];
        } while(currentVertex != smallestVertex);

        mateInduced.push_back(cycle);
        cycle.clear();

        for(i = 0; i < numScores; ++i){
            if(checked[i] == 0){
                smallestVertex = i;
                break;
            }
        }


    } while(smallestVertex != currentVertex);
    cycle.clear(); //clear cycle vector again for next instance

    totalCycles = mateInduced.size();
    for(i = 0; i < mateInduced.size(); ++i){
        lengthMateInduced.push_back(mateInduced[i].size());
    }

    cout << "Mate-Induced Structure:\n";
    for(i = 0; i < mateInduced.size(); ++i){
        for(j = 0; j < mateInduced[i].size(); ++j){
            cout << mateInduced[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;

    cout << "Number of cycles in mate-induced structure: " << totalCycles << endl;

    cout << "Number of vertices in each cycle of the mate-induced structure:\n";
    for(i = 0; i < lengthMateInduced.size(); ++i){
        cout << "Cycle " << i+1 << ": " << lengthMateInduced[i] << " vertices" << endl;
    }
    cout << endl;

    //If the mate-induced structure only consists of one cycle, then the problem has been solved and is feasible (just remove one matching edge to find feasible path)
    if(lengthMateInduced[0] == numScores){ //if all of the vertices are in the first (and only) cycle of the mate-induced structure
        cout << "Instance is feasible, mate-induced structure only consists of one cycle.\n";
        cout << "Feasible order of scores:\n";
        for(i = 0; i < mateInduced[0].size()-1; ++i){
            cout << allScores[mateInduced[0][i]] << " -> ";
        }
        cout << allScores[mateInduced[0][mateInduced[0].size()-1]] << endl;
        ++feasible;
        goto End;
        //continue;
    }



    //FCA
    //create list cycleVertex that contains for each vertex the cycle that each edge belongs to
    for(i = 0; i < mateInduced.size(); ++i){
        for(j = 0; j < mateInduced[i].size(); ++j){
            cycleVertex[mateInduced[i][j]] = i;
        }
    }
    cout << "Cycle Vertex:\n";
    for(i = 0; i < cycleVertex.size(); ++i){
        cout << cycleVertex[i] << " ";
    }
    cout << endl;

    //create list of edges without empty edges (those generated by mate swap)
    currentEdge = 0;
    for(i = 0; i < matchSize; ++i){ //matchSize should = numBox, as the edges vector is of size numBox
        edge[i] = vacant;
    }

    for(i = 0; i < matchSize; ++i){
        while(cycleVertex[i] == vacant){
            ++i;
        }
        edge[currentEdge] = i;
        ++currentEdge;
    }
    numEdges = currentEdge;

    cout << "Edges vector:\n";
    for(i = 0; i < edge.size(); ++i){
        cout << edge[i] << " ";
    }
    cout << endl;
    cout << "Number of Edges: " << numEdges << endl;




































    End:
    cout << "End of Program.\n"; //for poor matching

	endTime = clock();
	int totalTime = (int)(((endTime - startTime) / double(CLOCKS_PER_SEC)) * 100);
	cout << "CPU Time = " << totalTime << " milliseconds.\n";

}