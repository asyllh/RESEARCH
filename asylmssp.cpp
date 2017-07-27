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




int main(int argc, char **argv){
    if(argc < 5){
        cout << "Minimum Score Separation Problem: MBAHRA.\n";
        cout << "Arguments are the following:\n";
        cout << "- Number of instances (integer)\n";
        cout << "- Number of boxes (integer)\n";
        cout << "- Minimum width of scores (millimeters, min = 1)\n";
        cout << "- Maximum width of scores (millimeters, max = 70)\n";
        cout << "- Random Seed (integer)\n";
        exit(1);
    }

    //Variables from arguments
    int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
    int numBox = atoi(argv[2]) + 1; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
    int minWidth = atoi(argv[3]); //minimum width of scores (millimeters)
    int maxWidth = atoi(argv[4]); //maximum width of scores (millimeters)
    int randomSeed = atoi(argv[5]); // random seed


    //Variables
    int i, j, k, q, qstar;
    int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
    int numComp = (numBox + (numBox % 2)) / 2;
    int threshold = 70; //adjacency threshold of scores, minimum knife distance
    int mateMatch; //vertex number for matching algorithm, mate takes the value of the index of the vertex that the current vertex is mates with
    int lastMatch;
    int matchSize; //size/cardinality of the matching list
    int vacant = 999;
    int verticesNotMatched; //counts number of vertices that have not been matched to another vertex in MTGMA
    int smallestVertex;
    int currentVertex;
    int numCycles; //number of cycles in the mate-induced structure, i.e mateInduced.size() (each row in the mate-induced structure matrix represents a cycle)
    int feasible = 0; //number of feasible instances
    int infeasible = 0; //number of infeasible instances
    int currentEdge; //counter used for list of edges
    int numEdges; //number of (non-empty) edges
    int vacantFlag;
    int SSum; //number of MIS cycles already glued together
    int SqIntS; // == 0 iff Sq intersection S == emptyset


    vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4
    vector<int> matchList(numScores, vacant); //contains vertex index for matching vertices, e.g. if vertex 0 is matched with vertex 9, then matchList[0] = 9
    vector<int> allScores(numScores, 0); //vector containing all score widths
    vector<int> randOrder; //vector used to shuffle and assign mates
    vector<vector<int> > adjMatrix(numScores, vector<int>(numScores ,0)); //adjacency matrix, 0 if width sum < threshold, 1 if width sum >= threshold, 2 if scores are mates (either side of same box)
    vector<int> checked(numScores, 0); //contains 0 if vector i has not yet been included in the mate-induced structure, 1 if vector i has been placed in MIS
    vector<vector<int> > mateInduced; //size numscore by noComp, i.e. number of rows = numScores, number of columns = noComp
    vector<int> cycle; //used in building the mate-induced structure
    vector<int> lengthMateInduced; //each element holds the value corresponding to the length of the relative cycle in the mate-induced structure
    vector<int> cycleVertex(numScores, 1); //contains the number of the cycle of the mate-induced structure that the vertex i belongs to
    // (e.g. if vertex 4 is in the first cycle of the MIS, then cycleVertex[4] = 0 (0 = first cycle, 1 = second cycle etc))
    vector<int> edge; //contains the number of lower vertex of each (non-empty) edge
    vector<vector<int> > T; //T-cycles
    vector<vector<int> > S(numComp, vector<int>(numComp, 0)); // == 1 if edge from cycle is used in T-cycle
    vector<int> t;
    vector<int> s;
    vector<int> QSet(numComp, 0); // Tq-cycles already used for gluing
    vector<int> SSet; //MIS cycles already glued together

    srand(randomSeed); //seed

    //DON'T FORGET TO CLEAR VECTORS AND VARIABLES FOR NEXT INSTANCE
    //FUNCTION TO RESET ALL VECTORS

    cout << "Minimum Score Separation Problem\nMatching-Based Alternating Hamiltonicity Recognition Algorithm\n\n";

    time_t startTime, endTime; //start clock
    startTime = clock();

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for(i = 0; i < numScores - 2; ++i){
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


    //Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
    for(i = 0; i < numScores; ++i){
        randOrder.push_back(i);
    }

    //Randomly shuffle all values in randOrder vector EXCEPT the last two values (dominating vertices, must stay as mates)
    random_shuffle(randOrder.begin(), randOrder.begin()+(numScores-2));

    //Print out randOrder vector
   /* cout << "Random Order:\n";
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
    /*for(i = 0; i < numScores; ++i){
        matchList[i] = vacant;
        //cycleVertex[i] = 1; //will be set to vacant if resulting from mate swap
    }*/
    matchSize = 0;
    lastMatch = vacant;
    verticesNotMatched = 0;

    //THE MTGMA ALGORITHM
    for(i = 0; i < numScores; ++i){ //check all vertices
        vacantFlag = 0;
        if(matchList[i] == vacant){ //if vertex has not yet been matched
            for(j = numScores -1; j > i; --j){ //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
                if(adjMatrix[i][j] == 1 && matchList[j] == vacant){ //if vertices i and j are adjacent, and if vertex j has not yet been matched
                    matchList[i] = j;
                    matchList[j] = i;
                    lastMatch = i;
                    ++matchSize;
                    if(vacantFlag == 1){ //delete edge for FCA if matching was not with highest vertex due to the highest vertex being its mate
                        cycleVertex[i] = vacant;
                        cycleVertex[j] = vacant;
                    }
                    break;
                }
                else if(adjMatrix[i][j] == 2 && matchList[j] == vacant){ //if potential match == mate
                    vacantFlag = 1;
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
                    cycleVertex[lastMatch] = vacant; //edge from mate swap will not count for FCA
                    cycleVertex[mateMatch] = vacant; //edge from mate swap will not count for FCA
                    lastMatch = i;
                    ++matchSize;
                }
            }//end if matchList == vacant
        }//end if matchList[i] == i
    }//end for i

    cout << "Cycle Vertex vector after MTGMA:\n";
    for(i = 0; i < cycleVertex.size(); ++i){
        cout << cycleVertex[i] << " ";
    }
    cout << endl;


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
        cout << "Size of M (matchSize): " << matchSize << endl;
    }
    else {
        cout << "Number of unmatched vertices: " << verticesNotMatched << endl;
    }

    //If the number of matches (i.e. the size of the matching list M) is less than the number of boxes (n), then instance is infeasible ( |M| < n )
    if(matchSize < numBox){
        ++infeasible;
        cout << "Instance is infeasible, not enough matching edges available (|M| < n)." << endl;
        goto End;
        //continue;
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

    /*for(i = 0; i < numScores; ++i){
        checked[i] = 0;
    }*/

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

    numCycles = mateInduced.size(); //number of cycles in the mate-induced structure

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

    cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

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
            if(cycleVertex[mateInduced[i][j]] != vacant){ //if edge is not deleted for FCA
                cycleVertex[mateInduced[i][j]] = i;
            }
        }
    }
    cout << "Cycle Vertex:\n";
    for(i = 0; i < cycleVertex.size(); ++i){
        cout << cycleVertex[i] << " ";
    }
    cout << endl;

    //create list of edges without empty edges (those generated by mate swap)
    currentEdge = 0;

    for(i = 0; i < matchSize; ++i){
        while(cycleVertex[i] == vacant){
            ++i;
        }
        edge.push_back(i);
    }
    numEdges = edge.size();

    cout << "Edges vector:\n";
    for(i = 0; i < edge.size(); ++i){
        cout << edge[i] << " ";
    }
    cout << endl;
    cout << "Number of Edges: " << numEdges << endl;

    //FCA Algorithm
    qstar = -1;
    k = 0; //edge from matching that is under consideration

    do{
        while(k < numEdges - 2 && (adjMatrix[edge[k]][matchList[edge[k+1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k+1]])){
            ++k;
        }
        if(adjMatrix[edge[k]][matchList[edge[k+1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k+1]]){
            ++qstar;
            t.push_back(k);
            S[qstar][cycleVertex[edge[k]]] = 1;
            while(k < numEdges - 1 && adjMatrix[edge[k]][matchList[edge[k+1]]] == 1 && S[qstar][cycleVertex[edge[k+1]]] == 0){ //add more edges to current T-cycle
                ++k;
                t.push_back(k);
                S[qstar][cycleVertex[edge[k]]] = 1;
            }
            T.push_back(t);
            t.clear();
        } // end if
        ++k;
    } while (k < numEdges -1);

    t.clear();

    /*cout << "T matrix:\n";
    for(i = 0; i < T.size(); ++i){
        for(j = 0; j < T[i].size(); ++j){
            cout << T[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl << endl;

    cout << "S Matrix:\n";
    for(i = 0; i < S.size(); ++i){
        for(j = 0; j < S[i].size(); ++j){
            cout << S[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;*/

    //cout << "qstar: " << qstar << endl;

    //No family of T-cycle found
    if(qstar == -1){
        cout << "Instance is infeasible, no family of Tq-cycles found (q* = 0)." << endl;
        ++infeasible;
        goto End;
        //continue;
    }


    //CHECK IF PATCHING GRAPH IS CONNECTED
    //Setup

    /*for(q = 1; q <= qstar; ++q){
        QSet[q] = 0; // ==1 iff Tq-cycle number q has already been considered
    }*/

    q = 0; //Start with first Tq-cycle
    QSet[0] = 1;

    SSum = 0; //number of MIS-cycles that have been included
    for(i = 0; i < numCycles; ++i){
        SSet.push_back(S[q][i]); // ==1 if MIS cycle i has been included
    }
    for(i = 0; i < numCycles; ++i){
        SSum = SSum + SSet[i];
    }

    //Start connectivity check
    while(q <= qstar && SSum < numCycles){
        do{
            ++q;
            SqIntS = vacant;
            if(q <= qstar){
                for(j = 0; j < numCycles; ++j){ //is there a j such that S[q][j] = 1 and SSet[j] = 1?
                    if(S[q][j] == 1 && SSet[j] == 1){
                        SqIntS = 1;
                        //break here? no need to check all other j indices once one has been found such that S[q][j] =1 and SSet[j] = 1
                    }
                }
            }
        } while (q < qstar + 1 && (QSet[q] == 1 || SqIntS == vacant));

        if(q <= qstar){ //if Tq-cyce for enlargement has been found
            for(i = 0; i < numCycles; ++i){
                if(SSet[i] == 0 && S[q][i] == 1){
                    SSet[i] = 1;
                    ++SSum;
                }
            }
            QSet[q] = 1;
            q = 0;
        }
    }//end while


    //If patching graph is connected, then instance is feasible, else infeasible
    if(SSum == numCycles){
        cout << "FEASIBLE: Patching Graph Connected (SSum == numCycles).\n";
        ++feasible;
    }
    else if (SSum < numCycles){
        cout << "INFEASIBLE: Patching Graph Unconnected (SSum < numCycles).\n";
        ++infeasible;
    }
    else{
        cout << "PROBLEM: SSum > numCycles.\n";
        goto End;
        //continue
    }





    End:
    cout << "End of Program.\n";

	endTime = clock();
	int totalTime = (int)(((endTime - startTime) / double(CLOCKS_PER_SEC)) * 100);
	cout << "CPU Time = " << totalTime << " milliseconds.\n";

}



