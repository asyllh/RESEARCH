/*--------------------/
ALH
MSSP with output of solution.
17/07/2017
/--------------------*/

#include <iostream>
#include <vector>
#include <time.h>
#include <algorithm>
#include <iomanip> //header providing parametric manipulators
using namespace std;

void resetVectors(int vacant, int numScores, int numComp,  vector<vector<int> > &adjMatrix, vector<int> &allScores, vector<int> &checked, vector<int> &cycleVertex, vector<int> &matchList, vector<int> &mates, vector<int> &QSet, vector<vector<int> > &S){
	int i, j;

	for(i = 0; i < numScores; ++i){
		allScores[i] = 0;
		matchList[i] = vacant;
		cycleVertex[i] = 1;
		mates[i] = 0;
		checked[i] = 0;
	}

	for(i = 0; i < numScores; ++i){
		for(j = 0; j < numScores; ++j){
			adjMatrix[i][j] = 0;
		}
	}

	for(i = 0; i < numComp; ++i){
		QSet[i] = 0;
	}

	for(i = 0; i < numComp; ++i){
		for(j = 0; j < numComp; ++j){
			S[i][j] = 0;
		}
	}

}

void clearVectors(vector<int> &cycle, vector<int> &edge, vector<int> &lengthMateInduced, vector<vector<int> > &mateInduced, vector<int> &randOrder, vector<int> &t, vector<vector<int> > &T, vector<int> &SSet, vector<int> &fullCycle){
	cycle.clear();
	edge.clear();
	lengthMateInduced.clear();
	mateInduced.clear();
	randOrder.clear();
	t.clear();
	T.clear();
	SSet.clear();
    fullCycle.clear();

}

void loopCycle(int v, int full, int save, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &fullCycle, vector<vector<int> > &mateInduced, vector<vector<int> > &T){

    int i;


    //CASE ONE: if element matchList[T[full][v]] is before element T[full][v] in the mateInduced cycle
    //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save + 1" is T[full][v]
    if(mateInduced[cycleVertex[T[full][v]]][save + 1] == T[full][v]){
        for(i = save + 1; i-- > 0;){ //from element at position 'save' to the first element in the cycle
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
        for(i = mateInduced[cycleVertex[T[full][v]]].size(); i-- > save +1;){ //from end of cycle to element at position save+1
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
    }

    //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
    //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
    else if (mateInduced[cycleVertex[T[full][v]]][save - 1] == T[full][v]){
        for(i = save; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i){
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
        for(i = 0; i < save; ++i){
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
    }

    //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
    //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
    else if (save == 0 && mateInduced[cycleVertex[T[full][v]]][mateInduced[cycleVertex[T[full][v]]].size()-1] == T[full][v]){
        for(i = 0; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i){
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
    }

    //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
    //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
    else if(save == mateInduced[cycleVertex[T[full][v]]].size()-1 && mateInduced[cycleVertex[T[full][v]]][0] == T[full][v]){
        for(i = mateInduced[cycleVertex[T[full][v]]].size(); i-- > 0; ){
            fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        }
    }

}


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

	//VARIABLES FROM ARGUMENTS
	int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
	int numBox = atoi(argv[2]) + 1; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
	int minWidth = atoi(argv[3]); //minimum width of scores (millimeters)
	int maxWidth = atoi(argv[4]); //maximum width of scores (millimeters)
	int randomSeed = atoi(argv[5]); //random seed


	//VARIABLES
	int i, j, k, q;
	int instance; //counter for instances loop
	int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
	int numComp = (numBox + (numBox % 2)) / 2;
	int threshold = 70; //adjacency threshold of scores, minimum knife distance
	int vacant = 999; //large empty value
	int feasible = 0; //number of feasible instances
	int infeasible = 0; //number of infeasible instances


	//Creating scores and forming adjacency matrix:
	vector<int> allScores(numScores, 0); //vector containing all score widths
	vector<int> randOrder; //vector used to shuffle and assign mates
	vector<vector<int> > adjMatrix(numScores, vector<int>(numScores ,0)); //adjacency matrix, 0 if width sum < threshold, 1 if width sum >= threshold, 2 if scores are mates (either side of same box)


	//Modified Threshold Graph Matching Algorithm (MTGMA)
	int lastMatch; //takes value of the index of the last vertex to have been matched
	int mateMatch; //vertex number for matching algorithm, mate takes the value of the index of the vertex that the current vertex is mates with
	int matchSize; //size/cardinality of the matching list
	int noMatch = 0; //counter for number of instances with |M| < n (matchSize < numBox)
	int vacantFlag; //flag vertices that have not been matched with the highest vertex available due to that vertex being the current vertex's mate
	int verticesNotMatched; //counts number of vertices that have not been matched to another vertex in MTGMA
	vector<int> cycleVertex(numScores, 1); //contains the number of the cycle of the mate-induced structure that the vertex i belongs to
	//(e.g. if vertex 4 is in the first cycle of the MIS, then cycleVertex[4] = 0 (0 = first cycle, 1 = second cycle etc))
	vector<int> matchList(numScores, vacant); //contains vertex index for matching vertices, e.g. if vertex 0 is matched with vertex 9, then matchList[0] = 9


	//Mate-Induced Structure (MIS)
	int smallestVertex; //smallest vertex not yet checked (i.e. not yet included in a cycle in the MIS)
	int currentVertex; //current vertex being checked for placement in a cycle in the MIS
	int numCycles; //number of cycles in the mate-induced structure, i.e mateInduced.size() (each row in the mate-induced structure matrix represents a cycle)
	int oneCycle = 0; //counter for the number of instances where the MIS consists of one one cycle (therefore immediately feasible)
	vector<int> checked(numScores, 0); //contains 0 if vector i has not yet been included in the mate-induced structure, 1 if vector i has been placed in MIS
	vector<int> cycle; //used in building the mate-induced structure
	vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4
	vector<vector<int> > mateInduced; //size numscore by noComp, i.e. number of rows = numScores, number of columns = noComp
	vector<int> lengthMateInduced; //each element holds the value corresponding to the length of the relative cycle in the mate-induced structure


	//Familty Construction Algorithm (FCA)
	int qstar; //number of T-cycle families found (i.e T.size())
	int numEdges; //number of (non-empty) edges, (i.e. edge.size())
	int noFam = 0; //counter for number of instances where no family of T-cycles has been found
	vector<int> edge; //contains the number of lower vertex of each (non-empty) edge (max size = matchSize)
	vector<int> t; //used to build T-cycles row by row
	vector<vector<int> > T; //holds set of edges in T-cycles
	vector<vector<int> > S(numComp, vector<int>(numComp, 0)); // set of indices of cycles of MIS that have an edge in T, == 1 if edge from cycle j is used in T-cycle q (T[q][j])
	//FCA also uses the vector cycleVertex from MTGMA


	//Patching Graph Connectivity
	int SSum; //number of MIS cycles already glued together
	int SqIntS; // == 0 iff Sq intersection S == emptyset
	int problemInstance = 0; //counter for number of instances with issues (SSum > numCycles)
	vector<int> QSet(numComp, vacant); // Tq-cycles already used for gluing
	vector<int> SSet; //MIS cycles already glued together

    int x;
    int save = 0;
    int full = vacant; //the row of the T matrix which has the same number of edges as the number of T cycles
    int v;
    int poorT = 0;
    int fullT = 0;
    vector<int> fullCycle; //vector to hold vertices of final cycle in order

	//*************************************************************************


	srand(randomSeed); //seed


	cout << "Minimum Score Separation Problem - Matching-Based Alternating Hamiltonicity Recognition Algorithm\n\n";

	time_t startTime, endTime; //start clock
	startTime = clock();

	for(instance = 0; instance < numInstances; ++instance) {
        x = 0;
        full = vacant;

		resetVectors(vacant, numScores, numComp, adjMatrix, allScores, checked, cycleVertex, matchList, mates, QSet, S);
		clearVectors(cycle, edge, lengthMateInduced, mateInduced, randOrder, t, T, SSet, fullCycle);

		//Create random values to be used as score widths, put in allScores vector (except last two elements)
		for (i = 0; i < numScores - 2; ++i) {
			allScores[i] = rand() % (maxWidth - minWidth + 1) + minWidth;
		}
		//add two dominating vertices with score widths = 71 (these scores will be either side of same box, mates)
		allScores[numScores - 2] = 70;
		allScores[numScores - 1] = 70;

		//Sort all of the scores in the allScores vector in ascending order
		sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

		/*cout << "All scores - non-decreasing order:\n";
		for (i = 0; i < allScores.size(); ++i) {
			cout << allScores[i] << " ";
		}
		cout << endl << endl;*/

		//Filling in adjacency matrix - if sum of two scores >= threshold (70), then insert 1 into the matrix, else leave as 0
		for (i = 0; i < allScores.size() - 1; ++i) {
			for (j = i + 1; j < allScores.size(); ++j) {
				if (allScores[i] + allScores[j] >= threshold) {
					adjMatrix[i][j] = 1;
					adjMatrix[j][i] = 1;
				}
			}

		}

		//Initially, randOrder vector will contain elements in the order 0, ..., numScores -2, numScores -1
		for (i = 0; i < numScores; ++i) {
			randOrder.push_back(i);
		}

		//Randomly shuffle all values in randOrder vector EXCEPT the last two values (dominating vertices, must stay as mates)
		random_shuffle(randOrder.begin(), randOrder.begin() + (numScores - 2));

		//Assign mates to each score (i.e. pair up scores to define which scores are either side of the same box)
		//In the adjacency matrix, this will be represented by value 2
		//Therefore there will be a value of 2 in every row and every column, non repeating
		for (i = 0; i < numBox; ++i) {
			adjMatrix[randOrder[2 * i]][randOrder[2 * i + 1]] = 2;
			adjMatrix[randOrder[2 * i + 1]][randOrder[2 * i]] = 2;
		}


		//MATCHING ALGORITHM MTGMA
		matchSize = 0;
		lastMatch = vacant;
		verticesNotMatched = 0;

		//THE MTGMA ALGORITHM
		for (i = 0; i < numScores; ++i) { //check all vertices
			vacantFlag = 0;
			if (matchList[i] == vacant) { //if vertex has not yet been matched
				for (j = numScores - 1; j > i; --j) { //try match vertex i with largest unmatched vertex, start from largest vertex j, go down list of vertices in decreasing order of size
					if (adjMatrix[i][j] == 1 && matchList[j] == vacant) { //if vertices i and j are adjacent, and if vertex j has not yet been matched
						matchList[i] = j;
						matchList[j] = i;
						lastMatch = i;
						++matchSize;
						if (vacantFlag == 1) { //delete edge for FCA if matching was not with highest vertex due to the highest vertex being its mate
							cycleVertex[i] = vacant;
							cycleVertex[j] = vacant;
						}
						break;
					} else if (adjMatrix[i][j] == 2 && matchList[j] == vacant) { //if potential match == mate
						vacantFlag = 1;
					}
				}//end for j
				if (matchList[i] == vacant) { //if vertex has still not been matched
					for (k = 0; k < numScores; ++k) {
						if (adjMatrix[i][k] == 2) { //if vertex i and vertex k are mates
							mateMatch = k;
							break;
						}
					}
					if ((allScores[i] + allScores[mateMatch] >= threshold) //match with mate?
						&& (matchList[mateMatch] == vacant) //is mate unmatched?
						&& (lastMatch != vacant) //has the previous vertex been matched?
						&& (mateMatch > i) //is the mate larger? (sorted in increasing order of vertex weight, so index will be higher if vertex has larger value)
						&& (allScores[lastMatch] + allScores[mateMatch] >= threshold)) { //can mate be matched with last matched vertex?
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

		/*cout << "Cycle Vertex vector after MTGMA:\n";
		for (i = 0; i < cycleVertex.size(); ++i) {
			cout << cycleVertex[i] << " ";
		}
		cout << endl;*/


		/*cout << "Matching List:\n";
		for (i = 0; i < numScores; ++i) {
			cout << matchList[i] << " ";
		}
		cout << endl << endl;*/


		//If the number of matches (i.e. the size of the matching list M) is less than the number of boxes (n), then instance is infeasible ( |M| < n )
		if (matchSize < numBox) {
            cout << instance << ": Not enough matching edges.\n\n";
			++infeasible;
			++noMatch;
			//cout << "Instance is infeasible, not enough matching edges available (|M| < n)." << endl;
			continue;
		}



		//MATE-INDUCED STRUCTURE
		for (i = 0; i < numScores; ++i) {
			for (j = 0; j < numScores; ++j) {
				if (adjMatrix[i][j] == 2) {
					mates[i] = j;
					break;
				}
			}
		}

		/*cout << "Mates Vector:\n";
		for (i = 0; i < numScores; ++i) {
			cout << mates[i] << " ";
		}
		cout << endl << endl;*/

		//find the smallest vertex not yet checked for mate-induced structure - start with this vertex
		for (i = 0; i < numScores; ++i) {
			if (checked[i] == 0) {
				smallestVertex = i;
				break;
			}
		}

		//Building the mate-induced structure
		do {
			currentVertex = smallestVertex;
			do {
				cycle.push_back(currentVertex);
				checked[currentVertex] = 1;
				cycle.push_back(mates[currentVertex]);
				checked[mates[currentVertex]] = 1;
				currentVertex = matchList[mates[currentVertex]];
			} while (currentVertex != smallestVertex);

			mateInduced.push_back(cycle);
			cycle.clear();

			for (i = 0; i < numScores; ++i) {
				if (checked[i] == 0) {
					smallestVertex = i;
					break;
				}
			}


		} while (smallestVertex != currentVertex);

		cycle.clear(); //clear cycle vector again for next instance

		numCycles = mateInduced.size(); //number of cycles in the mate-induced structure

		/*cout << "Mate-Induced Structure:\n";
		for (i = 0; i < mateInduced.size(); ++i) {
			for (j = 0; j < mateInduced[i].size(); ++j) {
				cout << mateInduced[i][j] << "\t";
			}
			cout << endl;
		}
		cout << endl;*/

		//cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

		for (i = 0; i < mateInduced.size(); ++i) {
			lengthMateInduced.push_back(mateInduced[i].size());
		}


		//If the mate-induced structure only consists of one cycle, then the problem has been solved and is feasible (just remove one matching edge to find feasible path)
		if (lengthMateInduced[0] == numScores) { //if all of the vertices are in the first (and only) cycle of the mate-induced structure
			++feasible;
			++oneCycle;
            cout << instance << ": Full cycle (MIS):\n";
            for(j = 0; j < mateInduced[0].size(); ++j){
                cout << mateInduced[0][j] << " ";

            }
            cout << endl << endl;
			continue;
		}



		//FCA
		//create list cycleVertex that contains for each vertex the cycle that each edge belongs to
		for (i = 0; i < mateInduced.size(); ++i) {
			for (j = 0; j < mateInduced[i].size(); ++j) {
				if (cycleVertex[mateInduced[i][j]] != vacant) { //if edge is not deleted for FCA
					cycleVertex[mateInduced[i][j]] = i;
				}
			}
		}

		/*cout << "Cycle Vertex:\n";
		for (i = 0; i < cycleVertex.size(); ++i) {
			cout << cycleVertex[i] << " ";
		}
		cout << endl;*/



		//create list of edges without empty edges (those generated by mate swap)
		for (i = 0; i < matchSize; ++i) {
			while (cycleVertex[i] == vacant) {
				++i;
			}
			edge.push_back(i);
		}
		numEdges = edge.size();

		/*cout << "Edges vector:\n";
		for (i = 0; i < edge.size(); ++i) {
			cout << edge[i] << " ";
		}
		cout << endl;*/

		//cout << "Number of Edges: " << numEdges << endl;

		//FCA Algorithm
		qstar = -1;
		k = 0; //edge from matching that is under consideration

		do {
			while (k < numEdges - 2 && (adjMatrix[edge[k]][matchList[edge[k + 1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k + 1]])) {
				++k;
			}
			if (adjMatrix[edge[k]][matchList[edge[k + 1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k + 1]]) {
				++qstar;
				t.push_back(edge[k]);
				S[qstar][cycleVertex[edge[k]]] = 1;
				while (k < numEdges - 1 && adjMatrix[edge[k]][matchList[edge[k + 1]]] == 1 && S[qstar][cycleVertex[edge[k + 1]]] == 0) { //add more edges to current T-cycle
					++k;
					t.push_back(edge[k]);
					S[qstar][cycleVertex[edge[k]]] = 1;
				}
				T.push_back(t);
				t.clear();
			} // end if
			++k;
		} while (k < numEdges - 1);

		t.clear();

		/*cout << "T matrix:\n";
		for(i = 0; i < T.size(); ++i){
			for(j = 0; j < T[i].size(); ++j){
				cout << T[i][j] << "  ";
			}
			cout << endl;
		}
		cout << endl;*/

		/*cout << "S Matrix:\n";
		for(i = 0; i < T.size(); ++i){
			for(j = 0; j < numCycles; ++j){
				cout << S[i][j] << "  ";
			}
			cout << endl;
		}
		cout << endl;
		 */

		//cout << "qstar: " << qstar << endl;



		//No family of T-cycle found
		if (qstar == -1) {
			//cout << "Instance is infeasible, no family of Tq-cycles found (q* = 0)." << endl;
            cout << instance << ": Infeasible, qstar = -1\n\n";
			++infeasible;
			++noFam;
			continue;
		}


		//CHECK IF PATCHING GRAPH IS CONNECTED
		//Setup
		q = 0; //Start with first Tq-cycle
		QSet[0] = 1;


		SSum = 0; //number of MIS-cycles that have been included
		for (i = 0; i < numCycles; ++i) {
			SSet.push_back(S[q][i]); // ==1 if MIS cycle i has been included
		}
		for (i = 0; i < numCycles; ++i) {
			SSum = SSum + SSet[i];
		}

		//Start connectivity check
		while (q <= qstar && SSum < numCycles) {
			do {
				++q;
				SqIntS = vacant;
				if (q <= qstar) {
					for (j = 0; j < numCycles; ++j) { //is there a j such that S[q][j] = 1 and SSet[j] = 1?
						if (S[q][j] == 1 && SSet[j] == 1) {
							SqIntS = 1;
							//break here? no need to check all other j indices once one has been found such that S[q][j] =1 and SSet[j] = 1
						}
					}
				}
			} while (q < qstar + 1 && (QSet[q] == 1 || SqIntS == vacant));

			if (q <= qstar) { //if Tq-cyce for enlargement has been found
				for (i = 0; i < numCycles; ++i) {
					if (SSet[i] == 0 && S[q][i] == 1) {
						SSet[i] = 1;
						++SSum;
					}
				}
				QSet[q] = 1;
				q = 0;
			}
		}//end while


		//If patching graph is connected, then instance is feasible, else infeasible
		if (SSum == numCycles) {
			//cout << "FEASIBLE: Patching Graph Connected (SSum == numCycles).\n";
			++feasible;
			goto End;
		}
		else if (SSum < numCycles) {
			//cout << "INFEASIBLE: Patching Graph Unconnected (SSum < numCycles).\n";
            cout << instance << ": Infeasible SSum < numCycles\n\n";
			++infeasible;
			continue;
		}
		else {
            cout << instance << ": Problem.\n\n";
			++problemInstance;
			continue;
		}


        End:
        for(i = 0; i < T.size(); ++i){
            if(T[i].size() == numCycles){
                full = i;
                break;
            }
        }

        if(full != vacant) {

            //cout << "Full: " << full << endl;

            for (v = 0; v < T[full].size(); ++v) {

                for (j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j) {
                    if (mateInduced[cycleVertex[T[full][v]]][j] == matchList[T[full][v]]) {
                        save = j;
                        break;
                    }
                }
                loopCycle(v, full, save, matchList, cycleVertex, fullCycle, mateInduced, T);

            }

            cout << instance << ": Full cycle after T-cycle analysis:\n";
            for (i = 0; i < fullCycle.size(); ++i) {
                cout << fullCycle[i] << " ";
            }
            cout << endl << endl;
            ++fullT;
            continue;

        }
        else if(full == vacant){
            cout << instance << ": NO SOLN AVAILABLE: multiple t cycles required.\n";
            ++poorT;
            continue;
        }


	} //end of for loop instances

    cout << "\n------------------------------------------------------------------\n\n";
	cout << "Number of feasible instances: " << feasible << endl;
	cout << "Number of infeasible instances: " << infeasible << endl;
    cout << "Number of full T cycles: " << fullT << endl;
    cout << "Number of split T cycles: " << poorT << endl;
    cout << "Number of poor matchings: " << noMatch << endl;
    cout << "Number of single cycles: " << oneCycle << endl;




	//matchlist of any of the values in T (i.e the matching mate of 16, 17 and 18) is always linked back to the edge before it, unless it's the first matching mate, in which case
	//it links to the last edge in T
	/*i.e.
	16 is the value in T, its matching mate is 25
	17 is the value in T, its matching mate is 24
	18 is the value in T, its matching mate is 23

	in this case
	16 will link with the matching mate of the next edge, so 24
	17 will link with the matching mate of the next edge, so 23
	18 is the last edge in T, so it will link with the matching mate of the first edge, so 25

	16-24
	17-23
	18-25

	*/


	endTime = clock();
	int totalTime = (int)(((endTime - startTime) / double(CLOCKS_PER_SEC)) * 100);
	cout << "CPU Time = " << totalTime << " milliseconds.\n";


}//END INT MAIN


//********************************************************************************

//Print out allScores vector
/*cout << "All scores:\n";
for (i = 0; i < allScores.size(); ++i) {
    cout << allScores[i] << " ";
}
cout << endl << endl;*/

//Print out allScores vector (scores now in ascending order)
/*if(instance == 0 || instance == 1 || instance == 2 || instance == 3) {
    cout << "All scores - non-decreasing order:\n";
    for (i = 0; i < allScores.size(); ++i) {
        cout << allScores[i] << " ";
    }
    cout << endl << endl;
}*/

//Print out randOrder vector
/* cout << "Random Order:\n";
 for(i = 0; i < randOrder.size(); ++i){
     cout << randOrder[i] << endl;
 }
 cout << endl;*/


//Print out adjacency matrix inc threshold and mates
/*cout << "Adjacency Matrix:\n";
for (i = 0; i < adjMatrix.size(); ++i) {
    for (j = 0; j < adjMatrix[i].size(); ++j) {
        cout << adjMatrix[i][j] << "\t";
    }
    cout << endl;
}
cout << endl;*/

//Fill matchingList vector with values 0,..., numScores-1 (i.e. the index of each element)
/*for(i = 0; i < numScores; ++i){
    matchList[i] = vacant;
    //cycleVertex[i] = 1; //will be set to vacant if resulting from mate swap
}*/

/*cout << "Cycle Vertex vector after MTGMA:\n";
        for (i = 0; i < cycleVertex.size(); ++i) {
            cout << cycleVertex[i] << " ";
        }
        cout << endl;*/


/*cout << "Matching List:\n";
for (i = 0; i < numScores; ++i) {
    cout << matchList[i] << " ";
}
cout << endl << endl;*/

/*cout << "Matching Vertices Values:\n";
for (i = 0; i < numScores; ++i) {
    if (matchList[i] == vacant) {
        cout << i << "\t" << allScores[i] << "\t" << "No Match" << endl;
        ++verticesNotMatched;
    } else {
        cout << i << "\t" << allScores[i] << "\t" << allScores[matchList[i]] << "\t" << matchList[i] << endl;
    }

}*/

/*if (verticesNotMatched == 0) {
    cout << "All vertices have been matched.\n\n";
    cout << "Size of M (matchSize): " << matchSize << endl;
} else {
    cout << "Number of unmatched vertices: " << verticesNotMatched << endl;
}*/

/*cout << "Mates Vector:\n";
        for (i = 0; i < numScores; ++i) {
            cout << mates[i] << " ";
        }
        cout << endl << endl;*/

/*for(i = 0; i < numScores; ++i){
    checked[i] = 0;
}*/

/*cout << "Mate-Induced Structure:\n";
       for (i = 0; i < mateInduced.size(); ++i) {
           for (j = 0; j < mateInduced[i].size(); ++j) {
               cout << mateInduced[i][j] << "\t";
           }
           cout << endl;
       }
       cout << endl;*/

//cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

/*cout << "Number of vertices in each cycle of the mate-induced structure:\n";
for (i = 0; i < lengthMateInduced.size(); ++i) {
    cout << "Cycle " << i + 1 << ": " << lengthMateInduced[i] << " vertices" << endl;
}
cout << endl;*/

//If the mate-induced structure only consists of one cycle, then the problem has been solved and is feasible (just remove one matching edge to find feasible path)
//if (lengthMateInduced[0] == numScores) { //if all of the vertices are in the first (and only) cycle of the mate-induced structure
//cout << "FEASIBLE: mate-induced structure only consists of one cycle.\n";
//cout << "Feasible order of scores:\n";
/*for (i = 0; i < mateInduced[0].size() - 1; ++i) {
    cout << allScores[mateInduced[0][i]] << " -> ";
}
cout << allScores[mateInduced[0][mateInduced[0].size() - 1]] << endl;*/
//++feasible;
//++oneCycle;
//continue;
//}

/*cout << "Cycle Vertex:\n";
        for (i = 0; i < cycleVertex.size(); ++i) {
            cout << cycleVertex[i] << " ";
        }
        cout << endl;*/

/*cout << "Edges vector:\n";
        for (i = 0; i < edge.size(); ++i) {
            cout << edge[i] << " ";
        }
        cout << endl;*/
//cout << "Number of Edges: " << numEdges << endl;


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


/*for(q = 1; q <= qstar; ++q){
    QSet[q] = 0; // ==1 iff Tq-cycle number q has already been considered
}*/

//********************************************************************************************************************************************

/*for(j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j){
		if(mateInduced[cycleVertex[T[full][v]]][j] != T[full][v] && mateInduced[cycleVertex[T[full][v]]][j] != matchList[T[full][v]]){
			fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][j]);
		}
		else if (mateInduced[cycleVertex[T[full][v]]][j] == T[full][v]){
			brokenVertex1 = j;
			break;
		}
		else if (mateInduced[cycleVertex[T[full][v]]][j] == matchList[T[full][v]]){
			brokenVertex2 = j;
			break;
		}
	}*/
//visited[cycleVertex[T[full][v]]] = 1;


/*cout << "full cycle\n";
for(i = 0; i < fullCycle.size(); ++i){
    cout << fullCycle[i] << endl;
}*/

/*cout << "broken1: " << brokenVertex1 << endl;
cout << "broken2: " << brokenVertex2 << endl;
if(brokenVertex1 != vacant){
    fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][brokenVertex1]);
     ++v;
    fullCycle.push_back(matchList[T[full][v]]);
    brokenVertex1 = vacant;
}*/


/*if(brokenVertex2 != vacant){
    fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][brokenVertex2]);
    v = T[full].size()-1;
    cout << v << endl;
    fullCycle.push_back(T[full][v]);
    brokenVertex2 = vacant;
}*/

/*cout << "full cycle\n";
for(i = 0; i < fullCycle.size(); ++i){
    cout << fullCycle[i] << endl;
}*/

// v = 2
/*for(j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j){
    if(mateInduced[cycleVertex[T[full][v]]][j] == T[full][v]){
        save = j; //save = 3
        break;
    }
}*/

//cout << "Save: "<<  save << endl;

//THIS NEEDS TO GO IN A FUNCTION ALONG WITH THE OTHER FIVE POSSIBILITIES
//CASE 5: if the next element in the cycle is the matching mate for the current element
/*if(mateInduced[cycleVertex[T[full][v]]][save+1] == matchList[T[full][v]]){
    for(i = save; i-- >0;){
        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
        //cout << mateInduced[cycleVertex[T[full][v]]][i] << endl;
    }
    for(i = mateInduced[cycleVertex[T[full][v]]].size(); i-- > save + 1;){
        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
    }

}
visited[cycleVertex[T[full][v]]] = 1;*/

/*i = 0;
save = 0;
--v; // v = 1

for(j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j) {
    if (mateInduced[cycleVertex[T[full][v]]][j] == T[full][v]) {
        save = j; //save = 18
        break;
    }
}*/

//cout << "Save: " << save << endl;

//CASE 6: if the previous element in the cycle is the matching mate for the current element
/*if(mateInduced[cycleVertex[T[full][v]]][save-1] == matchList[T[full][v]]){
    for(i = save; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i){
        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
    }
    for(i = 0; i < save; ++i){
        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
    }

}
visited[cycleVertex[T[full][v]]] = 1;*/

/*i = 0;
save = 0;
--v; //v = 0*/

/*cout << "Visited:\n";
for(i = 0; i < visited.size(); ++i){
    cout << visited[i] << endl;
}
cout << endl;*/

//CONNECTING BACK TO FIRST CYCLE
/*for(j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j) {
    if (mateInduced[cycleVertex[T[full][v]]][j] == T[full][v]) {
        save = j; //save = 4
        break;
    }
}*/

//CYCLE 0 HAS ALREADY BEEN VISITED, CREATE NEW FUNCTION FOR THIS TYPE
/*if(mateInduced[cycleVertex[T[full][v]]][save-1] == matchList[T[full][v]] && visited[cycleVertex[T[full][v]]] == 1){
    for(i = save; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i){
        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
    }
}*/



/*cout << "full cycle\n";
for(i = 0; i < fullCycle.size(); ++i){
    cout << fullCycle[i] << " ";
}
cout << endl;*/
