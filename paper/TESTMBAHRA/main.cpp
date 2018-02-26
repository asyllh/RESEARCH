/*--------------
ALH
Fixing bug in mbahra from SCSPP
16/02/2018
---------------*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <random>
#include <cmath>

using namespace std;

int main(int argc, char **argv) {

	//int randseed = atoi(argv[1]);
	int instance;
    int i, j, k;
    int feasible = 0;
    int threshold = 70;
    int vacant = 999;
    //vector<int> scores = {1, 69, 8, 3, 70, 13, 68, 10, 68, 6, 70, 30, 62, 45, 57, 5, 19, 43, 70, 70};
    vector<int> scores;
    int minWidth = 1;
    int maxWidth = atoi(argv[1]);
    int nScores = atoi(argv[2]);
    int nBox = nScores / 2;
    int nComp = (nBox + (nBox%2)) /2;
    vector<int> randOrder;
    vector<int> order;
    vector<int> final;
    //vector<int> original = {4, 974, 118, 32, 989, 182, 951, 133, 968, 80, 991, 439, 876, 633, 794, 69, 279, 611};
    //int nScores = scores.size();
    //int nBox = nScores / 2;
    //int nComp = (nBox + (nBox % 2)) / 2;
    vector<int> invOrder(nScores);
    vector<vector<int> > adjMat(nScores, vector<int>(nScores, 0));
    vector<int> mates(nScores, vacant);
    vector<int> fullCycle;
    vector<int> completePath;
    //srand(randseed);

for(instance = 10; instance < 100; ++instance){
	//srand(1);
	scores.clear();
	randOrder.clear();
	for(i = 0; i < adjMat.size(); ++i){
		for(j = 0; j < adjMat[i].size(); ++j){
			adjMat[i][j] = 0;
		}
	}
	for(i = 0; i < mates.size(); ++i){
		mates[i] = vacant;
	}
	completePath.clear();
	fullCycle.clear();
	final.clear();
	feasible = 0;
    if (feasible == 0) {

        /*for (i = 0; i < nScores; ++i) {
            order.push_back(i);
        }

        for (i = 1; i < nScores; ++i) {
            for (j = i - 1; j >= 0; --j) {
                if (scores[i] < scores[order[j]]) {
                    order[j + 1] = order[j];
                    order[j] = i;
                }
            }
        }*/

        /*cout << "Order:\n";
        for (int v : order) {
            cout << v << " ";
        }
        cout << endl;

        for (i = 0; i < nScores; ++i) {
            invOrder[order[i]] = i;
        }

        cout << "invOrder:\n";
        for (int v : invOrder) {
            cout << v << " ";
        }
        cout << endl;*/

        /*for (i = 0; i < nScores; i += 2) {
            adjMat[invOrder[i]][invOrder[i + 1]] = 2;
            adjMat[invOrder[i + 1]][invOrder[i]] = 2;
        }*/

        //sort(scores.begin(), scores.end());

        /*for (i = 0; i < scores.size() - 1; ++i) {
            for (j = i + 1; j < scores.size(); ++j) {
                if (scores[i] + scores[j] >= threshold && adjMat[i][j] != 2) {
                    adjMat[i][j] = 1;
                    adjMat[j][i] = 1;
                }
            }
        }*/


        /*for (i = 0; i < nScores; ++i) {
            if (mates[i] == vacant) {
                for (j = 0; j < nScores; ++j) {
                    if (adjMat[i][j] == 2) {
                        mates[i] = j;
                        mates[j] = i;
                        break;
                    }
                }
            }
        }*/

        for(i = 0; i < nScores - 2; ++i){
    		scores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
	    }
	    scores.push_back(70);
	    scores.push_back(70);

	    if(instance == 22){
	    	cout << "scores:\n";
		    for(int v : scores){
		    	cout << v << " ";
		    }
		    cout << endl;
		}

	    sort(scores.begin(), scores.end());

	    for(i = 0; i < scores.size() - 1; ++i){
	    	for(j = i+1; j < scores.size(); ++j){
	    		if(scores[i] + scores[j] >= threshold){
	    			adjMat[i][j] = 1;
	    			adjMat[j][i] = 1;
	    		}
	    	}
	    }

	    for(i = 0; i < nScores - 2; ++i){
	    	randOrder.push_back(i);
	    }

	    random_shuffle(randOrder.begin(), randOrder.end());

	    /*cout << "randorder:\n";
	    for(int v : randOrder){
	    	cout << v << " ";
	    }
	    cout << endl;*/

	    for(i = 0; i < (nScores - 2) / 2; ++i){
	    	//cout << "hi\n";
	    	adjMat[randOrder[2*i]][randOrder[2*i+1]] = 2;
	    	adjMat[randOrder[2*i+1]][randOrder[2*i]] = 2;
	    }
	    //cout << "hi\n";
	    adjMat[nScores-1][nScores-2] = 2;
	    adjMat[nScores - 2][nScores - 1] = 2;

		/*cout << "adjMat:\n";
	    for(i = 0; i < adjMat.size(); ++i){
	       	for(j = 0; j < adjMat[i].size(); ++j){
	       		cout << adjMat[i][j] << " ";
	       	}
	       	cout << endl;
	    }
	    cout << endl;*/

	    for(i = 0; i < nScores; ++i){
	    	for(j = 0; j < nScores; ++j){
	    		if(adjMat[i][j] == 2){
	    			mates[i] = j;
	    			break;
	    		}
	    	}
	    }

	    if(instance == 90){
		    cout << "Mates:\n";
		    for (int v : mates) {
		        cout << v << " ";
		    }
		    cout << endl;
		}



        //MTGMA
        int vacantFlag = 0;
        int matchSize = 0;
        int lastMatch = vacant;
        int mateMatch = vacant;
        vector<int> cycleVertex(nScores, 1);
        vector<int> matchList(nScores, vacant);

        for (i = 0; i < nScores; ++i) {
            vacantFlag = 0;
            if (matchList[i] == vacant) {
                for (j = nScores - 1; j > i; --j) {
                    if (adjMat[i][j] == 1 && matchList[j] == vacant) {
                        matchList[i] = j;
                        matchList[j] = i;
                        lastMatch = i;
                        ++matchSize;
                        if (vacantFlag == 1) {
                            cycleVertex[i] = vacant;
                            cycleVertex[j] = vacant;
                        }
                        break;
                    }
                    else if (adjMat[i][j] == 2 && matchList[j] == vacant) {
                        vacantFlag = 1;
                    }
                }
                if (matchList[i] == vacant) {
                    mateMatch = mates[i]; //check if matches line 588 in SCSPP.
                    if ((scores[i] + scores[mateMatch] >= threshold) && (matchList[mateMatch] == vacant) &&
                        (lastMatch != vacant)
                        && (mateMatch > i) && (scores[lastMatch] + scores[mateMatch] >= threshold)) {
                        matchList[i] = matchList[lastMatch];
                        matchList[lastMatch] = mateMatch;
                        matchList[mateMatch] = lastMatch;
                        matchList[matchList[i]] = i;
                        cycleVertex[lastMatch] = vacant;
                        cycleVertex[mateMatch] = vacant;
                        lastMatch = i;
                        ++matchSize;
                    }
                }
            }
            //cout << "hi2\n";
        }
        //cout << matchSize << endl;

        if(instance == 90){
        	cout << "matchlist:\n";
        	for(int c: matchList){
        		cout << c << " ";
        	}
        	cout << endl;
     	}


        if (matchSize < nBox) {
            feasible = 0;
            goto End;
        }

        /**MIS**/
        int numCycles = 0;
        int smallestVertex;
        int currentVertex;
        vector<int> lengthMateInduced;
        vector<int> tempMIS;
        vector<int> checked(nScores, 0);
        vector<vector<int> > mateInduced;

        for (i = 0; i < nScores; ++i) {
            if (checked[i] == 0) {
                smallestVertex = i;
                break;
            }
        }

        do {
            currentVertex = smallestVertex;
            do {
                tempMIS.push_back(currentVertex);
                checked[currentVertex] = 1;
                tempMIS.push_back(mates[currentVertex]);
                checked[mates[currentVertex]] = 1;
                currentVertex = matchList[mates[currentVertex]];
            } while (currentVertex != smallestVertex);
            mateInduced.push_back(tempMIS);
            tempMIS.clear();
            for (i = 0; i < nScores; ++i) {
                if (checked[i] == 0) {
                    smallestVertex = i;
                    break;
                }
            }
        } while (smallestVertex != currentVertex);
        tempMIS.clear();
        numCycles = mateInduced.size();

        for (i = 0; i < mateInduced.size(); ++i) {
            lengthMateInduced.push_back(mateInduced[i].size());
        }

        if(mateInduced.size() > 4){
	        cout << "MPS - instance = " << instance << ": \n";
	        for(i = 0; i < mateInduced.size(); ++i){
	        	for(j = 0; j < mateInduced[i].size(); ++j){
	        		cout << mateInduced[i][j] << " ";
	        	}
	        	cout << endl;
	        }
	        cout << endl << endl;
    	}
    	else{
    		continue;
    	}
        //cout << "hi\n";
        if (lengthMateInduced[0] == nScores) {
            for (j = 0; j < mateInduced[0].size(); ++j) {
                fullCycle.push_back(mateInduced[0][j]);
            }

            /**MakePath**/
            for (i = 0; i < fullCycle.size() - 1; ++i) {
                if ((fullCycle[i] == nScores - 1 && fullCycle[i + 1] == nScores - 2) ||
                    (fullCycle[i] == nScores - 2 && fullCycle[i + 1] == nScores - 1)) {
                    if (i == 0) {
                        for (j = 2; j < fullCycle.size(); ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                    else if (i == fullCycle.size() - 2) {
                        for (j = 0; j < fullCycle.size() - 2; ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                    else {
                        for (j = i + 2; j < fullCycle.size(); ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                }
            }

            /*for (i = 0; i < completePath.size(); ++i) {
                final.push_back(original[order[completePath[i]]]);
            }*/
            feasible = 1;
            goto End;

        }

        /**FCA**/
        int qstar, numEdges;
        vector<int> edge;
        vector<int> t;
        vector<vector<int> > T;
        vector<vector<int> > S(nComp, vector<int>(nComp, 0));

        for (i = 0; i < mateInduced.size(); ++i) {
            for (j = 0; j < mateInduced[i].size(); ++j) {
                if (cycleVertex[mateInduced[i][j]] != vacant) {
                    cycleVertex[mateInduced[i][j]] = i;
                }
            }
        }

        for (i = 0; i < matchSize; ++i) {
            while (cycleVertex[i] == vacant) {
                ++i;
            }
            edge.push_back(i);
        }
        numEdges = edge.size();

        qstar = -1;
        k = 0;

        do {
            while (k < numEdges - 2 &&
                   (adjMat[edge[k]][matchList[edge[k + 1]]] != 1 || cycleVertex[edge[k]] == cycleVertex[edge[k + 1]])) {
                ++k;
            }
            if (adjMat[edge[k]][matchList[edge[k + 1]]] == 1 && cycleVertex[edge[k]] != cycleVertex[edge[k + 1]]) {
                ++qstar;
                t.push_back(edge[k]);
                S[qstar][cycleVertex[edge[k]]] = 1;
                while (k < numEdges - 1 && adjMat[edge[k]][matchList[edge[k + 1]]] == 1 &&
                       S[qstar][cycleVertex[edge[k + 1]]] == 0) {
                    ++k;
                    t.push_back(edge[k]);
                    S[qstar][cycleVertex[edge[k]]] = 1;
                }
                T.push_back(t);
                t.clear();
            }
            ++k;
        } while (k < numEdges - 1);
        t.clear();

        if(instance == 90){
	        cout << "T - instance = " << instance << ": \n";
	        for(i = 0; i < T.size(); ++i){
	        	for(j = 0; j < T[i].size(); ++j){
	        		cout << T[i][j] << " ";
	        	}
	        	cout << endl;
	        }
	        cout << endl << endl;
    	}

    	if(instance == 90){
	        cout << "S - instance = " << instance << ": \n";
	        for(i = 0; i < T.size(); ++i){
	        	for(j = 0; j < numCycles; ++j){
	        		cout << S[i][j] << " ";
	        	}
	        	cout << endl;
	        }
	        cout << endl << endl;
    	}


        if (qstar == -1) {
            feasible = 0;
            goto End;
        }

        /**PatchGraph**/
        int q, u, v, save, SSum, SqIntS;
        int a;
        bool found = false;
        int v1 = 0;
        int v2 = 0;
        int full = vacant;
        vector<int> SSet;
        vector<int> tempPG;
        vector<int> QSet(nComp, 0);
        vector<int> patchCycle(nComp, vacant);
        vector<int> patchVertex(nScores, vacant);
        vector<vector<int> > Tpatch;


        for (i = 0; i < T.size(); ++i) {
            if (T.size() == numCycles) {
                full = i;
                break;
            }
        }

        if (full != vacant) {
            save = 0;
            //cout << "Full: " << full << endl;
            for (v = 0; v < T[full].size(); ++v) {
                for (j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j) {
                    if (mateInduced[cycleVertex[T[full][v]]][j] == matchList[T[full][v]]) {
                        save = j;
                        break;
                    }
                }
                //region oneTCyclePatch
                /****************************oneTCyclePatch algorithm: ***********************************/
                //CASE ONE: if element matchList[T[full][v]] is before element T[full][v] in the mateInduced cycle
                //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save + 1" is T[full][v]
                if (mateInduced[cycleVertex[T[full][v]]][save + 1] == T[full][v]) {
                    for (i = save + 1; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                    for (i = mateInduced[cycleVertex[T[full][v]]].size();
                         i-- > save + 1;) { //from end of cycle to element at position save+1
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                }

                    //CASE TWO: if element matchList[T[full][v]] is after element T[full][v] in the mateInduced cycle
                    //i.e. if the element at position 'save' in the cycle is matchList[T[full][v]] and the element at position "save - 1" is T[full][v]
                else if (mateInduced[cycleVertex[T[full][v]]][save - 1] == T[full][v]) {
                    for (i = save; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i) {
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                    for (i = 0; i < save; ++i) {
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                }

                    //CASE THREE: if element matchList[T[full][v]] is the first element in the cycle, and T[full][v] is the last element in the cycle
                    //i.e. if save = 0 and T[full][v] is at position mateInduced[cycleVertex[T[full][v]]].size()-1
                else if (save == 0 &&
                         mateInduced[cycleVertex[T[full][v]]][mateInduced[cycleVertex[T[full][v]]].size() -
                                                              1] == T[full][v]) {
                    for (i = 0; i < mateInduced[cycleVertex[T[full][v]]].size(); ++i) {
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                }

                    //CASE FOUR: if element matchList[T[full][v]] is the last element in the cycle, and T[full][v] is the first element in the cycle
                    //i.e. if save = mateInduced[cycleVertex[T[full][v]]].size()-1 and T[full][v] is at position 0
                else if (save == mateInduced[cycleVertex[T[full][v]]].size() - 1 &&
                         mateInduced[cycleVertex[T[full][v]]][0] == T[full][v]) {
                    for (i = mateInduced[cycleVertex[T[full][v]]].size(); i-- > 0;) {
                        fullCycle.push_back(mateInduced[cycleVertex[T[full][v]]][i]);
                    }
                }
                //END ONETCYCLEPATCH FUNCTION
                //endregion
            }

            //region MakePath
            /**MakePath**/
            for (i = 0; i < fullCycle.size() - 1; ++i) {
                if ((fullCycle[i] == nScores - 1 && fullCycle[i + 1] == nScores - 2) ||
                    (fullCycle[i] == nScores - 2 && fullCycle[i + 1] == nScores - 1)) {
                    if (i == 0) {
                        for (j = 2; j < fullCycle.size(); ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                    else if (i == fullCycle.size() - 2) {
                        for (j = 0; j < fullCycle.size() - 2; ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                    else {
                        for (j = i + 2; j < fullCycle.size(); ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        for (j = 0; j < i; ++j) {
                            completePath.push_back(fullCycle[j]);
                        }
                        break;
                    }
                }
            }


            /*for (i = 0; i < completePath.size(); ++i) {
                final.push_back(original[order[completePath[i]]]);
            }*/
            feasible = 1;
            goto End;
            //endregion


        }

        else {

            for(a = 0; a < T.size() - 1; ++a){
                for(q = a+1; q < T.size(); ++q){
                    SqIntS = 0;
                    for(i = 0; i < numCycles; ++i){
                        if(S[a][i] + S[q][i] == 0){
                            SqIntS = 0;
                            break;
                        }
                        if(S[a][i] + S[q][i] == 2){
                            ++SqIntS;
                        }
                    }
                    if(SqIntS == 1){
                        v1 = a;
                        v2 = q;
                        SSum = numCycles;
                        found = true;
                        break;
                    }
                }
                if(found){
                    break;
                }
            }
            /*q = 0;
            QSet[0] = 1;
            SSum = 0;
            for (i = 0; i < numCycles; ++i) {
                SSet.push_back(S[q][i]);
            }
            for (i = 0; i < numCycles; ++i) {
                SSum = SSum + SSet[i];
            }
            if (SSum >= 1) {
                patchCycle[q] = 1;
            }

            while (q <= qstar && SSum < numCycles) {
                do {
                    ++q;
                    SqIntS = vacant;
                    if (q <= qstar) {
                        for (j = 0; j < numCycles; ++j) {
                            if (S[q][j] == 1 && SSet[j] == 1) {
                                SqIntS = 1;
                            }
                        }
                    }
                } while (q < qstar + 1 && (QSet[q] == 1 || SqIntS == vacant));

                if (q <= qstar) {
                    for (i = 0; i < numCycles; ++i) {
                        if (SSet[i] == 0 && S[q][i] == 1) {
                            SSet[i] = 1;
                            ++SSum;
                            patchCycle[q] = 1;
                        }
                    }
                    QSet[q] = 1;
                    q = 0;
                }
            }//endwhile*/

            //If patching graph is connected, then instance is feasible, else infeasible
            if (SSum == numCycles) {
                /*for (i = 0; i < patchCycle.size(); ++i) {
                    if (patchCycle[i] == 1) {
                        for (j = 0; j < T[i].size(); ++j) {
                            tempPG.push_back(T[i][j]);
                        }
                        Tpatch.push_back(tempPG);
                        tempPG.clear();
                    }
                }
                tempPG.clear();*/

                for(j = 0; j < T[v1].size(); ++j){
                    tempPG.push_back(T[v1][j]);
                }
                Tpatch.push_back(tempPG);
                tempPG.clear();

                for(j = 0; j < T[v2].size(); ++j){
                    tempPG.push_back(T[v2][j]);
                }
                Tpatch.push_back(tempPG);
                tempPG.clear();

                if(instance == 90){
                    cout << "Tpatch:\n";
                    for(i =0; i < Tpatch.size(); ++i){
                        for(j = 0; j < Tpatch[i].size(); ++j){
                            cout << Tpatch[i][j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                }

                vector<int> patchML;
                vector<int> inCycle(nScores, 0);

                copy(matchList.begin(), matchList.end(), back_inserter(patchML));

                /*cout << "Patch MatchList:\n";
                for(int w : patchML){
                    cout << w << " ";
                }
                cout << endl;*/

                for(u = 0; u < Tpatch.size(); ++u){
                    for(v = 0; v < Tpatch[u].size() - 1; ++v){
                        patchML[Tpatch[u][v]] = matchList[Tpatch[u][v+1]];
                        patchML[matchList[Tpatch[u][v+1]]] = Tpatch[u][v];
                    }
                    patchML[Tpatch[u][Tpatch[u].size()-1]] = matchList[Tpatch[u][0]];
                    patchML[matchList[Tpatch[u][0]]] = Tpatch[u][Tpatch[u].size()-1];
                }

                /*cout << "Patch MatchList:\n";
                for(int w : patchML){
                    cout << w << " ";
                }
                cout << endl;*/

                int current = nScores - 2;

                do{
                    fullCycle.push_back(current);
                    inCycle[current] = 1;
                    fullCycle.push_back(mates[current]);
                    inCycle[mates[current]] = 1;
                    if(inCycle[patchML[mates[current]]] == 0){
                        current = patchML[mates[current]];
                    }
                    else{
                        current = matchList[mates[current]];
                    }
                } while(fullCycle.size() < nScores);


                //region MAKEPATH
                /**MAKEPATH FUNCTION**/
                for (i = 0; i < fullCycle.size() - 1; ++i) {
                    if ((fullCycle[i] == nScores - 1 && fullCycle[i + 1] == nScores - 2) ||
                        (fullCycle[i] == nScores - 2 && fullCycle[i + 1] == nScores - 1)) {
                        if (i == 0) { //if the dominating vertices are at the beginning of the fullCycle vector
                            for (j = 2; j < fullCycle.size(); ++j) {
                                completePath.push_back(fullCycle[j]);
                            }
                            break;
                        }

                        else if (i == fullCycle.size() -
                                      2) { //if the dominating vertices are at the end of the fullCycle vector
                            for (j = 0; j < fullCycle.size() - 2; ++j) {
                                completePath.push_back(fullCycle[j]);
                            }
                            break;
                        }
                        else { //if the dominating vertices are in the middle of the fullCycle vector
                            for (j = i + 2; j < fullCycle.size(); ++j) {
                                completePath.push_back(fullCycle[j]);
                            }
                            for (j = 0; j < i; ++j) {
                                completePath.push_back(fullCycle[j]);
                            }
                            break;

                        }

                    }
                }

                /*for (i = 0; i < completePath.size(); ++i) {
                    final.push_back(original[order[completePath[i]]]);
                }*/
                feasible = 1;
                goto End;
                //END MAKE PATH FUNCTION
                //endregion

            }
            else if (SSum < numCycles) {
                feasible = 0;
                goto End;

            }
            else {
                feasible = 0;
                goto End;
            }


        }
    }


    End:
    
    //cout << "end:\n";
    if (feasible == 1) {
        //cout << "Final:\n";
        /*for (int v : final) {
            cout << v << " ";
        }*/
        cout << endl;
    }

}

cout << "COMPLETE.\n";


}

/*
    for (i = 0; i < Tpatch.size(); ++i) {
        for (j = 0; j < Tpatch[i].size(); ++j) {
            patchVertex[Tpatch[i][j]] = i;
            patchVertex[matchList[Tpatch[i][j]]] = i;
        }
    }

    u = 0;
    v = 0;
    save = 0;
    for (j = 0; j < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++j) {
        if (mateInduced[cycleVertex[Tpatch[u][v]]][j] == matchList[Tpatch[u][v]]) {
            save = j;
            break;
        }
    }
    //region multipleTCyclePatch
    //MULTIPLE T CYCLE PATCH FUNCTION ***********************
    int x = 0;

    //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
    if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == Tpatch[u][v]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if (x == 0) {
            for (i = 0; i < save; ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
    else if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == matchList[Tpatch[u][v]]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if (x == 0) {
            for (i = 0; i < save; ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
    else if (mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == matchList[Tpatch[u][v]]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = save; i-- > 0;) { //from element at position 'save' to the first element in the cycle
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if (x == 0) {
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size();
                 i-- > save + 1;) { //from end of cycle to element at position save+1
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
    else if (mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == Tpatch[u][v]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = save; i-- > 0;) { //from element at position 'save' to the first element in the cycle
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if (x == 0) {
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size();
                 i-- > save + 1;) { //from end of cycle to element at position save+1
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
        //the last element in the cycle is Tpatch[u][v]
    else if (save == 0 &&
             mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1] ==
             Tpatch[u][v]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
        //the first element in the cycle is Tpatch[u][v]
    else if (save == mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1 &&
             mateInduced[cycleVertex[Tpatch[u][v]]][0] == Tpatch[u][v]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1; i-- > 0;) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
        //and the last element in the cycle is matchList[Tpatch[u][v]]
    else if (save == 0 &&
             mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1] ==
             matchList[Tpatch[u][v]]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
        //the first element in the cycle is matchList[Tpatch[u][v]]
    else if (save == mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1 &&
             mateInduced[cycleVertex[Tpatch[u][v]]][0] == matchList[Tpatch[u][v]]) {
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1; i-- > 0;) {
            if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }
    //END MULTIPLE T CYCLE PATCH FUNCTION
    //endregion
    */

   /* while (fullCycle.size() < nScores) {
        save = 0;
        for (i = 0; i < Tpatch[u].size(); ++i) {
            if (Tpatch[u][i] == fullCycle.back()) {
                if (i == Tpatch[u].size() - 1) {
                    v = 0;
                }
                else {
                    v = ++i;
                }
                for (j = 0; j < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++j) {
                    if (mateInduced[cycleVertex[Tpatch[u][v]]][j] == matchList[Tpatch[u][v]]) {
                        save = j;
                        break;
                    }
                }
                break;
            }
            else if (matchList[Tpatch[u][i]] == fullCycle.back()) {
                if (i == 0) {
                    v = Tpatch[u].size() - 1;
                }
                else {
                    v = --i;
                }
                for (j = 0; j < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++j) {
                    if (mateInduced[cycleVertex[Tpatch[u][v]]][j] == Tpatch[u][v]) {
                        save = j;
                        break;
                    }
                }
                break;
            }
        }
        //region multipleTCyclePatch
        //MULTIPLE T CYCLE PATCH FUNCTION
        int x = 0;

        //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
        if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == Tpatch[u][v]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    x = 1;
                    break;
                }
            }
            if (x == 0) {
                for (i = 0; i < save; ++i) {
                    if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    }
                    else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                        u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                        break;
                    }
                }
            }
        }

            //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
        else if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == matchList[Tpatch[u][v]]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    x = 1;
                    break;
                }
            }
            if (x == 0) {
                for (i = 0; i < save; ++i) {
                    if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    }
                    else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                        u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                        break;
                    }
                }
            }
        }

            //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
        else if (mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == matchList[Tpatch[u][v]]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = save; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    x = 1;
                    break;
                }
            }
            if (x == 0) {
                for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size();
                     i-- > save + 1;) { //from end of cycle to element at position save+1
                    if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    }
                    else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                        u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                        break;
                    }
                }
            }
        }

            //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
        else if (mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == Tpatch[u][v]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = save; i-- > 0;) { //from element at position 'save' to the first element in the cycle
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    x = 1;
                    break;
                }
            }
            if (x == 0) {
                for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size();
                     i-- > save + 1;) { //from end of cycle to element at position save+1
                    if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    }
                    else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                        u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                        break;
                    }
                }
            }
        }

            //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
            //the last element in the cycle is Tpatch[u][v]
        else if (save == 0 && mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1] == Tpatch[u][v]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }

            //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
            //the first element in the cycle is Tpatch[u][v]
        else if (save == mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1 &&
                 mateInduced[cycleVertex[Tpatch[u][v]]][0] == Tpatch[u][v]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1; i-- > 0;) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }

            //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
            //and the last element in the cycle is matchList[Tpatch[u][v]]
        else if (save == 0 && mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1] == matchList[Tpatch[u][v]]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }

            //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
            //the first element in the cycle is matchList[Tpatch[u][v]]
        else if (save == mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1 &&
                 mateInduced[cycleVertex[Tpatch[u][v]]][0] == matchList[Tpatch[u][v]]) {
            fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size() - 1; i-- > 0;) {
                if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if (patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
        //END MULTIPLE T CYCLE PATCH FUNCTION
        //endregion

    }*/





