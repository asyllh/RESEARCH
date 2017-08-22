/*--------------------/
ALH
mbahra.cpp
18/08/2017
/--------------------*/
#include <algorithm>
#include "mbahra.h"
using namespace std;

void resetVectors(int vacant, int numScores, int numComp, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList, vector<int> &mates, vector<vector<int> > &S, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes){

    int i, j;

    for(i = 0; i < numScores; ++i){
        //allScores[i] = 0;
        cycleVertex[i] = 1;
        matchList[i] = vacant;
        mates[i] = 0;
        for(j = 0; j < numScores; ++j){
            adjMatrix[i][j] = 0;
            boxWidths[i][j] = 0;
            allBoxes[i][j] = 0;
        }
    }

    for(i = 0; i < numComp; ++i){
        for(j = 0; j < numComp; ++j){
            S[i][j] = 0;
        }
    }
    /*for(i = 0; i < numScores - 2; ++i){
        matchList[i] = vacant;
    }*/

}

void clearVectors(vector<int> &allScores, vector<vector<int> > &mateInduced, vector<int> &lengthMateInduced, vector<vector<int> > &T, vector<int> &fullCycle, vector<int> &completePath){

    mateInduced.clear();
    lengthMateInduced.clear();
    T.clear();
    fullCycle.clear();
    completePath.clear();
    allScores.clear();
}

void createInstance(int threshold, int minWidth, int maxWidth, int minBoxWidth, int maxBoxWidth, int numScores, int numBox, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes){

    int i, j, k;
    vector<int> randOrder;

    //Create random values to be used as score widths, put in allScores vector (except last two elements)
    for (i = 0; i < numScores - 2; ++i) {
        allScores.push_back(rand() % (maxWidth - minWidth + 1) + minWidth);
    }
    //add two dominating vertices with score widths = 71 (these scores will be either side of same box, mates)
    allScores.push_back(70);
    allScores.push_back(70);



    //Sort all of the scores in the allScores vector in ascending order
    sort(allScores.begin(), allScores.end()); //sorts elements of vector in ascending order

    cout << "All scores:\n";
    for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl <<endl;

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

    cout << "AdjMatrix:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;



    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                mates[i] = j;
                break;
            }
        }
    }
    cout << "Mates Vector:\n";
    for(i = 0; i < mates.size(); ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            if(adjMatrix[i][j] == 2 && boxWidths[i][j] == 0){
                boxWidths[i][j] = rand() % (maxBoxWidth - minBoxWidth + 1) + minBoxWidth;
                boxWidths[j][i] = boxWidths[i][j];
                break;
            }

        }
    }

    boxWidths[numScores - 1][numScores - 2] = 0;
    boxWidths[numScores - 2][numScores - 1] = 0;

    cout << "Box Widths:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << boxWidths[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    k = 1;
    for(i = 0; i < numScores; ++i){
        for(j = i+1; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                allBoxes[i][j] = k;
                allBoxes[j][i] = k * -1;
                ++k;
                break;
            }
        }
    }
    cout << endl;

    /*cout << "allBoxes:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allBoxes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}

void createInstanceUser(int threshold, int numScores, vector<int> &allScores, vector<vector<int> > &userInput, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<vector<int> > &boxWidths, vector<vector<int> > &allBoxes){

    int i, j, m1, m2, k;

    for(i = 0; i < userInput.size(); ++i){
        for(j = 0; j < userInput[i].size() - 1; ++j){
            allScores.push_back(userInput[i][j]);
        }
    }

    allScores.push_back(70);
    allScores.push_back(70);

    sort(allScores.begin(), allScores.end());

    cout << "AllScores User:\n";
    for(i = 0; i < allScores.size(); ++i){
        cout << allScores[i] << " ";
    }
    cout << endl << endl;

    for (i = 0; i < allScores.size() - 1; ++i) {
        for (j = i + 1; j < allScores.size(); ++j) {
            if (allScores[i] + allScores[j] >= threshold) {
                adjMatrix[i][j] = 1;
                adjMatrix[j][i] = 1;
            }
        }

    }

    vector<int>::iterator it1;
    vector<int>::iterator it2;

    for(i = 0; i < userInput.size(); ++i) {
        it1 = find(allScores.begin(), allScores.end(), userInput[i][0]);
        m1 = it1 - allScores.begin();
        it2 = find(allScores.begin(), allScores.end(), userInput[i][1]);
        m2 = it2 - allScores.begin();
        //cout << "position of pairs " << userMates[i][0] << " and " << userMates[i][1] << ": " << m1 << "-" << m2 << endl;
        adjMatrix[m1][m2] = 2;
        adjMatrix[m2][m1] = 2;
        boxWidths[m1][m2] = userInput[i][2];
        boxWidths[m2][m1] = userInput[i][2];
    }
    adjMatrix[allScores.size()-1][allScores.size()-2] = 2;
    adjMatrix[allScores.size()-2][allScores.size()-1] = 2;

    cout << "AdjMatrix:\n";
    for(i = 0; i < adjMatrix.size(); ++i){
        for(j = 0; j < adjMatrix[i].size(); ++j){
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "BoxWidths:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << boxWidths[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;

    for (i = 0; i < numScores; ++i) {
        for (j = 0; j < numScores; ++j) {
            if (adjMatrix[i][j] == 2) {
                mates[i] = j;
                break;
            }
        }
    }

    cout << "Mates:\n";
    for(i = 0; i < mates.size(); ++i){
        cout << mates[i] << " ";
    }
    cout << endl << endl;

    k = 1;
    for(i = 0; i < numScores; ++i){
        for(j = i+1; j < numScores; ++j){
            if(adjMatrix[i][j] == 2){
                allBoxes[i][j] = k;
                allBoxes[j][i] = k * -1;
                ++k;
                break;
            }
        }
    }

    /*cout << "allBoxes:\n";
    for(i = 0; i < numScores; ++i){
        for(j = 0; j < numScores; ++j){
            cout << allBoxes[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;*/

}

void MTGMA(int vacant, int threshold, int numScores, int &matchSize, vector<int> &allScores, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList){

    int i, j, k;
    int lastMatch = vacant;
    int mateMatch = vacant;
    int vacantFlag = 0;
    matchSize = 0;

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
                }
                else if (adjMatrix[i][j] == 2 && matchList[j] == vacant) { //if potential match == mate
                    vacantFlag = 1;
                }
            }//end for j
            if (matchList[i] == vacant) { //if vertex has still not been matched
                for (k = 0; k < numScores - 2; ++k) {
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


    cout << "Cycle Vertex vector after MTGMA:\n";
    for(i = 0; i < cycleVertex.size(); ++i){
        cout << cycleVertex[i] << " ";
    }
    cout << endl << endl;

    cout << "Matching List:\n";
    for(i = 0; i < matchList.size(); ++i){
        cout << matchList[i] << " ";
    }
    cout << endl << endl;


}

void MIS(int numScores, int &numCycles, vector<vector<int> > &adjMatrix, vector<int> &mates, vector<int> &matchList, vector<vector<int> > &mateInduced, vector<int> &lengthMateInduced){

    int i, j;
    int smallestVertex;
    int currentVertex;
    vector<int> tempMIS;
    vector<int> checked(numScores, 0);
    numCycles = 0;


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
            tempMIS.push_back(currentVertex);
            checked[currentVertex] = 1;
            tempMIS.push_back(mates[currentVertex]);
            checked[mates[currentVertex]] = 1;
            currentVertex = matchList[mates[currentVertex]];
        } while (currentVertex != smallestVertex);

        mateInduced.push_back(tempMIS);
        tempMIS.clear();

        for (i = 0; i < numScores; ++i) {
            if (checked[i] == 0) {
                smallestVertex = i;
                break;
            }
        }


    } while (smallestVertex != currentVertex);

    tempMIS.clear(); //clear cycle vector again for next instance

    numCycles = mateInduced.size(); //number of cycles in the mate-induced structure

    cout << "Mate-Induced Structure:\n";
    for(i = 0; i < mateInduced.size(); ++i){
        for(j = 0; j < mateInduced[i].size(); ++j){
            cout << mateInduced[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl << endl;
    //cout << "Number of cycles in mate-induced structure: " << numCycles << endl;

    for (i = 0; i < mateInduced.size(); ++i) {
        lengthMateInduced.push_back(mateInduced[i].size());
    }


}

void FCA(int &qstar, int vacant, int matchSize, vector<vector<int> > &adjMatrix, vector<int> &cycleVertex, vector<int> &matchList, vector<vector<int> > &mateInduced, vector<vector<int> > &S, vector<vector<int> > &T){

    int i, j, k;
    int numEdges;
    vector<int> edge;
    vector<int> t;

    //FCA
    //create list cycleVertex that contains for each vertex the cycle that each edge belongs to
    for (i = 0; i < mateInduced.size(); ++i) {
        for (j = 0; j < mateInduced[i].size(); ++j) {
            if (cycleVertex[mateInduced[i][j]] != vacant) { //if edge is not deleted for FCA
                cycleVertex[mateInduced[i][j]] = i;
            }
        }
    }
    /*cout << "Cycle Vertex:\n";*/
    //create list of edges without empty edges (those generated by mate swap)
    for (i = 0; i < matchSize; ++i) {
        while (cycleVertex[i] == vacant) {
            ++i;
        }
        edge.push_back(i);
    }
    numEdges = edge.size();
    /*cout << "Edges vector:\n";*/
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



}

void oneTCyclePatch(int v, int full, int save, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &fullCycle, vector<vector<int> > &mateInduced, vector<vector<int> > &T){

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

void multipleTCyclePatch(int &u, int v, int save, int vacant, vector<int> &matchList, vector<int> &cycleVertex, vector<int> &fullCycle, vector<vector<int> > &mateInduced, vector<vector<int> > &Tpatch, vector<int> &patchVertex){

    int i;
    int x = 0;

    //CASE ONE: if the current value is matchList[Tpatch[u][v]] and the previous value is Tpatch[u][v]
    if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == Tpatch[u][v]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if(x == 0) {
            for (i = 0; i < save; ++i) {
                if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE TWO: if the current value is Tpatch[u][v] and the previous value is matchList[Tpatch[u][v]]
    else if (mateInduced[cycleVertex[Tpatch[u][v]]][save - 1] == matchList[Tpatch[u][v]]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = save + 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if(x == 0) {
            for (i = 0; i < save; ++i) {
                if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE THREE: if the current vertex is Tpatch[u][v] and the next element is matchList[Tpatch[u][v]]
    else if(mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == matchList[Tpatch[u][v]]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = save; i-- > 0;){ //from element at position 'save' to the first element in the cycle
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if(x == 0) {
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size(); i-- > save + 1;) { //from end of cycle to element at position save+1
                if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE FOUR: if the current vertex is matchList[Tpatch[u][v]] and the next element is Tpatch[u][v]
    else if(mateInduced[cycleVertex[Tpatch[u][v]]][save + 1] == Tpatch[u][v]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = save; i-- > 0;){ //from element at position 'save' to the first element in the cycle
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                x = 1;
                break;
            }
        }
        if(x == 0) {
            for (i = mateInduced[cycleVertex[Tpatch[u][v]]].size(); i-- > save + 1;) { //from end of cycle to element at position save+1
                if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                }
                else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                    fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                    u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                    break;
                }
            }
        }
    }

        //CASE FIVE: if the current vertex is matchList[Tpatch[u][v]] and is the first element in the cycle, and
        //the last element in the cycle is Tpatch[u][v]
    else if(save == 0 && mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size()-1] == Tpatch[u][v]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE SIX: if the current vertex is matchList[Tpatch[u][v]] and is last element in the cycle, and
        //the first element in the cycle is Tpatch[u][v]
    else if(save == mateInduced[cycleVertex[Tpatch[u][v]]].size()-1 && mateInduced[cycleVertex[Tpatch[u][v]]][0] == Tpatch[u][v]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = mateInduced[cycleVertex[Tpatch[u][v]]].size()-1; i-- > 0; ){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE SEVEN: if the current vertex is Tpatch[u][v] and is the first element in the cycle, and
        //and the last element in the cycle is matchList[Tpatch[u][v]]
    else if(save == 0 && mateInduced[cycleVertex[Tpatch[u][v]]][mateInduced[cycleVertex[Tpatch[u][v]]].size()-1] == matchList[Tpatch[u][v]]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = 1; i < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++i){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

        //CASE EIGHT: if the current vertex is Tpatch[u][v] and is last element in the cycle, and
        //the first element in the cycle is matchList[Tpatch[u][v]]
    else if(save == mateInduced[cycleVertex[Tpatch[u][v]]].size()-1 && mateInduced[cycleVertex[Tpatch[u][v]]][0] == matchList[Tpatch[u][v]]){
        fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][save]);
        for(i = mateInduced[cycleVertex[Tpatch[u][v]]].size()-1; i-- > 0; ){
            if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] == vacant) {
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
            }
            else if(patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]] != vacant){
                fullCycle.push_back(mateInduced[cycleVertex[Tpatch[u][v]]][i]);
                u = patchVertex[mateInduced[cycleVertex[Tpatch[u][v]]][i]];
                break;
            }
        }
    }

}

void makePath(int numScores, vector<int> &fullCycle, vector<int> &completePath, vector<vector<int> > &boxWidths, vector<int> &allScores, vector<vector<int> > &allBoxes){

    int i, j;
    int totalLength = 0;
    vector<int> completeScoresPath;
    vector<int> boxOrder;
    vector<int> boxWidthsOrder;


    for(i = 0; i < fullCycle.size()-1; ++i){
        if((fullCycle[i] == numScores - 1 && fullCycle[i+1] == numScores - 2) || (fullCycle[i] == numScores - 2 && fullCycle[i+1] == numScores - 1)){
            if(i == 0){ //if the dominating vertices are at the beginning of the fullCycle vector
                for(j = 2; j < fullCycle.size(); ++j){
                    completePath.push_back(fullCycle[j]);
                }
                break;
            }

            else if(i == fullCycle.size()-2){ //if the dominating vertices are at the end of the fullCycle vector
                for(j = 0; j < fullCycle.size()-2; ++j){
                    completePath.push_back(fullCycle[j]);
                }
                break;
            }
            else{ //if the dominating vertices are in the middle of the fullCycle vector
                for(j = i+2; j < fullCycle.size(); ++j){
                    completePath.push_back(fullCycle[j]);
                }
                for(j = 0; j < i; ++j){
                    completePath.push_back(fullCycle[j]);
                }
                break;

            }

        }
    }

    cout << "Complete Path:\n";
    for(i = 0; i < completePath.size(); ++i){
        cout << completePath[i] << " ";
    }
    cout << endl << endl;

    for(i = 0; i < completePath.size(); ++i){
        completeScoresPath.push_back(allScores[completePath[i]]);
    }

    /*cout << "Complete Path - Score Widths:\n";
    for(i = 0; i < completeScoresPath.size(); ++i){
        cout << completeScoresPath[i] << " ";
    }
    cout << endl << endl;*/



    for(i = 0; i < completePath.size() - 1; ++i){
        if(boxWidths[completePath[i]][completePath[i+1]] != 0){
            totalLength += boxWidths[completePath[i]][completePath[i+1]];
            boxOrder.push_back(allBoxes[completePath[i]][completePath[i+1]]);
            boxWidthsOrder.push_back(boxWidths[completePath[i]][completePath[i+1]]);
        }
    }

    cout << "Order of Boxes:\n";
    for(i = 0; i < boxOrder.size(); ++i){
        cout << boxOrder[i] << " ";
    }
    cout << endl << endl;

    cout << "Box Widths in Order:\n";
    for(i = 0; i < boxWidthsOrder.size(); ++i){
        cout << boxWidthsOrder[i] << " ";
    }
    cout << endl;

    cout << "Total Length of Path: " << totalLength << " millimeters.\n\n";

}

void patchGraph(int qstar, int vacant, int instance, int numScores, int numCycles, int &feasible, int &infeasible, int &fullT, int &splitT, int &noPatch, int &problem, vector<int> &matchList, vector<int> &cycleVertex, vector<vector<int> > &mateInduced, vector<vector<int> > &S, vector<vector<int> > &T, vector<int> &fullCycle, vector<int> &completePath, vector<vector<int> > &boxWidths, vector<int> &allScores, vector<vector<int> > &allBoxes){

    int i, j, q, u, v, save, SSum, SqIntS;
    int full = vacant;
    vector<int> QSet(qstar, 0);
    vector<int> SSet;
    vector<int> patchCycle(qstar, vacant);
    vector<int> tempPG;
    vector<int> patchVertex(numScores, vacant);
    vector<vector<int> > Tpatch;

    for (i = 0; i < T.size(); ++i) {
        if (T[i].size() == numCycles) {
            full = i;
            break;
        }
    }

    if(full != vacant){
        save = 0;
        //cout << "Full: " << full << endl;
        for (v = 0; v < T[full].size(); ++v) {
            for (j = 0; j < mateInduced[cycleVertex[T[full][v]]].size(); ++j) {
                if (mateInduced[cycleVertex[T[full][v]]][j] == matchList[T[full][v]]) {
                    save = j;
                    break;
                }
            }
            oneTCyclePatch(v, full, save, matchList, cycleVertex, fullCycle, mateInduced, T);

        }
        /*cout << instance << ": Full Cycle (one T-cycle):\n";
        for (i = 0; i < fullCycle.size(); ++i) {
            cout << fullCycle[i] << " ";
        }
        cout << endl << endl;*/

        makePath(numScores, fullCycle, completePath, boxWidths, allScores, allBoxes);

        ++feasible;
        ++fullT;

    }

    else {
        q = 0; //Start with first Tq-cycle
        QSet[0] = 1;
        SSum = 0; //number of MIS-cycles that have been included
        for (i = 0; i < numCycles; ++i) {
            SSet.push_back(S[q][i]); // ==1 if MIS cycle i has been included
        }
        for (i = 0; i < numCycles; ++i) {
            SSum = SSum + SSet[i];
        }

        if (SSum >= 1) {
            patchCycle[q] = 1;
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
                        patchCycle[q] = 1;
                    }
                }
                QSet[q] = 1;
                q = 0;
            }
        }//end while


        //If patching graph is connected, then instance is feasible, else infeasible
        if (SSum == numCycles) {
            for (i = 0; i < patchCycle.size(); ++i) {
                if (patchCycle[i] == 1) {
                    for (j = 0; j < T[i].size(); ++j) {
                        tempPG.push_back(T[i][j]);
                    }
                    Tpatch.push_back(tempPG);
                    tempPG.clear();
                }
            }
            tempPG.clear();

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
            multipleTCyclePatch(u, v, save, vacant, matchList, cycleVertex, fullCycle, mateInduced, Tpatch, patchVertex);

            while (fullCycle.size() < numScores) {
                save = 0;
                for (i = 0; i < Tpatch[u].size(); ++i) {
                    if (Tpatch[u][i] == fullCycle.back()) {
                        if (i == Tpatch[u].size() - 1) {
                            v = 0;
                        } else {
                            v = ++i;
                        }
                        for (j = 0; j < mateInduced[cycleVertex[Tpatch[u][v]]].size(); ++j) {
                            if (mateInduced[cycleVertex[Tpatch[u][v]]][j] == matchList[Tpatch[u][v]]) {
                                save = j;
                                break;
                            }
                        }
                        break;
                    } else if (matchList[Tpatch[u][i]] == fullCycle.back()) {
                        if (i == 0) {
                            v = Tpatch[u].size() - 1;
                        } else {
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
                multipleTCyclePatch(u, v, save, vacant, matchList, cycleVertex, fullCycle, mateInduced, Tpatch, patchVertex);
            }

            /*cout << instance << ": Full Cycle (Mutiple T-cycles):\n";
            for (i = 0; i < fullCycle.size(); ++i) {
                cout << fullCycle[i] << " ";
            }
            cout << endl << endl;*/

            makePath(numScores, fullCycle, completePath, boxWidths, allScores, allBoxes);

            ++feasible;
            ++splitT;
        }
        else if (SSum < numCycles) {
            //cout << instance << ": Infeasible SSum < numCycles\n\n";
            ++noPatch;
            ++infeasible;
        }
        else {
            //cout << instance << ": Problem.\n\n";
            ++problem;
        }


    }

}

