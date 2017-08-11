	
**//VARIABLE FROM ARGUMENTS:**
* int numInstances = atoi(argv[1]); //number of instances of mssp, use in main for loop
* int numBox = atoi(argv[2]) + 1; //number of boxes in mssp plus 1 extra box (scores on either side of extra box will be dominating vertices, score widths = 71)
* int minWidth = atoi(argv[3]); //minimum width of scores (millimeters)
* int maxWidth = atoi(argv[4]); //maximum width of scores (millimeters)
* int minBoxWidth = atoi(argv[5]); //min box width (mm)
* int maxBoxWidth = atoi(argv[6]); //max box width (mm)
* int randomSeed = atoi(argv[7]); //random seed


**//VARIABLES**
* int i, j, k, q;
* int instance; //counter for instances loop
* int numScores = numBox * 2; //number of scores, 2 per box (1 either side), last two scores are dominating vertices
* int numComp = (numBox + (numBox % 2)) / 2;
* int threshold = 70; //adjacency threshold of scores, minimum knife distance
* int vacant = 999; //large empty value
* int feasible = 0; //number of feasible instances
* int infeasible = 0; //number of infeasible instances
* int noMatch = 0; //number of instances with |M| < n (therefore immediately infeasible)
* int oneCycle = 0; //number of instances where the MIS consists of only one cycle (therefore immediately feasible)
* int noFam = 0; //number of instances with no family of T-cycles (qstar = -1, therefore immediately infeasible)
* int noPatch = 0; //number of instances where T-cycles do not produce a connected patching graph (SSum < numCycles, therefore infeasible)
* int fullT = 0; //number of instances where only one T-cycle is required to connected all cycles in MIS (patching graph is connected using only one T-cycle, therefore feasible)
* int splitT = 0; //number of instances where multiple T-cycles are required to connected all cycles in MIS (patching graph is connected using multiple T-cycles, therefore feasible)
* int problem = 0; //number of problematic instances, SSum > numCycles (ERROR)
* int matchSize; //size (cardinality) of the matching list (matchList.size()) (&MTGMA, FCA)
* int numCycles; //number of cycles in the MIS (mateInduced.size()) (&MIS, patchGraph)
* int qstar; //number of T-cycles (&FCA, patchGraph)
* vector<int> allScores(numScores, 0); //vector containing all score widths (createInstance, MTGMA)
* vector<vector<int> > adjMatrix(numScores, vector<int>(numScores, 0)); //adjaceny matrix (createInstance, MTGMA, MIS, FCA)
* vector<int> mates(numScores, 0); //contains vertex index for mates, e.g if vertex 0 is mates with vertex 4, then mates[0] = 4 (createInstance, MIS)
* vector<int> matchList(numScores, vacant); //contains vertex index for matching vertices, e.g. if vertex 0 is matched with vertex 9, then matchList[0] = 9 (MTGMA, MIS, FCA, patchGraph)
* vector<int> cycleVertex(numScores, 1); //contains the number of the cycle of the mate-induced structure that the vertex i belongs to (MTGMA, FCA, patchGraph)
* vector<vector<int> > mateInduced; //each row of the matrix corresponds to one cycle, and contains the indices of the score widths (MIS, FCA, patchGraph)
* vector<int> lengthMateInduced; //each elements holds the value of the length of the corresponding cycle in the MIS (MIS)
* vector<vector<int> > S(numComp, vector<int>(numComp, 0)); // == 1 if edge from cycle j is used in T-cycle q (T[q][j]) (FCA, patchGraph)
* vector<vector<int> > T; //each row hold the lower vertex of the edges that make up one T-cycle
* vector<vector<int> > boxWidths(numScores, vector<int>(numScores, 0)); // holds widths in mm for the widths of boxes
* vector<int> fullCycle; //holds initial cycle of final solution
* vector<int> completePath; //holds final path, i.e. fullCycle without dominating vertices
* vector<vector<int> > allBoxes(numScores, vector<int>(numScores, vacant)); //contains values from 1 to numBox to denote which box the scores belong to
* //if value is positive, i.e "3", then the 3rd box has the smallest score on the LHS, and the larger score on the RHS
* //if value is negative, i.e "-3", then the 3rd box is rotated, i.e. smallest score on RHS and larger score on LHS.

//VARIABLES IN FUNCTIONS:

**makePath(numScores, fullCycle, completePath, boxWidths, allScores, allBoxes)**
* int totalLength = 0; //holds total length of completePath, i.e. total of all box widths in the path
* vector<int> completeScoresPath; // output actual score widths in order of completePath
* vector<int> boxOrder; // output order of boxes from 1 to numBox -1
* vector<int> boxWidthsOrder; // outputs all box Widths in order (in millimeters)
* **Method:**
* creates completePath
* creates completeScoresPath
* creates boxOrder, boxWidthsOrder, calculates totalLength

**createInstance(threshold, minWidth, maxWidth, minBoxWidth, maxBoxWidth, numScores, numBox, allScores, adjMatrix, mates, boxWidths, allBoxes)**
* vector<int> randOrder; //vector that creates random permutations of the scores that will be used to assign pairs of scores as mates
* **Method:**
* creates allScores vector, 
* sorts allScores vector
* fills in adjMatrix for sum >= threshold (=1)
* creates randOrder and shuffles randOrder
* assigns pairs of scores as mates
* fills in adjMatrix with mates (=2)
* creates mates vector
* creates boxWidths matrix
* creates allBoxes matrix

**MTGMA(vacant, threshold, numScores, &matchSize, allScores, adjMatrix, cycleVertex, matchList)**
* int lastMatch = vacant; 
* int mateMatch = vacant;
* **Method:** MTGMA only

**MIS(numScores, &numCycles, adjMatrix, mates, matchList, mateInduced, lengthMateInduced)**
* int smallestVertex;
* int currentVertex;
* vector<int> cycle; //temp vertex
* vector<int> checked(numScores, 0);
* **Method:** creates MIS

**FCA(&qstar, vacant, matchSize, adjMatrix, cycleVertex, matchList, mateInduced, S, T)**
* int numEdges;
* vector<int> edge;
* vector<int> t;
* **Method:**
* fill cycleVertex with number of MIS cycle
* create edges vector
* FCA - attempt to find T-cycles

**patchGraph(qstar, vacant, instance, numScores, numCycles, &feasible, &infeasible, &fullT, &splitT, &noPatch, &problem, matchList, cycleVertex, mateInduced, S, T, fullCycle, completePath, boxWidths, allScores, allBoxes)**
* int save, SSum, SqIntS;
* int full = vacant;
* vector<int> QSet(qstar, 0);
* vector<int> SSet;
* vector<int> patchCycle(qstar, vacant);
* vector<int> temp;
* vector<int> patchVertex(numScores, vacant);
* vector<vector<int> > Tpatch;
* **Method:**
* Patch Graph if one full T cycle exists
* Patch Graph if there are multiple T-cycles









































