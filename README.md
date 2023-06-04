# JavaTourOptimizer
This repository is a Java exact solver for the Travelling Salesman Problem ([TSP](https://en.wikipedia.org/wiki/Travelling_salesman_problem)). The core solver use a branch-and-bound with the a [1-tree based Lagrangian relaxation]([The Traveling-Salesman Problem and Minimum Spanning Trees](https://www.jstor.org/stable/169411)) as presented by M. Held and R. M. Karp. Other algorithms and heuristics are use with the goal of speed up the search of an optimum TSP tour or at least find a good tour rapidly.

## Installation

The [build.gradle](build.gradle) file is given and can be used to open the project in your IDE (Intellij for exemple) .

## Usage
Simple creation of an instance with the best parameter:
```java
// Creation of a TSP instance from an xml file (format : http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/XML-TSPLIB/Description.pdf)
TSP tsp = new TSP("pla85900.xml");

// From an EUC_2D tsp file (format : http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf)
TSP tsp = new TSP("pla85900.tsp");

// Or from a distance matrice
double[] distanceMatrix;
TSP tsp = new TSP(distanceMatrix);
```
Solve the TSP instance :
```java
tsp.solve();
```
A builder class *TSPBuilder* allow to give specific solving parameters a TSP instance. These parameters will change the background behavior of the solver. Here a the most interesting parameters
```java
TSP.TSPBuilder tspSolver = new TSP.TSPBuilder()
                .setMargCostFilter(boolean margCostFilter)                        // apply marginal cost filtering
                .setRepCostFilter(boolean repCostFilter)                          // apply replacement cost filtering
                .setBranchingStrategy(String branchingStrategy)                   // branching strategy (maxCost, minCost, nearest neighbour, minCost path)
                .setSearchStrategy(String searchStrategy)                         // search strategy (DFS, BFS, DFS on BFS)
                .setLastConflictSearch(boolean lastConflictSearch)                // use last conflict search
                .setIncrementalPi(boolean incrementalPi)                          // reuse the parent Lagrangian multiplier
                .setInitialHeuristic(boolean initialHeuristic)                    // heuristic to find a initial solution used a upper bound
                .setTimeout(int timeout);                                         // timeout in second
                .setOnBetterSolutionFound(Consumer<int[]> onBetterSolutionFound)  // callback function all each time a new tour is found
```
## License

The source code for the site is licensed under the MIT license, which you can find in the [LICENSE.txt](LICENSE.txt) file.
