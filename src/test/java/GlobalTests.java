import static org.junit.jupiter.api.Assertions.*;

import core.TSP;
import org.junit.jupiter.api.Test;
import experiment.soTSP;

public class GlobalTests {
    @Test
    void AllParameterTests() {
        // Symetric instances of the TSPlib
        //String[] instances = {"burma14", "ulysses16", "gr17", "gr21", "gr24", "fri26", "bayg29", "bays29"}; // , "ulysses22"
        String[] instances = {"fri26", "bays29", "ulysses16"}; // , "bayg29"
        boolean[] booleanOptions = {true, false};
        String[] branchingStrategies = {"nearest neighbour", "maxCost", "minCost path"};
        String[] searchStrategies = {"BFS", "DFS", "DFS on BFS"};
        double[] TSPtour = {937, 2020, 6859};

        int i = -1;
        for (String instance : instances) {
            i++;
            for (boolean filter : booleanOptions) {
                for (String branchingStrategy : branchingStrategies) {
                    for (String searchStrategy : searchStrategies) {
                        for (boolean lastConflictSearch : booleanOptions) {
                            for (boolean incrementalPi : booleanOptions) {
                                for (boolean keepBestPi : booleanOptions) {
                                    for (boolean initialHeuristic : booleanOptions) {
                                        try {
                                            if (searchStrategy == "DFS" && branchingStrategy == "maxCost")
                                                continue; // terrible performances
                                            TSP TSP = new TSP.TSPBuilder()
                                                    .setMargCostFilter(filter)
                                                    .setRepCostFilter(filter)
                                                    .setBranchingStrategy(branchingStrategy)
                                                    .setSearchStrategy(searchStrategy)
                                                    .setLastConflictSearch(lastConflictSearch)
                                                    .setIncrementalPi(incrementalPi)
                                                    .setKeepBestPi(keepBestPi)
                                                    .setInitialHeuristic(initialHeuristic)
                                                    .build("../TSPlib/xml files/" + instance + ".xml");

                                            System.out.println(instance + " : "
                                                    + "filter=" + filter
                                                    + ", branchingStrategy=" + branchingStrategy
                                                    + ", searchStrategy=" + searchStrategy
                                                    + ", lastConflictSearch=" + lastConflictSearch
                                                    + ", incrementalPi=" + incrementalPi
                                                    + ", keepBestPi=" + keepBestPi
                                                    + ", initialHeuristic=" + initialHeuristic);
                                            TSP.solve();
                                            assertEquals(TSPtour[i], TSP.getLowerBound(), 10e-6);
                                            // Print the parameter values after the successful test
                                        } catch (Exception e) {
                                            fail(e);
                                        }

                                    }
                                }

                            }
                        }
                    }
                }

            }
        }
    }


    @Test
    void outputVerySmall() {
        // Symetric instances of the TSPlib
        String[] instances = {"burma14", "ulysses16", "gr17", "gr21", "gr24", "fri26", "bayg29", "bays29"}; // , "ulysses22"

        for (String instance : instances) {
            soTSP trueTSP = new soTSP();
            trueTSP.xmlReader("../TSPlib/xml files/" + instance + ".xml");
            trueTSP.verbose = false;
            trueTSP.solve();

            try {
                TSP TSP = new TSP("../TSPlib/xml files/" + instance + ".xml");
                TSP.verbose = false;
                TSP.solve();
                System.out.println(instance + " : " + TSP.getLowerBound());
                assertEquals(trueTSP.getLB(), TSP.getLowerBound(), 10e-6);
            } catch (Exception e) {
                fail(e);
            }
        }
    }


    //@Test
    void outputSmall() {
        // Symetric instances of the TSPlib
        String[] instances = {"swiss42", "berlin52", "st70", "eil76"}; // , "gr48", "brazil58"

        double simpleTSPTime = 0; // ms
        long simpleNode = 0;
        double filterTSPTime = 0; // ms
        long filterNode = 0;
        long start;

        for (String instance : instances) {
            soTSP trueTSP = new soTSP();
            trueTSP.xmlReader("../TSPlib/xml files/" + instance + ".xml");
            trueTSP.verbose = false;
            start = System.currentTimeMillis();
            trueTSP.solve();
            simpleTSPTime += (System.currentTimeMillis() - start);
            simpleNode += trueTSP.visited_node;

            try {
                TSP TSP = new TSP("../TSPlib/xml files/" + instance + ".xml");
                TSP.verbose = false;
                start = System.currentTimeMillis();
                TSP.solve();
                filterTSPTime += (System.currentTimeMillis() - start);
                filterNode += TSP.visitedNodeCounter;

                System.out.println(instance + "  :  " + TSP.getLowerBound());
                assertEquals(trueTSP.getLB(), TSP.getLowerBound(), 10e-6);
            } catch (Exception e) {
                fail(e);
            }
        }
    }

    //@Test
    void outputMedium() {
        // Symetric instances of the TSPlib
        String[] instances = {"dantzig42", "att48", "gr48", "rat99"}; // , "gr48", "brazil58"

        double simpleTSPTime = 0; // ms
        long simpleNode = 0;
        double filterTSPTime = 0; // ms
        long filterNode = 0;
        long start;

        for (String instance : instances) {
            soTSP trueTSP = new soTSP();
            trueTSP.xmlReader("../TSPlib/xml files/" + instance + ".xml");
            trueTSP.verbose = false;
            start = System.currentTimeMillis();
            trueTSP.solve();
            simpleTSPTime += (System.currentTimeMillis() - start);
            simpleNode += trueTSP.visited_node;

            try {
                TSP TSP = new TSP("../TSPlib/xml files/" + instance + ".xml");
                TSP.verbose = false;
                start = System.currentTimeMillis();
                TSP.solve();
                filterTSPTime += (System.currentTimeMillis() - start);
                filterNode += TSP.visitedNodeCounter;

                System.out.println(instance + " : " + TSP.getLowerBound());
                assertEquals(trueTSP.getLB(), TSP.getLowerBound(), 10e-6);
            } catch (Exception e) {
                fail(e);
            }
        }
    }

    //@Test
    void outputBigger() {
        // Symetric instances of the TSPlib
        String[] instances = {"rat99", "eil101", "lin105", "kroA100", "ch130"};
        // http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/STSP.html
        double[] TSPtour = {1211, 629, 14379, 21282, 6110};
        //double[] TSPtour = {1219.2437684, 640.2115908, 14382.9959333, 21282, 6110}; // !!! solution may varie depending of the rounding use

        double simpleTSPTime = 0; // ms
        long simpleNode = 0;
        double filterTSPTime = 0; // ms
        long filterNode = 0;
        long start;

        int index = 0;
        for (String instance : instances) {
            try {
                TSP TSP = new TSP("../TSPlib/tsp files/" + instance + ".tsp");
                TSP.verbose = false;
                start = System.currentTimeMillis();
                TSP.solve();
                filterTSPTime += (System.currentTimeMillis() - start);
                filterNode += TSP.visitedNodeCounter;

                System.out.println(instance + " : " + TSP.getLowerBound());
                assertEquals(TSPtour[index], TSP.getLowerBound(), 10e-6);
            } catch (Exception e) {
                fail(e);
            }
            index++;
        }
    }

    void outputSmallRandom() {
        // Symetric instances of the TSPlib
        String[] instances = {"random/rd020_mmhrP", "random/rd021_iAW5s", "random/rd022_cM9gb", "random/rd023_5tp4r", "random/rd024_mki7C", "random/rd025_g60Yb", "random/rd026_7XDE1", "random/rd027_IE0Vj", "random/rd028_6akpg", "random/rd029_Sfk9O", "random/rd030_wRKSe"}; // , "ulysses22"

        double simpleTSPTime = 0; // ms
        double filterTSPTime = 0; // ms
        long start;

        for (String instance : instances) {
            soTSP trueTSP = new soTSP();
            trueTSP.xmlReader("../TSPlib/xml files/" + instance + ".xml");
            trueTSP.verbose = false;
            start = System.currentTimeMillis();
            trueTSP.solve();
            simpleTSPTime += (System.currentTimeMillis() - start);

            try {
                TSP TSP = new TSP("../TSPlib/xml files/" + instance + ".xml");
                TSP.verbose = false;
                start = System.currentTimeMillis();
                TSP.solve();
                filterTSPTime += (System.currentTimeMillis() - start);

                assertEquals(trueTSP.getLB(), TSP.getLowerBound(), 10e-6);
            } catch (Exception e) {
                fail(e);
            }
        }
    }
}
