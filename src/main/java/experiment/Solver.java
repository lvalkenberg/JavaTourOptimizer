package experiment;

import core.TSP;

public class Solver {
    public static void main(String[] args) {
        // Check argument count
        if (args.length != 11) {
            System.out.println("Incorrect number of arguments! Please provide 11 arguments in the following order:");
            System.out.println("margCostFilter, repCostFilter, branchingStrategy, searchStrategy, lastConflictSearch, incrementalPi, fixedPi, keepBestPi, initialHeuristic, filePath, timeout");
            return;
        }
        // Convert string arguments to appropriate types
        boolean margCostFilter = Boolean.parseBoolean(args[0]);
        boolean repCostFilter = Boolean.parseBoolean(args[1]);
        String branchingStrategy = args[2];
        String searchStrategy = args[3];
        boolean lastConflictSearch = Boolean.parseBoolean(args[4]);
        boolean incrementalPi = Boolean.parseBoolean(args[5]);
        boolean fixedPi = Boolean.parseBoolean(args[6]);
        boolean keepBestPi = Boolean.parseBoolean(args[7]);
        boolean initialHeuristic = Boolean.parseBoolean(args[8]);
        double timeout = Double.parseDouble(args[9]);
        String filePath = args[10];

        // Create TSP solver using builder pattern
        TSP.TSPBuilder tspSolver = new TSP.TSPBuilder()
                .setMargCostFilter(margCostFilter)
                .setRepCostFilter(repCostFilter)
                .setBranchingStrategy(branchingStrategy)
                .setSearchStrategy(searchStrategy)
                .setLastConflictSearch(lastConflictSearch)
                .setIncrementalPi(incrementalPi)
                .setFixedPi(fixedPi)
                .setKeepBestPi(keepBestPi)
                .setInitialHeuristic(initialHeuristic)
                .setTimeout(timeout);

        // Create TSP instance with solver and file path
        try {
            Long time = System.currentTimeMillis();
            TSP tsp = new TSP(tspSolver, filePath);
            tsp.solve();
            System.out.print(System.currentTimeMillis() - time);
            System.out.print(' ');
            System.out.print(tsp.visitedNodeCounter);
        } catch (Exception e) {
            System.out.println(e);
        }

        // Add the rest of your TSP solving code here
    }
}
