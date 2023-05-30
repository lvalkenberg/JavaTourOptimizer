package heuristics;
import java.util.Arrays;
import java.util.Random;

public class Heuristic {

    public static Random rd = new Random();

    /**
     * Compute an approximation of a TSP tour using 2-opt local search
     * with restarting.
     *
     * @param distanceMatrix    distance matrix of the graph
     * @return an int[] containing the order of visit of each node
     */
    public static int[] solve(double[][] distanceMatrix) {
        int n = distanceMatrix.length;
        int[] tour = generateInitialTour(distanceMatrix);
        int[] bestTour = Arrays.copyOf(tour, n);
        double bestCost = tourCost(tour, distanceMatrix);

        int maxRestart = n*5;
        int nbRestart = 0;
        int unsuccesfulReastart = 0;

        boolean improvement = true;
        while (improvement) {
            improvement = false;
            for (int i = 0; i < n; i++) {
                for (int j = i + 2; j < n; j++) {
                    double delta = distanceMatrix[tour[i]][tour[j]] + distanceMatrix[tour[i + 1]][tour[(j + 1) % n]]
                            - distanceMatrix[tour[i]][tour[i + 1]] - distanceMatrix[tour[j]][tour[(j + 1) % n]];
                    if (delta < -10e-9 * bestCost) {
                        reverse(tour, i + 1, j);
                        improvement = true;
                    }
                }
            }
            if(!improvement && nbRestart < maxRestart){
                unsuccesfulReastart ++;
                restartRandomTour(tour);
                nbRestart ++;
                improvement = true;
            }

            if (tourCost(tour, distanceMatrix) < tourCost(bestTour, distanceMatrix)) {
                bestTour = Arrays.copyOf(tour, n);
                unsuccesfulReastart = 0;
                bestCost = tourCost(bestTour, distanceMatrix);
            }
        }
        return bestTour;
    }

    /**
     * Given a tour shuffle it
     *
     * @param arr   tour to suffle
     */
    public static void restartRandomTour(int[] arr) {
        Random random = new Random();
        int n = arr.length;

        for (int i = 0; i < n; i++) {
            // Generate a random index between i and n-1
            int randomIndex = i + random.nextInt(n - i);

            // Swap arr[i] and arr[randomIndex]
            int temp = arr[i];
            arr[i] = arr[randomIndex];
            arr[randomIndex] = temp;
        }
    }


    /**
     * Generate an initial tour based on the nearest neighbour heuristic
     *
     * @param distanceMatrix    distance matrix representing the graph
     * @return a tour as int[] containing the node in the order of visit
     */
    private static int[] generateInitialTour(double[][] distanceMatrix) {
        int n = distanceMatrix.length;
        int[] tour = new int[n];
        boolean[] visited = new boolean[n];
        boolean tourCompleted = false;
        int current = 0;
        tour[0] = 0;
        visited[0] = true;

        int closest = -1;
        double closestDist = Double.MAX_VALUE;
        int k = 0;
        while(!tourCompleted){
            closest = -1;
            closestDist = Double.MAX_VALUE;
            for(int i = 0; i < n; i++){
                if(!visited[i] && (closest == -1 || distanceMatrix[current][i] < closestDist)){
                    closest = i;
                    closestDist = distanceMatrix[current][i];
                }
            }
            k ++;
            if(closest != -1){ current = closest; visited[closest] = true; tour[k] = closest;}
            else tourCompleted = true;
        }
        return tour;
    }

    /**
     * reverse the element if @tour between i and j (i<=j)
     *
     * @param tour  tour to apply the transformation
     * @param i     first index
     * @param j     second index
     */
    private static void reverse(int[] tour, int i, int j) {
        while (i < j) {
            int temp = tour[i];
            tour[i] = tour[j];
            tour[j] = temp;
            i++;
            j--;
        }
    }

    /**
     * Compute the cost of a tour
     *
     * @param tour              tour as an ordered list of node
     * @param distanceMatrix    distance matrix of the graph
     * @return  the cost of the tour on the graph
     */
    public static double tourCost(int[] tour, double[][] distanceMatrix) {
        double cost = 0;
        for (int i = 0; i < tour.length; i++) {
            cost += distanceMatrix[tour[i]][tour[(i + 1) % tour.length]];
        }
        return cost;
    }
}
