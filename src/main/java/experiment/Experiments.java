package experiment;

import core.TSP;

import java.io.IOException;
import java.util.Arrays;

public class Experiments {

    public static void main(String[] args) throws IOException {
        TSP.TSPBuilder builder = new TSP.TSPBuilder().setOnBetterSolutionFound(arr -> System.out.println(Arrays.toString(arr)));
        TSP tsp = new TSP(builder, "../TSPlib/xml files/fri26.xml");
        //TSP tsp = new TSP("../TSPlib/tsp files/bays29.tsp");
        tsp.verbose = true;
        tsp.timeout = -1;

        long start = System.nanoTime();
        System.out.println();
        tsp.solve();
        //System.out.println(Arrays.toString(tsp.solve()));
        System.out.println("\n Solved in " + (System.nanoTime() - start) / 1e6 + " ms");
    }

}
