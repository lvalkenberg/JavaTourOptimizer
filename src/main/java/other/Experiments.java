package other;

import core.TSP;

import java.io.IOException;
import java.util.Arrays;

public class Experiments {

    public static void main(String[] args) throws IOException {
        //TSP tsp = new TSP("../TSPlib/xml files/lin105.xml");
        TSP tsp = new TSP("../TSPlib/tsp files/ch130.tsp");
        tsp.verbose = true;
        tsp.timeout = -1;

        long start = System.nanoTime();
        System.out.println();
        System.out.println(Arrays.toString(tsp.solve()));
        System.out.println("\n Solved in " + (System.nanoTime() - start) / 1e6 + " ms");
    }

}
