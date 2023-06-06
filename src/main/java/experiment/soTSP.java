package experiment;// simple exact Solver.TSP solver based on branch-and-bound/Held--Karp

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * From on https://stackoverflow.com/questions/7159259/optimized-tsp-algorithms/7163961#7163961
 */
public class soTSP {
    // instance name
    private String name;
    // number of cities
    private int n;
    // city locations
    private double[] x;
    private double[] y;
    // cost matrix
    private double[][] cost;
    // matrix of adjusted costs
    private double[][] costWithPi;
    public Node bestNode = new Node();
    public boolean verbose = true;
    public long visited_node;


    public static void main(String[] args) throws IOException {
        // read the input in TSPLIB format
        // assume TYPE: Solver.TSP, EDGE_WEIGHT_TYPE: EUC_2D
        // no error checking
        System.err.printf("%n");
        soTSP tsp = new soTSP();
        //tsp.readInput(new InputStreamReader(System.in));
        tsp.readInput("../TSPlib/tsp files/lin105.tsp");

        long start = System.nanoTime();
        tsp.solve();
        if (tsp.verbose) System.out.println("\nSolved in " + (System.nanoTime() - start) / 1e6 + " ms");
    }

    /**
     * EUC_2D .tsp reader
     *
     * @param filePath path the the file to read .tsp format EUC_2D
     * @throws IOException Error with the provided path
     */
    public void readInput(String filePath) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(filePath));
        Pattern specification = Pattern.compile("\\s*([A-Z_]+)\\s*(:\\s*([0-9]+))?\\s*");
        Pattern data = Pattern.compile("\\s*([0-9]+)\\s+([-+.0-9Ee]+)\\s+([-+.0-9Ee]+)\\s*");
        String line;
        while ((line = in.readLine()) != null) {
            Matcher m = specification.matcher(line);
            if (!m.matches()) continue;
            String keyword = m.group(1);
            if (keyword.equals("DIMENSION")) {
                n = Integer.parseInt(m.group(3));
                cost = new double[n][n];
            } else if (keyword.equals("NODE_COORD_SECTION")) {
                x = new double[n];
                y = new double[n];
                for (int k = 0; k < n; k++) {
                    line = in.readLine();
                    m = data.matcher(line);
                    m.matches();
                    int i = Integer.parseInt(m.group(1)) - 1;
                    x[i] = Double.parseDouble(m.group(2));
                    y[i] = Double.parseDouble(m.group(3));
                }
                // TSPLIB distances are rounded to the nearest integer to avoid the sum of square roots problem
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        double dx = x[i] - x[j];
                        double dy = y[i] - y[j];
                        //cost[i][j] = Math.rint(Math.sqrt(dx * dx + dy * dy));
                        cost[i][j] = Math.round(Math.sqrt(dx * dx + dy * dy) * 1e7) / 1e7;
                    }
                }
            }
        }
    }

    /**
     * Read Solver.TSP Instance from xml
     * See http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/XML-TSPLIB/Description.pdf
     *
     * @param xmlPath path to the file
     */
    public void xmlReader(String xmlPath) {
        String[] elem = xmlPath.split("/");
        this.name = elem[elem.length - 1].split("[.]")[0];
        // Instantiate the Factory
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        try {

            // optional, but recommended
            // process XML securely, avoid attacks like XML External Entities (XXE)
            dbf.setFeature(XMLConstants.FEATURE_SECURE_PROCESSING, true);

            // parse XML file
            DocumentBuilder db = dbf.newDocumentBuilder();

            Document doc = db.parse(new File(xmlPath));
            doc.getDocumentElement().normalize();

            NodeList list = doc.getElementsByTagName("vertex");

            n = list.getLength();
            this.cost = new double[n][n];

            for (int i = 0; i < n; i++) {
                NodeList edgeList = list.item(i).getChildNodes();
                for (int v = 0; v < edgeList.getLength(); v++) {

                    org.w3c.dom.Node node = edgeList.item(v);
                    if (node.getNodeType() == org.w3c.dom.Node.ELEMENT_NODE) {
                        Element element = (Element) node;
                        String cost = element.getAttribute("cost");
                        String adjacentNode = element.getTextContent();
                        int j = Integer.parseInt(adjacentNode);
                        //distanceMatrix[i][j] = Math.rint(Double.parseDouble(cost)); // Rounded !
                        this.cost[i][j] = Double.parseDouble(cost);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void solve() {
        bestNode.lowerBound = Double.MAX_VALUE;
        Node currentNode = new Node();
        currentNode.excluded = new boolean[n][n];
        costWithPi = new double[n][n];
        computeHeldKarp(currentNode);
        PriorityQueue<Node> pq = new PriorityQueue<Node>(11, new NodeComparator());
        visited_node = 0;
        do {
            do {
                visited_node++;
                int i = -1;
                for (int j = 0; j < n; j++) {
                    if (currentNode.degree[j] > 2 && (i < 0 || currentNode.degree[j] < currentNode.degree[i]))
                        i = j; // branch sur le plus petit degrée supérieur à 2
                }
                if (i < 0) {
                    if (currentNode.lowerBound < bestNode.lowerBound) {
                        bestNode = currentNode;
                        if (verbose) System.out.println("upper bound : " + bestNode.lowerBound);
                        if (verbose) System.err.printf("%f", bestNode.lowerBound);
                    }
                    break;
                }
                if (verbose) System.err.printf(".");
                PriorityQueue<Node> children = new PriorityQueue<Node>(11, new NodeComparator());
                children.add(exclude(currentNode, i, currentNode.parent[i])); // branching exclude each node
                for (int j = 0; j < n; j++) {
                    if (currentNode.parent[j] == i) children.add(exclude(currentNode, i, j));
                }
                currentNode = children.poll();
                pq.addAll(children); // optimised ?
            } while (currentNode.lowerBound < bestNode.lowerBound); // mix DFS / BFS
            if (verbose) System.err.printf("%n");
            currentNode = pq.poll();
        } while (currentNode != null && currentNode.lowerBound < bestNode.lowerBound);
        // output suitable for gnuplot
        // set style data vector
        if (verbose) System.out.println("---------------------------");
        if (verbose) System.out.println("-> optimum found : " + bestNode.lowerBound);
        if (verbose) System.out.println("-> " + visited_node + " nodes visisted");
        int j = 0;
        do {
            int i = bestNode.parent[j];
            //System.out.printf("%f\t%f\t%f\t%f%n", x[j], y[j], x[i] - x[j], y[i] - y[j]);
            if (verbose) System.out.printf("(%d,%d) ", j, i);
            j = i;
        } while (j != 0);
    }

    /**
     * Exclude a edge from the Solver.TSP
     *
     * @param node
     * @param i
     * @param j
     * @return
     */
    private Node exclude(Node node, int i, int j) {
        Node child = new Node();
        child.excluded = node.excluded.clone();
        child.excluded[i] = node.excluded[i].clone();
        child.excluded[j] = node.excluded[j].clone();
        child.excluded[i][j] = true;
        child.excluded[j][i] = true;
        computeHeldKarp(child);
        return child;
    }

    private void computeHeldKarp(Node node) {
        node.pi = new double[n];
        node.lowerBound = Double.MIN_VALUE;
        node.degree = new int[n];
        node.parent = new int[n];
        double lambda = 0.1;
        double eps = 1e-06;
        while (lambda > eps) {
            double previousLowerBound = node.lowerBound;
            computeOneTree(node);
            if (!(node.lowerBound < bestNode.lowerBound)) return; // retour direct si meilleure LB TODO : impact
            if (!(node.lowerBound < previousLowerBound)) lambda *= 0.9;
            int denom = 0; // distance from the Solver.TSP
            for (int i = 1; i < n; i++) {
                int d = node.degree[i] - 2;
                denom += d * d;
            }
            if (denom == 0) return;
            double t = lambda * node.lowerBound / denom; // TODO: impact de cette step
            for (int i = 1; i < n; i++) node.pi[i] += t * (node.degree[i] - 2);
        }
    }

    private void computeOneTree(Node node) {
        // compute adjusted costs
        node.lowerBound = 0.0;
        Arrays.fill(node.degree, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                costWithPi[i][j] = node.excluded[i][j] ? Double.MAX_VALUE : cost[i][j] + node.pi[i] + node.pi[j];
        }
        int firstNeighbor;
        int secondNeighbor;
        // find the two cheapest edges from 0
        if (costWithPi[0][2] < costWithPi[0][1]) {
            firstNeighbor = 2;
            secondNeighbor = 1;
        } else {
            firstNeighbor = 1;
            secondNeighbor = 2;
        }
        for (int j = 3; j < n; j++) {
            if (costWithPi[0][j] < costWithPi[0][secondNeighbor]) {
                if (costWithPi[0][j] < costWithPi[0][firstNeighbor]) {
                    secondNeighbor = firstNeighbor;
                    firstNeighbor = j;
                } else {
                    secondNeighbor = j;
                }
            }
        }
        addEdge(node, 0, firstNeighbor);
        Arrays.fill(node.parent, firstNeighbor);
        node.parent[firstNeighbor] = 0;
        // compute the minimum spanning tree on nodes 1..n-1
        double[] minCost = costWithPi[firstNeighbor].clone();
        for (int k = 2; k < n; k++) {
            int i;
            for (i = 1; i < n; i++) {
                if (node.degree[i] == 0) break;
            }
            for (int j = i + 1; j < n; j++) {
                if (node.degree[j] == 0 && minCost[j] < minCost[i]) i = j;
            }
            addEdge(node, node.parent[i], i);
            for (int j = 1; j < n; j++) {
                if (node.degree[j] == 0 && costWithPi[i][j] < minCost[j]) {
                    minCost[j] = costWithPi[i][j];
                    node.parent[j] = i;
                }
            }
        }
        addEdge(node, 0, secondNeighbor);
        node.parent[0] = secondNeighbor;
        //node.lowerBound = Math.rint(node.lowerBound); // TODO: 2nd round ?
        //node.lowerBound = node.lowerBound;
    }

    private void addEdge(Node node, int i, int j) {
        double q = node.lowerBound;
        node.lowerBound += costWithPi[i][j];
        node.degree[i]++;
        node.degree[j]++;
    }

    // getters

    public double getLB() {
        return bestNode.lowerBound;
    }
}

class Node {
    public boolean[][] excluded;
    // Held--Karp solution
    public double[] pi;
    public double lowerBound;
    public int[] degree;
    public int[] parent;
}

class NodeComparator implements Comparator<Node> {
    public int compare(Node a, Node b) {
        return Double.compare(a.lowerBound, b.lowerBound);
    }
}

