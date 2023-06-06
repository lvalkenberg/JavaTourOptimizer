package core;

import heuristics.Heuristic;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.*;
import java.sql.SQLOutput;
import java.util.*;
import java.util.regex.*;
import java.util.function.Consumer;

/**
 * Inspired from : https://stackoverflow.com/questions/7159259/optimized-tsp-algorithms/7163961#7163961
 * for the base of the solver.
 */
public class TSP {
    public boolean verbose = false;
    // best solution found so far
    private Node bestNode;
    // instance name
    public String name;
    // number of cities
    private int n;
    // city locations
    private double[] x;
    private double[] y;
    // cost (or distance) matrix
    private double[][] cost;
    // matrix of adjusted costs considering the nodes' potential
    private double[][] costWithPi;
    // matrix of cost without the 1Tree root
    private double[][] subGraph;
    private double globalCostWithPi[][];
    // best pi found in HK
    private double[] bestHKPi;

    // counter on the visited node
    public long visitedNodeCounter;
    // counter on the number of inconsistent node
    public int inconsistentCounter = 0;
    public int iterTillImprovement;
    public int uselessSubgradientCounter;
    // optimality gap
    public double optimalityGap = 1;

    // data structure storing the unvisited nodes
    private Stack<Node> stackNodes = new Stack();
    private PriorityQueue<Node> pqNodes = new PriorityQueue<Node>(11, new NodeComparator());
    private PriorityQueue<Node> pqChildren = new PriorityQueue<Node>(11, new NodeComparator());

    // mst algo
    public PrimCCT mst;

    // last conflict
    private int lastConflictNode;
    private int getLastConflictNode2;
    // true if the currentNode is a direct children of the last one
    private boolean directChildren = false;

    private static ArrayList<String> supportedSearchStrategy = new ArrayList<>(Arrays.asList("DFS", "BFS", "DFS on BFS"));
    private static ArrayList<String> supportedBranchStrategy = new ArrayList<>(Arrays.asList("exclude", "nearest neighbour", "maxCost", "minCost path", "minCost"));

    // default solver parameters
    public boolean margCostFilter = true;                   // marginal cost filtering
    public boolean repCostFilter = true;                    // replacement cost filtering
    public String branchingStrategy = "maxCost";       // branching strategy available : nearest neighbour, exclude, maxCost, minCost & minCost path
    public String searchStrategy = "DFS on BFS";            // search strategy available : BFS, DFS, DFS on BFS
    public boolean lastConflictSearch = true;               // last conflict search heuristic
    public boolean incrementalPi = true;                    // hot restart of the Lagrangian multipliers
    public boolean fixedPi = false;                         // fixed multiplier
    public boolean keepBestPi = false;                      // keep the best lower bound found in the subgradient methode
    public boolean initialHeuristic = true;                // initial heuristic
    public double timeout = -1;                             // timeout in second, None by default
    private long startTime = System.currentTimeMillis();

    public Consumer<int[]> onBetterSolutionFound;    // callback function on better solution found

    public TSP(double[][] distanceMatrix) {
        this.cost = distanceMatrix;
        this.n = distanceMatrix.length;
    }

    public TSP(String tspFile) throws IOException {
        readInstance(tspFile);
    }

    public TSP(TSPBuilder builder, String tspFile) throws IOException {
        this.margCostFilter = builder.margCostFilter;
        this.repCostFilter = builder.repCostFilter;
        this.branchingStrategy = builder.branchingStrategy;
        this.searchStrategy = builder.searchStrategy;
        this.lastConflictSearch = builder.lastConflictSearch;
        this.incrementalPi = builder.incrementalPi;
        this.fixedPi = builder.fixedPi;
        this.keepBestPi = builder.keepBestPi;
        this.initialHeuristic = builder.initialHeuristic;
        this.onBetterSolutionFound = builder.onBetterSolutionFound;
        readInstance(tspFile);
    }

    public static class TSPBuilder {
        private boolean margCostFilter;
        private boolean repCostFilter;
        private String branchingStrategy;
        private String searchStrategy;
        private boolean lastConflictSearch;
        private boolean incrementalPi;
        private boolean fixedPi;
        private boolean keepBestPi;
        private boolean initialHeuristic;
        private double timeout;
        private Consumer<int[]> onBetterSolutionFound = null;

        public TSPBuilder() {
            // Initialize default values
            this.margCostFilter = true;
            this.repCostFilter = true;
            this.branchingStrategy = "nearest neighbour";
            this.searchStrategy = "DFS on BFS";
            this.lastConflictSearch = true;
            this.incrementalPi = true;
            this.fixedPi = false;
            this.keepBestPi = true;
            this.initialHeuristic = true;
            this.timeout = -1;
        }

        public TSPBuilder setMargCostFilter(boolean margCostFilter) {
            this.margCostFilter = margCostFilter;
            return this;
        }

        public TSPBuilder setRepCostFilter(boolean repCostFilter) {
            this.repCostFilter = repCostFilter;
            return this;
        }

        public TSPBuilder setBranchingStrategy(String branchingStrategy) {
            if (!supportedBranchStrategy.contains(branchingStrategy)) {
                throw new IllegalArgumentException("Unsupported branching strategy. Use : " + String.join(", ", supportedBranchStrategy));
            }
            this.branchingStrategy = branchingStrategy;
            return this;
        }

        public TSPBuilder setSearchStrategy(String searchStrategy) {
            if (!supportedSearchStrategy.contains(searchStrategy)) {
                throw new IllegalArgumentException("Unsupported search strategy. Use : " + String.join(", ", supportedSearchStrategy));
            }
            this.searchStrategy = searchStrategy;
            return this;
        }

        public TSPBuilder setLastConflictSearch(boolean lastConflictSearch) {
            this.lastConflictSearch = lastConflictSearch;
            return this;
        }

        public TSPBuilder setIncrementalPi(boolean incrementalPi) {
            this.incrementalPi = incrementalPi;
            return this;
        }

        public TSPBuilder setFixedPi(boolean fixedPi) {
            this.fixedPi = fixedPi;
            return this;
        }

        public TSPBuilder setKeepBestPi(boolean keepBestPi) {
            this.keepBestPi = keepBestPi;
            return this;
        }

        public TSPBuilder setInitialHeuristic(boolean initialHeuristic) {
            this.initialHeuristic = initialHeuristic;
            return this;
        }

        public TSPBuilder setTimeout(double timeout) {
            this.timeout = timeout;
            return this;
        }

        public TSPBuilder setOnBetterSolutionFound(Consumer<int[]> onBetterSolutionFound) {
            this.onBetterSolutionFound = onBetterSolutionFound;
            return this;
        }

        public TSP build(String tspFile) throws IOException {
            return new TSP(this, tspFile);
        }
    }

    /**
     * Given the input file path use the adapted function to open it
     *
     * @param filePath file path of the instance (.xml or .tsp EUC_2D format)
     * @throws IOException Error with the provided path
     */
    private void readInstance(String filePath) throws IOException {
        if (filePath.toLowerCase().endsWith(".xml")) xmlReader(filePath);
        else if (filePath.toLowerCase().endsWith(".tsp")) readInput(filePath);
        else throw new IllegalArgumentException("File extension not supported");
    }

    /**
     * EUC_2D .tsp reader (http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp95.pdf)
     *
     * @param filePath complete path of the .tsp file
     * @throws IOException Error with the provided path
     */
    public void readInput(String filePath) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(filePath));
        Pattern specification = Pattern.compile("\\s*([A-Z_]+)\\s*(:\\s*([A-Z_0-9]+))?\\s*");
        Pattern data = Pattern.compile("\\s*([0-9]+)\\s+([-+.0-9Ee]+)\\s+([-+.0-9Ee]+)\\s*");
        String line;
        while ((line = in.readLine()) != null) {
            Matcher m = specification.matcher(line);
            if (!m.matches()) continue;
            String keyword = m.group(1);
            if (keyword.equals("DIMENSION")) {
                n = Integer.parseInt(m.group(3));
                cost = new double[n][n];
            } else if (keyword.equals("EDGE_WEIGHT_TYPE")) {
                if (!m.group(3).equals("EUC_2D"))
                    throw new IllegalArgumentException("Only TSP file with EDGE_WEIGHT_TYPE = EUC_2D are supported");
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
                        cost[i][j] = Math.rint(Math.sqrt(dx * dx + dy * dy));
                        //cost[i][j] = Math.round(Math.sqrt(dx * dx + dy * dy) * 1e7) / 1e7;
                    }
                }
            }
        }
    }

    /**
     * .xml TSP reader (http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/XML-TSPLIB/Description.pdf)
     *
     * @param xmlPath path to the .xml file
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

    /**
     * Solve the instance.
     *
     * @return a Integer array with the node in the order of the TSP tour. null if no tour found
     */
    public int[] solve() {
        mst = new PrimCCT(n - 1);
        bestNode = new Node(n);
        bestNode.lowerBound = Double.MAX_VALUE;
        if (initialHeuristic) {
            int[] tour = Heuristic.solve(this.cost);
            this.bestNode.lowerBound = Heuristic.tourCost(tour, this.cost);
            if (verbose) System.out.println("initial solution = " + this.bestNode.lowerBound);
            for (int i = 0; i < n; i++) {
                bestNode.V[i] = tour[i];
                bestNode.W[i] = tour[(i + 1) % n];
            }
            //this.bestNode.lowerBound += 1;
        }
        // initialization of the search tree
        Node currentNode = new Node(n);
        currentNode.excluded = new boolean[n][n];
        currentNode.mandatoryEdges = new int[n][2];
        for (int i = 0; i < n; i++) {
            currentNode.mandatoryEdges[i][0] = -1;
            currentNode.mandatoryEdges[i][1] = -1;
        }
        costWithPi = new double[n][n];
        subGraph = new double[n - 1][n - 1];
        visitedNodeCounter = 0;
        if (fixedPi) globalCostWithPi = new double[n][];
        computeHeldKarp(currentNode);
        do {
            if (fixedPi) visitedNodeCounter++;
            // timeout
            if (timeout != -1 && (System.currentTimeMillis() - startTime) / 1000.0 > timeout) {
                return currentBestTour();
            }
            int i = -1;
            for (int j = 0; j < n; j++) {
                if (currentNode.degree[j] != 2) {
                    i = j;
                    break;
                }
            }

            if (i < 0 && currentNode.lowerBound < bestNode.lowerBound) {
                bestNode = currentNode;
                if (verbose) System.err.printf("%f \n", bestNode.lowerBound);
                if (onBetterSolutionFound != null)
                    new Thread(() -> this.onBetterSolutionFound.accept(currentBestTour())).start();
            }
            PriorityQueue<Node> children = branch(currentNode);
            currentNode = nextNode(children);
            if (!pqNodes.isEmpty())
                optimalityGap = Math.min(optimalityGap, (bestNode.lowerBound - pqNodes.peek().lowerBound) / bestNode.lowerBound);
            //System.out.println(optimalityGap);
        } while (currentNode != null && currentNode.lowerBound < bestNode.lowerBound);

        if (verbose) {
            System.out.println("-> optimum found : " + bestNode.lowerBound);
            System.out.println("-> " + visitedNodeCounter + " nodes visisted");
        }
        return currentBestTour();
    }


    /**
     * Add the node to the priority queue if the node is not pruned or inconsistent
     *
     * @param pq   priority queue
     * @param node node to add in the priority queue
     */
    private void addPQ(PriorityQueue<Node> pq, Node node) {
        if (!(node.inconsistency) && node.lowerBound < bestNode.lowerBound)
            pq.add(node);
        else if (lastConflictSearch) {
            lastConflictNode = node.branchedVertices;
        }
    }

    /**
     * Branch following the branching strategy chose.
     *
     * @param currentNode node to branch on
     * @return A priority queue containing children Node or NULL if the currentnNode inconsistent
     */
    private PriorityQueue<Node> branch(Node currentNode) {
        pqChildren.clear();
        if (this.branchingStrategy == "exclude") {
            int i = -1;
            for (int j = 0; j < n; j++) { // Exclude each edge adjacent to a node with an degree higher than 2
                if (currentNode.degree[j] > 2 && (i < 0 || currentNode.degree[j] < currentNode.degree[i])) i = j;
            }
            // exclude each edge of the 1-tree adjacent to node i
            for (int j = 0; j < n; j++) {
                if (currentNode.V[j] == i) addPQ(pqChildren, exclude(currentNode, i, currentNode.W[j]));
                else if (currentNode.W[j] == i) addPQ(pqChildren, exclude(currentNode, i, currentNode.V[j]));
            }
        } else if (this.branchingStrategy == "maxCost") { // force higher weight edge adjacent to a node with a degree < 2
            int i = -1;
            int j = -1;
            // branch on the last conflict variable
            if (lastConflictSearch && currentNode.degree[lastConflictNode] < 2 && stillNonExcludeEdge(currentNode, lastConflictNode)) {
                i = lastConflictNode;
                for (int k = 0; k < n; k++) {
                    if (i != k && !currentNode.excluded[i][k] && k != currentNode.mandatoryEdges[i][0] && currentNode.mandatoryEdges[k][1] == -1 && currentNode.mandatoryEdges[i][1] == -1 && (j == -1 || cost[i][j] < cost[i][k])) {
                        j = k;
                    }
                }
            } else { // force first node on the min weight edge
                for (int l = 0; l < n; l++) {
                    if (currentNode.degree[l] < 2 && stillNonExcludeEdge(currentNode, l)) {
                        for (int k = 0; k < n; k++) {
                            if (l != k && !currentNode.excluded[l][k] && k != currentNode.mandatoryEdges[l][0] && currentNode.mandatoryEdges[k][1] == -1 && currentNode.mandatoryEdges[l][1] == -1 && (j == -1 || cost[i][j] < cost[l][k])) {
                                j = k;
                                i = l;
                            }
                        }
                    }
                }
            }
            if (i == -1) return pqChildren; // TSP solution
            if (j == -1) throw new IllegalArgumentException();

            addPQ(pqChildren, exclude(currentNode, i, j));
            addPQ(pqChildren, include(currentNode, i, j));
        } else if (this.branchingStrategy == "minCost") { // force higher weight edge adjacent to a node with a degree < 2
            int i = -1;
            int j = -1;
            // branch on the last conflict variable
            if (lastConflictSearch && currentNode.degree[lastConflictNode] < 2 && stillNonExcludeEdge(currentNode, lastConflictNode)) {
                i = lastConflictNode;
                for (int k = 0; k < n; k++) {
                    if (i != k && !currentNode.excluded[i][k] && k != currentNode.mandatoryEdges[i][0] && currentNode.mandatoryEdges[k][1] == -1 && currentNode.mandatoryEdges[i][1] == -1 && (j == -1 || cost[i][j] > cost[i][k])) {
                        j = k;
                    }
                }
            } else { // force first node on the min weight edge
                for (int l = 0; l < n; l++) {
                    if (currentNode.degree[l] < 2 && stillNonExcludeEdge(currentNode, l)) {
                        for (int k = 0; k < n; k++) {
                            if (l != k && !currentNode.excluded[l][k] && k != currentNode.mandatoryEdges[l][0] && currentNode.mandatoryEdges[k][1] == -1 && currentNode.mandatoryEdges[l][1] == -1 && (j == -1 || cost[i][j] > cost[l][k])) {
                                j = k;
                                i = l;
                            }
                        }
                    }
                }
            }
            if (i == -1) return pqChildren; // TSP solution
            if (j == -1) throw new IllegalArgumentException();

            addPQ(pqChildren, exclude(currentNode, i, j));
            addPQ(pqChildren, include(currentNode, i, j));
        } else if (this.branchingStrategy == "minCost path") { // force higher weight edge adjacent to a node with a degree < 2
            int i = -1;
            int j = -1;
            for (int k = 0; k < n; k++) { // find a path
                if (currentNode.mandatoryEdges[k][1] != -1) {
                    i = k;
                    break;
                }
            }
            if (i == -1) { //no path found -> create one
                for (int k = 0; k < n; k++) { // find a path
                    if (currentNode.mandatoryEdges[k][0] != -1) {
                        i = k;
                        break;
                    }
                }
                if (i == -1) {
                    //return fallbackBranching(currentNode);
                    i = 0;
                }
            }

            int prev = -1;
            int tmp;
            for (int k = 0; k <= n; k++) { // travel path
                tmp = i;
                if (currentNode.mandatoryEdges[i][1] == -1) break;
                if (prev != currentNode.mandatoryEdges[i][1]) i = currentNode.mandatoryEdges[i][1];
                else i = currentNode.mandatoryEdges[i][0];
                prev = tmp;
                if (k == n) return null; // subtour
            }
            for (int k = 0; k < n; k++) { // new node
                if (i != k && !currentNode.excluded[i][k] && currentNode.mandatoryEdges[i][0] != k && currentNode.mandatoryEdges[k][1] == -1 && (j == -1 || costWithPi[i][k] < costWithPi[i][j])) {
                    j = k;
                }
            }
            if (j == -1) return null; // inconsitency
            addPQ(pqChildren, exclude(currentNode, i, j));
            addPQ(pqChildren, include(currentNode, i, j));
        } else { // default : nearest neighbour = force the lightest weight edge adjacent to the first node with degree less than 2
            int i = -1;
            if (lastConflictSearch && currentNode.degree[lastConflictNode] < 2 && stillNonExcludeEdge(currentNode, lastConflictNode))
                i = lastConflictNode;
            else {
                for (int j = 0; j < n; j++) {
                    if (currentNode.degree[j] < 2 && stillNonExcludeEdge(currentNode, j)) {
                        i = j;
                        break;
                    }
                }
            }
            if (i == -1) {
                return null; // inconsistency
            }

            int k = -1;
            for (int j = 0; j < n; j++) {
                if (i != j && !currentNode.excluded[i][j] && j != currentNode.mandatoryEdges[i][0] && currentNode.mandatoryEdges[j][1] == -1 && currentNode.mandatoryEdges[i][1] == -1 && (k == -1 || cost[i][k] > cost[i][j]))
                    k = j;
            }
            if (k == -1) {
                throw new IllegalArgumentException();   // should not append
            }
            addPQ(pqChildren, exclude(currentNode, i, k));
            addPQ(pqChildren, include(currentNode, i, k));
        }
        return pqChildren;
    }


    /**
     * Return the number of edges adjacent to i that can still be exclude.
     *
     * @param node current search tree node
     * @param i    graph node
     * @return the number of edge adjacent to i than can still be excluded
     */
    private boolean stillNonExcludeEdge(Node node, int i) {
        int sum = 0;
        for (int j = 0; j < n; j++) {
            if (!node.excluded[i][j] && i != j && node.mandatoryEdges[j][1] == -1) sum++;
        }
        return sum >= 2;
    }

    /**
     * Return the next node to explore depending of the search strategy.
     *
     * @param children priority queue containing children of the node that we just explored
     * @return the next node to explore, or null if there is no interesting node left
     */
    private Node nextNode(PriorityQueue<Node> children) {
        Node child = null;
        if (searchStrategy == "DFS") {   // Depth-First Search
            if (children != null) {
                child = children.poll();
                directChildren = true;
                for (Node c : children) {
                    if (c.lowerBound < bestNode.lowerBound && !child.inconsistency) {
                        stackNodes.add(c);
                        pqNodes.add(c);
                    }
                }
            }
            while ((child == null || child.lowerBound >= bestNode.lowerBound) && !stackNodes.isEmpty()) {
                child = stackNodes.pop();
                pqNodes.remove(child); //O(size queue)
                directChildren = false;
            }
            if ((child != null && child.lowerBound >= bestNode.lowerBound)) return null;
            else return child;
        } else if (searchStrategy == "BFS") {  // Best-First Search
            if (children != null) pqNodes.addAll(children);
            return pqNodes.poll();
        } else {   // BFS mixed with DFS
            if (children != null) {
                child = children.poll();
                directChildren = true;
                pqNodes.addAll(children);
            }
            if (child == null || child.lowerBound >= bestNode.lowerBound) {
                child = pqNodes.poll();
                directChildren = false;
            }
            return child;
        }
    }

    /**
     * Exclude a edge (i,j) from the node TSP tour
     *
     * @param node parent node of the search tree
     * @param i    first end point of the edge
     * @param j    second end point of the edge
     * @return the children node of @node with (i,j) as exclude edge
     */
    private Node exclude(Node node, int i, int j) {
        if (node.mandatoryEdges[i][0] == j || node.mandatoryEdges[i][1] == j) node.inconsistency = true;
        Node child = new Node(node);
        child.exclude(i, j);

        if (!fixedPi) computeHeldKarp(child);
        else computeOneTree(child);

        child.branchedVertices = i;

        if (!child.inconsistency) reduce(child);
        return child;
    }


    /**
     * Force the edge (i,j) in the node's TSP tour
     *
     * @param node parent node of the search tree
     * @param i    first end point of the edge
     * @param j    second end point of the edge
     * @return the children node of @node with (i,j) as mandatory edge
     */
    private Node include(Node node, int i, int j) {
        if (node.excluded[i][j]) throw new IllegalArgumentException();
        Node child = new Node(node);
        child.forceEdge(i, j);

        if (!fixedPi) computeHeldKarp(child);
        else computeOneTree(child);

        child.branchedVertices = i;

        if (!child.inconsistency) reduce(child);
        return child;
    }

    /**
     * Use techniques to reduce the problem size.
     *
     * @param node current node
     */
    private void reduce(Node node) {
        if (bestNode.lowerBound != Double.MAX_VALUE && node.lowerBound <= bestNode.lowerBound) {
            if (margCostFilter) marginalCostFiltering(node);
            if (repCostFilter) replacementCostForcing(node);
        }
    }

    /**
     * Compute the Held-Karp relaxation
     *
     * @param node node of the search tree
     */
    private void computeHeldKarp(Node node) {
        visitedNodeCounter++;
        double lambda = 0.1;
        if (!incrementalPi) {
            for (int i = 0; i < n; i++) node.pi[i] = 0;
        }
        double bestPreviousLowerBound = Double.MIN_VALUE;
        boolean improvedPreviousPi = false;
        int improvementCounter = 0;
        int finalValue = -1;
        double eps = 1e-06;
        while (lambda > eps) {
            double previousLowerBound = node.lowerBound;
            computeOneTree(node);
            improvementCounter++;
            if (node.lowerBound > bestPreviousLowerBound) {
                if (bestPreviousLowerBound != Double.MIN_VALUE) improvedPreviousPi = true;
                if (improvedPreviousPi && finalValue == -1) finalValue = improvementCounter;
                bestPreviousLowerBound = node.lowerBound;
                bestHKPi = Arrays.copyOf(node.pi, n);
            }
            // reset the lagrangian multiplier if the lower bound is negatif
            if (node.lowerBound < 0) {
                for (int i = 0; i < n; i++) {
                    node.pi[i] = 0;
                }
                continue;
            }
            if (node.inconsistency) {
                inconsistentCounter++;
                return;
            }
            if (!(node.lowerBound < bestNode.lowerBound)) return; // node to prune
            if (!(node.lowerBound < previousLowerBound)) lambda *= 0.9;
            if (lambda <= eps) continue;

            int denom = 0; // distance from the a TSP tour
            for (int i = 1; i < n; i++) {
                int d = node.degree[i] - 2;
                denom += d * d;
            }
            if (denom == 0) return; // TSP tour found
            double t = lambda * bestPreviousLowerBound / denom; // step
            for (int i = 1; i < n; i++) node.pi[i] += t * (node.degree[i] - 2);
        }
        if (keepBestPi && node.lowerBound < bestPreviousLowerBound) {
            node.pi = Arrays.copyOf(bestHKPi, n);
            computeOneTree(node);
        }
        if (fixedPi && visitedNodeCounter <= 1) {
            if (fixedPi && visitedNodeCounter <= 1) {
                for (int i = 0; i < costWithPi.length; i++) {
                    globalCostWithPi[i] = costWithPi[i].clone();
                }
            }

            if (directChildren && !improvedPreviousPi && node.lowerBound <= bestPreviousLowerBound) {
                uselessSubgradientCounter++;
            }
        }
    }

    /**
     * Compute the 1-tree given the nodes' potential.
     *
     * @param node node of the search tree for which the 1-tree is calculated
     */
    private void computeOneTree(Node node) {
        // compute adjusted costs
        node.lowerBound = 0.0;
        Arrays.fill(node.degree, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                if (fixedPi && visitedNodeCounter > 1)
                    costWithPi[i][j] = node.excluded[i][j] ? Double.MAX_VALUE : globalCostWithPi[i][j];
                else costWithPi[i][j] = node.excluded[i][j] ? Double.MAX_VALUE : cost[i][j] + node.pi[i] + node.pi[j];
        }
        // root edges
        int firstNeighbor = -1;
        int secondNeighbor = -1;

        if (node.mandatoryEdges[0][1] == -1) {
            for (int i = 1; i < n; i++) {
                if (node.mandatoryEdges[i][1] == -1) {
                    if (firstNeighbor == -1 || costWithPi[0][i] < costWithPi[0][firstNeighbor]) {
                        secondNeighbor = firstNeighbor;
                        firstNeighbor = i;
                    } else if (secondNeighbor == -1 || costWithPi[0][i] < costWithPi[0][secondNeighbor]) {
                        secondNeighbor = i;
                    }
                }
            }
        }
        if (node.mandatoryEdges[0][1] != -1) {
            firstNeighbor = node.mandatoryEdges[0][0];
            secondNeighbor = node.mandatoryEdges[0][1];
        } else if (node.mandatoryEdges[0][0] != -1 && node.mandatoryEdges[0][0] != firstNeighbor) {
            secondNeighbor = node.mandatoryEdges[0][0];
        }

        if (firstNeighbor == -1 || secondNeighbor == -1) {
            node.inconsistency = true;
            return;
        }

        // compute the 1-tree
        for (int i = 0; i < n - 1; i++) {
            for (int j = 0; j < n - 1; j++) {
                subGraph[i][j] = costWithPi[i + 1][j + 1];
            }
        }

        addEdge(node, 0, firstNeighbor);
        addEdge(node, 0, secondNeighbor);

        // Case of MaxDegTwoMST
//        int[] rootEdges = {firstNeighbor, secondNeighbor};
//        mst.computeMST(subGraph, node.mandatoryEdges, rootEdges);
        // Classical MST
        mst.computeMST(subGraph, node.mandatoryEdges);
        node.inconsistency = mst.inconsistent;
        if (node.inconsistency) return;

        for (int i = 1; i < n - 1; i++) {
            // shift edges from the MST to consider the root of the 1tree
            int v = i + 1;
            int w = mst.parent[i] + 1;
            node.V[i] = v;
            node.W[i] = w;
            addEdge(node, v, w);
        }
        node.V[0] = 0;
        node.W[0] = firstNeighbor;
        node.V[n - 1] = 0;
        node.W[n - 1] = secondNeighbor;
    }

    /**
     * Add an edge to the minimum 1-tree
     *
     * @param node current node
     * @param i    first end point of the edge
     * @param j    second end point of the edge
     */
    private void addEdge(Node node, int i, int j) {
        node.lowerBound += costWithPi[i][j];
        node.degree[i]++;
        node.degree[j]++;
    }

    /**
     * Apply marginal cost filtering.
     *
     * @param node to apply the filtering
     */
    private void marginalCostFiltering(Node node) {
        // filtering on the min 1-tree root
        double suppWeight = Math.max(costWithPi[0][node.W[0]], costWithPi[0][node.W[n - 1]]); // max edge weight of the 2 root edges
        Edge toExclude = new Edge(-1, -1, -1);
        for (int i = 0; i < n; i++) {
            if (!node.excluded[0][i] && costWithPi[0][i] - suppWeight + node.lowerBound > bestNode.lowerBound) {
                toExclude.v = 0;
                toExclude.w = i;
                toExclude.weight = costWithPi[0][i];
                node.exclude(toExclude);
            }
        }
        // filtering on the MST
        HashMap<Edge, ArrayList<Edge>> supports = mst.computeAllSupportEdges();
        for (Map.Entry<Edge, ArrayList<Edge>> entry : supports.entrySet()) {
            for (Edge nonTreeEdge : entry.getValue()) {
                if (nonTreeEdge.weight - entry.getKey().weight + node.lowerBound >= bestNode.lowerBound) {
                    toExclude.v = nonTreeEdge.v + 1;
                    toExclude.w = nonTreeEdge.w + 1;
                    toExclude.weight = nonTreeEdge.weight;
                    node.exclude(toExclude);
                }
            }
        }
    }

    /**
     * Force the mandatory edges based on their replacement cost.
     * If an edge is mandatory but one of his vertices as already
     * 2 adjacent mandatory edges, the node is inconsistent.
     *
     * @param node to apply the filtering
     */
    private void replacementCostForcing(Node node) {
        // filtering on the min 1-tree root
        boolean newOT = false; // does the MST will change with the addition of mandatory edges
        double minRepWeight = Double.MAX_VALUE;
        for (int i = 1; i < n; i++) {
            if (!node.excluded[0][i] && node.mandatoryEdges[i][1] == -1 && node.W[0] != i && node.W[n - 1] != i && costWithPi[0][i] < minRepWeight) {
                minRepWeight = costWithPi[0][i];
            }
        }
        if (node.lowerBound - costWithPi[0][node.W[0]] + minRepWeight > bestNode.lowerBound) {
            node.forceEdge(0, node.W[0]);
            newOT = true;
        }
        if (node.lowerBound - costWithPi[0][node.W[n - 1]] + minRepWeight > bestNode.lowerBound) {
            node.forceEdge(0, node.W[n - 1]);
            newOT = true;
        }
        // filtering on the MST edges
        double[] repCost = mst.computeReplacementCost();
        for (int i = 0; i < repCost.length; i++) {
            if (node.lowerBound + repCost[i] > bestNode.lowerBound) {
                int v = i + 1;
                int w = mst.parent[i] + 1;
                if (node.mandatoryEdges[v][1] != -1 && (!node.isMandatoryEdge(v, w) || !node.isMandatoryEdge(w, v))) {
                    node.inconsistency = true; // the edge has already a degree 2
                    return;
                }
                if (node.mandatoryEdges[w][1] != -1 && (!node.isMandatoryEdge(v, w) || !node.isMandatoryEdge(w, v))) {
                    node.inconsistency = true;
                    return;
                }
                if (node.mandatoryEdges[v][1] != -1 || node.mandatoryEdges[w][1] != -1) newOT = true;
                node.forceEdge(v, w);
            }
        }
        if (newOT) computeOneTree(node);
    }

    /**
     * Return the current best tour found as a sequence of points to visit
     *
     * @return a int[] in the order of the tour, null if no tour found
     */
    public int[] currentBestTour() {
        if (bestNode.lowerBound == Double.MAX_VALUE) return null;
        HashMap<Integer, ArrayList<Integer>> tourDict = new HashMap<>();
        boolean[] visited = new boolean[n];
        for (int i = 0; i < n; i++) {
            if (!tourDict.containsKey(bestNode.V[i])) tourDict.put(bestNode.V[i], new ArrayList<>());
            if (!tourDict.containsKey(bestNode.W[i])) tourDict.put(bestNode.W[i], new ArrayList<>());
            tourDict.get(bestNode.V[i]).add(bestNode.W[i]);
            tourDict.get(bestNode.W[i]).add(bestNode.V[i]);
        }
        int[] tour = new int[n];
        tour[0] = 0;
        visited[0] = true;
        for (int i = 1; i < n; i++) {
            for (int next : tourDict.get(tour[i - 1])) {
                if (visited[next]) continue;
                tour[i] = next;
                visited[next] = true;
                break;
            }
        }
        return tour;
    }

    /**
     * Lower bound getter
     *
     * @return the cost of the best lower tour
     */
    public double getLowerBound() {
        return bestNode.lowerBound;
    }

}

class Node {
    // excluded edges
    public boolean[][] excluded;
    public boolean[] clonedExclu;
    public int[] nonEclAdjEdge;

    // mandatory edges
    public int[][] mandatoryEdges;
    public boolean[] clonedMandatory;

    // Held--Karp solution
    public double[] pi;
    public double lowerBound;
    public int[] degree;

    // edges in the 1tree
    public int[] V;
    public int[] W;

    boolean inconsistency;

    public int branchedVertices;

    /**
     * Initialise the parameter of a node
     *
     * @param n number of cities
     */
    public Node(int n) {
        this.lowerBound = Double.MIN_VALUE;
        this.degree = new int[n];
        this.pi = new double[n];
        //node.parent = new int[n];
        this.V = new int[n];
        this.W = new int[n];
        this.nonEclAdjEdge = new int[n];
        for (int i = 0; i < n; i++) {
            nonEclAdjEdge[i] = n;
        }
    }

    /**
     * Initialise a node based on an other node (paren node)
     *
     * @param node parent node
     */
    public Node(Node node) {
        int n = node.pi.length;
        this.lowerBound = Double.MIN_VALUE;
        this.degree = new int[n];
        this.pi = Arrays.copyOf(node.pi, n);
        //node.parent = new int[n];
        this.V = new int[n];
        this.W = new int[n];
        this.nonEclAdjEdge = new int[n];
        this.nonEclAdjEdge = Arrays.copyOf(node.nonEclAdjEdge, n);
        this.excluded = node.excluded.clone();
        this.clonedExclu = new boolean[n];
        this.mandatoryEdges = node.mandatoryEdges.clone();
        this.clonedMandatory = new boolean[n];
    }

    /**
     * Exclude a edge from the Solver.TSP in the node.
     *
     * @param edge edge that will be exclude
     */
    public void exclude(Edge edge) {
        int i = edge.v;
        int j = edge.w;

        nonEclAdjEdge[i]--;
        nonEclAdjEdge[j]--;
        if (nonEclAdjEdge[i] < 2 || nonEclAdjEdge[j] < 2) {
            this.inconsistency = true;
        }

        if (!clonedExclu[i]) {
            this.excluded[i] = this.excluded[i].clone();
            clonedExclu[i] = true;
        }
        if (!clonedExclu[j]) {
            this.excluded[j] = this.excluded[j].clone();
            clonedExclu[j] = true;
        }
        this.excluded[i][j] = true;
        this.excluded[j][i] = true;
    }

    /**
     * Exclude a edge from the Solver.TSP in the node.
     *
     * @param i first edge endpoint
     * @param j second edge endpoint
     */
    public void exclude(int i, int j) {
        nonEclAdjEdge[i]--;
        nonEclAdjEdge[j]--;
        if (nonEclAdjEdge[i] < 2 || nonEclAdjEdge[j] < 2) {
            this.inconsistency = true;
        }

        if (!clonedExclu[i]) {
            this.excluded[i] = this.excluded[i].clone();
            clonedExclu[i] = true;
        }
        if (!clonedExclu[j]) {
            this.excluded[j] = this.excluded[j].clone();
            clonedExclu[j] = true;
        }
        this.excluded[i][j] = true;
        this.excluded[j][i] = true;
    }

    /**
     * Add an edge to the mandatory edges of the node.
     *
     * @param i first end point of the edge
     * @param j second end point of the edge
     */
    public void forceEdge(int i, int j) {
        if (this.isMandatoryEdge(i, j)) return;

        if (!clonedMandatory[i]) {
            this.mandatoryEdges[i] = this.mandatoryEdges[i].clone();
            clonedMandatory[i] = true;
        }
        if (!clonedMandatory[j]) {
            this.mandatoryEdges[j] = this.mandatoryEdges[j].clone();
            clonedMandatory[j] = true;
        }

        if ((this.mandatoryEdges[i][1] != -1 && this.mandatoryEdges[i][0] != j && this.mandatoryEdges[i][1] != j) || (this.mandatoryEdges[j][1] != -1) && this.mandatoryEdges[j][0] != i && this.mandatoryEdges[j][1] != i) {
            throw new IllegalArgumentException();
        }

        if (this.mandatoryEdges[i][0] == -1) this.mandatoryEdges[i][0] = j;
        else this.mandatoryEdges[i][1] = j;

        if (this.mandatoryEdges[j][0] == -1) this.mandatoryEdges[j][0] = i;
        else this.mandatoryEdges[j][1] = i;
    }

    /**
     * Check if a edge is already mandatory.
     *
     * @param i first end point of the edge
     * @param j second end point of the edge
     * @return A boolean equals to true if the edge is already mandatory, false otherwise.
     */
    public boolean isMandatoryEdge(int i, int j) {
        if (this.mandatoryEdges[i][0] == j || this.mandatoryEdges[i][1] == j) return true;
        else return false;
    }
}

class NodeComparator implements Comparator<Node> {
    public int compare(Node a, Node b) {
        return Double.compare(a.lowerBound, b.lowerBound);
    }
}
