package core;

import util.UF;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Stack;

public class MaxDegTwoMST {
    double[][] distanceMatrix;
    int size;
    double cost;
    int[] nodesDegree;
    int[][] mandatoryEdges;
    ArrayList<Edge> mst;
    CCTree ccTree;
    int[] parent;
    boolean inconsistent = false;
    Edge[] edgeArray;
    UF uf;
    HashMap<Integer, ArrayList<Integer>> mstAdjList = new HashMap<>();
    Stack<Integer> stack = new Stack();
    boolean[] visited;

    /**
     * Create an instance to solve MST on graph of size n
     *
     * @param size size of the graph
     */
    public MaxDegTwoMST(int size) {
        this.size = size;
        this.nodesDegree = new int[size];
        this.cost = 0;
        this.parent = new int[size];
        this.visited = new boolean[size];
        this.edgeArray = new Edge[size * size];
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                edgeArray[i + j * size] = new Edge(i, j, 0);
            }
        }
        this.uf = new UF(size);
        this.mst = new ArrayList<>();
    }

    /**
     * Constructor of the class
     *
     * @param distanceMatrix distance matrix used to compute the MST
     * @param mandatoryEdges mandatory edges to include in the MST
     */
    public MaxDegTwoMST(double[][] distanceMatrix, int[][] mandatoryEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        this.size = distanceMatrix.length;
        this.nodesDegree = new int[size];
        this.cost = 0;
        this.parent = new int[size];
        this.visited = new boolean[size];
        this.edgeArray = new Edge[size * size];
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                edgeArray[i + j * size] = new Edge(i, j, 0);
            }
        }
        this.uf = new UF(size);
        this.mst = new ArrayList<>();
        computeMST(5);
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. Based on the Boruvka MST algorithm
     *
     * @param distanceMatrix distance matrix used to compute the MST
     * @param mandatoryEdges mandatory edges to include in the MST
     * @param rootEdges      Edges on the root of the 1-tree took in account to
     *                       maximise the number of degree-2 nodes
     */
    public MaxDegTwoMST(double[][] distanceMatrix, int[][] mandatoryEdges, int[] rootEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        this.size = distanceMatrix.length;
        this.nodesDegree = new int[size];
        this.cost = 0;
        this.parent = new int[size];
        this.visited = new boolean[size];
        this.edgeArray = new Edge[size * size];
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                edgeArray[i + j * size] = new Edge(i, j, 0);
            }
        }
        this.uf = new UF(size);
        this.mst = new ArrayList<>();

        nodesDegree[rootEdges[0]] += 1;
        nodesDegree[rootEdges[1]] += 1;

        computeMST(5);
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. Based on the Boruvka MST algorithm
     *
     * @param distanceMatrix distance matrix used to compute the MST
     * @param mandatoryEdges mandatory edges to include in the MST
     * @param rootEdges      Edges on the root of the 1-tree took in account to
     *                       maximise the number of degree-2 nodes
     */
    public void computeMST(double[][] distanceMatrix, int[][] mandatoryEdges, int[] rootEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        for (int i = 0; i < size; i++) {
            nodesDegree[i] = 0;
            parent[i] = 0;
            visited[i] = false;
        }
        this.inconsistent = false;
        this.cost = 0;
        this.uf.reset();

        nodesDegree[rootEdges[0] - 1] += 1;
        nodesDegree[rootEdges[1] - 1] += 1;

        computeMST(5);
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. Based on the Boruvka MST algorithm
     *
     * @param alt alternative of how the equally ranked safe
     *            edges should be sorted. Default use 5
     */
    private void computeMST(int alt) {
        mst.clear();
        uf.reset();

        // update edges weight
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                edgeArray[i + j * size].weight = distanceMatrix[i][j];
            }
        }

        addMandatoryEdges();

        for (int t = 1; t < size && mst.size() < size - 1; t *= 2) { // O(log(V))
            //Edge[] closest = new Edge[graphSize];
            ArrayList<Edge>[] closest = new ArrayList[size];
            for (int i = 0; i < size; i++) {
                closest[i] = new ArrayList<>();
            }
            Edge[] lowerWeight = new Edge[size];
            for (int v = 0; v < size; v++) { // O(V^2)
                for (int w = v + 1; w < size; w++) {
                    int i = uf.find(v), j = uf.find(w);
                    if (i == j) continue;
                    Edge e = edgeArray[v + w * size];
                    if (lowerWeight[i] == null || e.compareTo(lowerWeight[i]) < 0) lowerWeight[i] = e;
                    if (lowerWeight[j] == null || e.compareTo(lowerWeight[j]) < 0) lowerWeight[j] = e;
                }
            }

            for (int v = 0; v < size; v++) {// O(V^2)
                for (int w = v + 1; w < size; w++) {
                    int i = uf.find(v), j = uf.find(w); //posible O(1)
                    if (i == j) continue;
                    Edge e = edgeArray[v + w * size];
                    if (lowerWeight[i].weight == e.weight) closest[i].add(e);
                    if (lowerWeight[j].weight == e.weight) closest[j].add(e);
                }
            }

            ArrayList<ArrayList<Edge>> CC_edges = new ArrayList<>();
            for (int i = 0; i < size; i++) {
                if (closest[i].size() != 0) CC_edges.add(closest[i]);
            }

            if (alt == 1) CC_edges.sort((a, b) -> Integer.compare(a.size(), b.size()));
            if (alt == 2) CC_edges.sort((a, b) -> Integer.compare(-a.size(), -b.size()));
            if (alt == 3) Collections.shuffle(CC_edges);
            if (alt == 4) {
                CC_edges = new ArrayList<>();
                for (int i = size - 1; i >= 0; i--) {
                    if (closest[i].size() != 0) CC_edges.add(closest[i]);
                }
            }
            if (alt == 5) CC_edges.sort((a, b) -> Integer.compare(Math.min(a.size(), 2), Math.min(b.size(), 2)));

            // add newly discovered edges to MST
            //for (int i = 0; i < graphSize; i++) { //O(V)
            for (ArrayList<Edge> candidates : CC_edges) { //O(V)
                //Edge e = closest[i];
                //ArrayList<Edge> candidates = closest[i];
                if (candidates.size() != 0) {
                    //Collections.sort(candidates); // O(V'*log(V')) -> O(V') en ne gardant que la maximum
                    Edge e = candidates.get(0);
                    for (Edge o : candidates) {
                        if (uf.find(o.v) != uf.find(o.w) && nodesDegree[e.v] + nodesDegree[e.w] < nodesDegree[o.v] + nodesDegree[o.w]) {
                            e = o;
                        }
                    }

                    int v = e.v, w = e.w;
                    // don't add the same edge twice
                    if (uf.find(v) != uf.find(w)) {
                        addEdge(e);
                        //uf.union(v, w);
                        uf.union_complete(v, w);
                    }
                }
            }
        }
        // create the parent[] array
        mstAdjList.clear();
        for (int i = 0; i < size; i++) visited[i] = false;
        for (Edge e : mst) {
            if (!mstAdjList.containsKey(e.v)) mstAdjList.put(e.v, new ArrayList<>());
            if (!mstAdjList.containsKey(e.w)) mstAdjList.put(e.w, new ArrayList<>());
            mstAdjList.get(e.v).add(e.w);
            mstAdjList.get(e.w).add(e.v);
        }
        stack.clear();
        visited[0] = true;
        stack.add(0);
        while (!stack.empty()) {
            int current = stack.pop();
            for (int next : mstAdjList.get(current)) {
                if (visited[next]) parent[current] = next;
                else stack.add(next);
            }
            visited[current] = true;
        }
    }

    /**
     * Add edge (i,j) to the mst
     *
     * @param i first end point of the edge
     * @param j second end point of the edge
     */
    private void addEdge(Edge e) {
        nodesDegree[e.v] += 1;
        nodesDegree[e.w] += 1;
        mst.add(e);
        cost += e.weight;
    }

    /**
     * Add mandatory edges of the MST.
     */
    private void addMandatoryEdges() {
        for (int i = 1; i < mandatoryEdges.length; i++) {
            for (int j = 0; j < 2; j++) {
                if (mandatoryEdges[i][j] != -1) {
                    int k = i - 1;
                    int l = mandatoryEdges[i][j] - 1;
                    if (k > l) continue;
                    if (uf.connected(k, l)) {
                        inconsistent = true;
                        return;
                    }
                    uf.union(k, l);
                    addEdge(edgeArray[k + l * size]);
                }
            }
        }
    }

    /**
     * Compute the support edge for each non minimum 1-tree edge of the current node.
     *
     * @return a hashmap : support edge -> non MST edges
     */
    public HashMap<Edge, ArrayList<Edge>> computeAllSupportEdges() {
        // construct the CCT (add  edges in increasing weight order)
        ccTree = new CCTree(distanceMatrix.length);
        ArrayList<Edge> mst = new ArrayList<>();
        for (int i = 1; i < size; i++) {
            mst.add(new Edge(i, parent[i], distanceMatrix[i][parent[i]]));
        }
        Collections.sort(mst);
        for (Edge e : mst) {
            ccTree.updateCCTree(e);
        }
        ccTree.travelCCTree();
        ccTree.precomputeRMQ();
        HashMap<Edge, ArrayList<Edge>> supports = new HashMap<>(); // ! only contain edges (i,j) st i<j
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                if (distanceMatrix[i][j] == Double.MAX_VALUE) continue;
                int lca = ccTree.LCA(i, j);
                Edge sup = ccTree.tree[lca].edge;
                if (sup.equals(new Edge(i, j, distanceMatrix[i][j])))
                    continue; // the support of a MST edge is not useful
                if (sup == null) throw new IllegalArgumentException();
                if (!supports.containsKey(sup)) supports.put(sup, new ArrayList<>());
                supports.get(sup).add(new Edge(i, j, distanceMatrix[i][j]));
            }
        }
        return supports;
    }

    /**
     * Compute the replacement cost of all minimum 1-tree edges of the current node.
     *
     * @return a double list of the replacement cost
     */
    public double[] computeReplacementCost() {
        // get all non tree edges ordered depending of the pi O(n log(n))
        ArrayList<Edge> nonTreeEdges = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            for (int j = i + 1; j < size; j++) {
                if (distanceMatrix[i][j] < Double.MAX_VALUE && parent[i] != j && parent[j] != i)
                    nonTreeEdges.add(new Edge(i, j, distanceMatrix[i][j]));
            }
        }
        Collections.sort(nonTreeEdges);

        ArrayList<ArrayList<Integer>> children = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            children.add(new ArrayList<>());
        }
        for (int i = 1; i < size; i++) {
            children.get(parent[i]).add(i);
        }
        // tree rooted in 0
        int[] depth = new int[size];
        Stack<Integer> rootedNodes = new Stack<>();
        rootedNodes.add(0);
        while (!rootedNodes.isEmpty()) {
            int current = rootedNodes.pop();
            for (int child : children.get(current)) {
                depth[child] = depth[current] + 1;
                rootedNodes.add(child);
            }
        }
        // compute the replacement edges
        double[] replacementCost = new double[size];
        boolean[] makred = new boolean[size];
        int[] pointer = new int[size];
        for (int i = 0; i < size; i++) {
            pointer[i] = i;
        }
        UF uft = new UF(size);
        for (Edge e : nonTreeEdges) {
            if (uft.count() == 1) break;

            int i = e.v;
            int j = e.w;
            i = pointer[uft.find(i)];
            j = pointer[uft.find(j)];
            while (i != j) {
                if (depth[i] >= depth[j]) {
                    replacementCost[i] = e.weight - distanceMatrix[i][parent[i]];
                    pointer[uft.find(i)] = pointer[uft.find(parent[i])];
                    uft.union(i, parent[i]);
                    i = pointer[uft.find(i)];
                }
                if (depth[i] < depth[j]) {
                    replacementCost[j] = e.weight - distanceMatrix[j][parent[j]];
                    pointer[uft.find(j)] = pointer[uft.find(parent[j])];
                    uft.union(j, parent[j]);
                    j = pointer[uft.find(j)];
                }
            }
        }
        return replacementCost;
    }

}

