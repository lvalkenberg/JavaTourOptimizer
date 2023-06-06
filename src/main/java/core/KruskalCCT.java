package core;

import util.UF;

import java.util.*;

/**
 * Algorithhms 4th Edition by Robert Sedgewick, Kevin Wayne
 */
public class KruskalCCT {
    double[][] distanceMatrix;
    int[][] mandatoryEdges;
    ArrayList<Edge> mst;
    double cost;
    int[] nodesDegree;
    CCTree ccTree;
    int size;
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
    public KruskalCCT(int size) {
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
    public KruskalCCT(double[][] distanceMatrix, int[][] mandatoryEdges) {
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
        computeMST();
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. Based on the Kruskal MST algorithm
     *
     * @param distanceMatrix    matrix distance
     * @param mandatoryEdges    mandatory edges include in every MST
     */
    public void computeMST(double[][] distanceMatrix, int[][] mandatoryEdges) {
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
        computeMST();
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. Based on the Kruskal MST algorithm
     */
    public void computeMST() {
        // new object each time !
        mst.clear();
        PriorityQueue<Edge> pq = new PriorityQueue<>();
        //add each edge in the pq
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                edgeArray[i + j * size].weight = distanceMatrix[i][j];
                pq.add(edgeArray[i + j * size]);
            }
        }

        addMandatoryEdges();
        if (inconsistent) return;

        while (!pq.isEmpty() && mst.size() < distanceMatrix.length - 1) {
            Edge e = pq.poll();
            if (uf.connected(e.v, e.w)) continue;
            //ccTree.updateCCTree(uf.find(e.v), uf.find(e.w), e);
            uf.union(e.v, e.w);
            addEdge(e);
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
     * Add edge e to the mst
     *
     * @param e new edge of the MST
     */
    private void addEdge(Edge e) {
        nodesDegree[e.v] += 1;
        nodesDegree[e.w] += 1;
        cost += e.weight;
        mst.add(e);
    }

    /**
     * Add mandatory (store in mandatoryEdges) edges adjacent with node i.
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
     * @return a hashmap = support edge : non MST edges
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
