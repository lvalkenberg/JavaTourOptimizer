package core;

import util.UF;

import java.util.*;

public class PrimCCT {
    private int size;
    private double[][] distanceMatrix;
    private int[][] mandatoryEdges;
    private double weight;
    private boolean[][] edgesAdded;
    private int[] nodesDegree;
    private boolean[] mandatoryLock;
    private CCTree ccTree;
    public int[] parent;
    private double[] minCost;
    public boolean inconsistent = false;
    private UF uf;

    /**
     * Create an instance to solve MST on graph of size n
     *
     * @param size size of the graph
     */
    public PrimCCT(int size) {
        this.size = size;
        this.nodesDegree = new int[size];
        this.mandatoryLock = new boolean[size];
        this.weight = 0;
        this.parent = new int[size];
        this.edgesAdded = new boolean[size][2];
        uf = new UF(size);
    }

    /**
     * Constructor of the class
     *
     * @param distanceMatrix distance matrix used to compute the MST
     * @param mandatoryEdges mandatory edges to include in the MST
     */
    public PrimCCT(double[][] distanceMatrix, int[][] mandatoryEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        this.size = distanceMatrix.length;
        this.edgesAdded = new boolean[size][2];
        this.nodesDegree = new int[size];
        this.mandatoryLock = new boolean[size];
        this.weight = 0;
        this.parent = new int[size];
        uf = new UF(size);
        computeMST();
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. This is a O(V^2) implementation of
     * the Prim minimum spanning tree algorithm.
     *
     * @param distanceMatrix distance matrix used to compute the MST
     * @param mandatoryEdges mandatory edges to include in the MST
     */
    public void computeMST(double[][] distanceMatrix, int[][] mandatoryEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        for (int i = 0; i < size; i++) {
            nodesDegree[i] = 0;
            mandatoryLock[i] = false;
            parent[i] = 0;
            weight = 0;
            edgesAdded[i][0] = false;
            edgesAdded[i][1] = false;
        }
        this.inconsistent = false;
        this.weight = 0;
        computeMST();
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. This is a O(V^2) implementation of
     * the Prim minimum spanning tree algorithm.
     */
    public void computeMST() {
        // add fist node
        this.minCost = Arrays.copyOf(distanceMatrix[0], size);
        // 2 mandatory edges nodes
        for (int i = 0; i < size; i++) {
            if (mandatoryEdges[i + 1][1] != -1) {
                minCost[i] = Double.MAX_VALUE;
                mandatoryLock[i] = true;
            }
        }

        for (int i = 0; i < size; i++) {
            parent[i] = 0;
        }
        addMandatoryEdgesFrom(0);
        if (inconsistent) return;

        for (int k = 1; k < size; k++) {
            int i;
            for (i = 1; i < size; i++) {
                if (nodesDegree[i] == 0) break;
            }
            if (i == size) return; // MST completed with the mandatory edges
            for (int j = i + 1; j < size; j++) {
                if (nodesDegree[j] == 0 && minCost[j] < minCost[i]) i = j;
            }
            if (minCost[i] == Double.MAX_VALUE || nodesDegree[i] != 0) { // forbidden edge selected -> inconsistency
                inconsistent = true;
                return;
            }
            addEdge(parent[i], i);
            for (int j = 1; j < size; j++) {
                if (nodesDegree[j] == 0 && !mandatoryLock[j] && distanceMatrix[i][j] < minCost[j]) {
                    minCost[j] = distanceMatrix[i][j];
                    parent[j] = i;
                }
            }
            addMandatoryEdgesFrom(i);
            if (inconsistent) return;
        }
    }

    /**
     * Add mandatory (store in mandatoryEdges) edges adjacent with node i.
     *
     * @param i node to which adjacent mandatory edges are added
     */
    private void addMandatoryEdgesFrom(int i) {
        int offset = -1;
        int v = i + 1;
        for (int w : mandatoryEdges[i + 1]) {
            int j = w - 1;  // i,j are the node id in the mst and v,w are the shifted node from the 1-tree
            offset++;
            if (w == -1) return;    // no mandatory edge
            if (w == 0 || edgesAdded[i][offset]) continue;  // root edge or already added edge
            if (nodesDegree[j] != 0 && !edgesAdded[i][offset]) {    // subtour detected
                inconsistent = true;
                return;
            }
            //w--; // shift from the 1tree

            addEdge(i, j);
            parent[j] = i;

            edgesAdded[i][offset] = true;
            if (mandatoryEdges[w][0] == v) edgesAdded[j][0] = true;
            else edgesAdded[j][1] = true;

            // update the prim's cut weight
            for (int k = 1; k < size; k++) {
                if (nodesDegree[k] == 0 && !mandatoryLock[k] && distanceMatrix[j][k] < minCost[k]) {
                    minCost[k] = distanceMatrix[j][k];
                    parent[k] = j;
                }
            }

            addMandatoryEdgesFrom(j);
        }
    }

    /**
     * Add edge (i,j) to the mst
     *
     * @param i first end point of the edge
     * @param j second end point of the edge
     */
    private void addEdge(int i, int j) {
        weight += distanceMatrix[i][j];
        nodesDegree[i]++;
        nodesDegree[j]++;
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
        Edge compEdge = new Edge(-1, -1, -1);
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                if (distanceMatrix[i][j] == Double.MAX_VALUE) continue;
                int lca = ccTree.LCA(i, j);
                Edge sup = ccTree.tree[lca].edge;
                if (sup == null) throw new IllegalArgumentException();
                compEdge.v = i;
                compEdge.w = j;
                compEdge.weight = distanceMatrix[i][j];
                if (sup.equals(compEdge)) continue; // the support of a MST edge is not useful
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
        int[] pointer = new int[size];
        for (int i = 0; i < size; i++) {
            pointer[i] = i;
        }
        uf.reset();
        for (Edge e : nonTreeEdges) {
            if (uf.count() == 1) break;

            int i = e.v;
            int j = e.w;
            i = pointer[uf.find(i)];
            j = pointer[uf.find(j)];
            while (i != j) {
                if (depth[i] >= depth[j]) {
                    replacementCost[i] = e.weight - distanceMatrix[i][parent[i]];
                    pointer[uf.find(i)] = pointer[uf.find(parent[i])];
                    uf.union(i, parent[i]);
                    i = pointer[uf.find(i)];
                }
                if (depth[i] < depth[j]) {
                    replacementCost[j] = e.weight - distanceMatrix[j][parent[j]];
                    pointer[uf.find(j)] = pointer[uf.find(parent[j])];
                    uf.union(j, parent[j]);
                    j = pointer[uf.find(j)];
                }
            }
        }
        return replacementCost;
    }
}
