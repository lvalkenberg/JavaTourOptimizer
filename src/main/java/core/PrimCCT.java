package core;

import util.UF;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Stack;

public class PrimCCT {
    double[][] distanceMatrix;
    int[][] mandatoryEdges;
    boolean[][] edgesAdded;
    int size;
    double cost;
    int[] nodesDegree;
    CCTree ccTree;
    int[] parent;
    double[] minCost;
    boolean inconsistent = false;

    public PrimCCT(double[][] distanceMatrix) {
        this.distanceMatrix = distanceMatrix;
        this.size = distanceMatrix.length;
        this.nodesDegree = new int[size];
        this.cost = 0;
        this.parent = new int[size];
        computeMST();
    }

    public PrimCCT(double[][] distanceMatrix, int[][] mandatoryEdges) {
        this.distanceMatrix = distanceMatrix;
        this.mandatoryEdges = mandatoryEdges;
        this.size = distanceMatrix.length;
        this.edgesAdded = new boolean[size][2];
        this.nodesDegree = new int[size];
        this.cost = 0;
        this.parent = new int[size];
        computeMST();
    }

    /**
     * Compute the mst based on the distance matrix, including
     * the mandatory edges. This is a O(V^2) implementation of
     * the Prim minimum spanning tree algorithm.
     */
    public void computeMST() {
        // add fist node
        this.minCost = distanceMatrix[0].clone();
        for (int i = 0; i < size; i++) {
            parent[i] = 0;
        }
        addMandatoryEdgesFrom(0);
        if (inconsistent) return;

        for (int k = 1; k < size; k++) {
            int i;
            for (i = 1; i < size; i++) { if (nodesDegree[i] == 0) break; }
            if (i == size) return; // MST completed with the mandatory edges
            for (int j = i + 1; j < size; j++) {
                if (nodesDegree[j] == 0 && minCost[j] < minCost[i]) i = j;
            }
            addEdge(parent[i], i);
            for (int j = 1; j < size; j++) {
                if (nodesDegree[j] == 0 && distanceMatrix[i][j] < minCost[j]) {
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
        int o = -1;
        for (int k : mandatoryEdges[i + 1]) {
            o++;
            if (k == -1) return;
            if (k == 0 || edgesAdded[i][o]) continue;
            if (nodesDegree[k - 1] != 0 && !edgesAdded[i][o]) {
                inconsistent = true; // loop detected
                return;
            }
            k--; // shift from the 1tree

            addEdge(i, k);
            parent[k] = i;

            edgesAdded[i][o] = true;
            if (mandatoryEdges[k + 1][0] == i + 1) edgesAdded[k][0] = true;
            else edgesAdded[k][1] = true;

            // update the prim's cut weight
            for (int j = 1; j < size; j++) {
                if (nodesDegree[j] == 0 && distanceMatrix[k][j] < minCost[j]) {
                    minCost[j] = distanceMatrix[k][j];
                    parent[j] = k;
                }
            }

            addMandatoryEdgesFrom(k);
        }
    }

    /**
     * Add edge (i,j) to the mst
     *
     * @param i first end point of the edge
     * @param j second end point of the edge
     */
    private void addEdge(int i, int j) {
        cost += distanceMatrix[i][j];
        nodesDegree[i]++;
        nodesDegree[j]++;
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
            ccTree.updateCCTree(e.v, e.w, e);
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
     * @return  a double list of the replacement cost
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
        for(int i=0;i<size;i++){
            children.add(new ArrayList<>());
        }
        for(int i=1;i<size;i++){
            children.get(parent[i]).add(i);
        }
        // tree rooted in 0
        int[] depth = new int[size];
        Stack<Integer> rootedNodes = new Stack<>();
        rootedNodes.add(0);
        while(!rootedNodes.isEmpty()){
            int current = rootedNodes.pop();
            for(int child : children.get(current)){
                depth[child] = depth[current] + 1;
                rootedNodes.add(child);
            }
        }
        // compute the replacement edges
        double[] replacementCost = new double[size];
        boolean[] makred = new boolean[size];
        int[] pointer = new int[size];
        for(int i=0;i<size;i++){
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
