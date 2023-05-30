package core;

import util.UF;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.PriorityQueue;

/**
 * Algorithhms 4th Edition by Robert Sedgewick, Kevin Wayne
 */
public class KruskalCCT {
    double[][] distanceMatrix;
    ArrayList<Edge> mst;
    double cost;
    int[] nodesDegree;
    CCTree ccTree;

    public KruskalCCT(double[][] distanceMatrix) {
        this.distanceMatrix = distanceMatrix;
        computeMST();
    }

    public void computeMST() {
        // new object each time !
        mst = new ArrayList<>();
        nodesDegree = new int[distanceMatrix.length];
        PriorityQueue<Edge> pq = new PriorityQueue<>();
        UF uf = new UF(distanceMatrix.length);
        ccTree = new CCTree(distanceMatrix.length);
        cost = 0;

        //add each edge in the pq
        for (int i = 0; i < distanceMatrix.length; i++) {
            for (int j = i + 1; j < distanceMatrix.length; j++) {
                pq.add(new Edge(i, j, distanceMatrix[i][j]));
            }
        }

        while (!pq.isEmpty() && mst.size() < distanceMatrix.length - 1) {
            Edge e = pq.poll();
            if (uf.connected(e.v, e.w)) continue;
            ccTree.updateCCTree(uf.find(e.v), uf.find(e.w), e);
            uf.union(e.v, e.w);
            addEdge(e);
        }

        ccTree.travelCCTree(); // inorder travelling of the ccTree
    }

    private void addEdge(Edge e) {
        nodesDegree[e.v] += 1;
        nodesDegree[e.w] += 1;
        cost += e.weight;
        mst.add(e);
    }

    public ArrayList<Edge> getMST() {
        return this.mst;
    }

    public double getCost() {
        return this.cost;
    }

    // filtering on the MST

    public HashMap<Edge, ArrayList<Edge>> computeAllSupport() {
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

}
