package core;

import util.UF;

/**
 * Connected Components Tree class
 * ( Simpler and incremental consistency checking
 * and arc consistency filtering algorithms for the
 * weighted spanning tree constraint, J-C Régin.)
 */
public class CCTree {
    int elementCounter;
    int size;
    int[] inorder;
    int[] pos;
    int[] depth;
    CCTNode[] tree; // contain the nodes of the ccTree
    int[] p; // pointer the root node of the connected component (index of vallue[])
    int globalCounter = 0;

    int[] log2Array;
    int[] pow2Array;
    int[][] M;

    UF uf;

    /**
     * ccTree constructor.
     *
     * @param numberOfLeaf  Number of node in the graph
     */
    public CCTree(int numberOfLeaf) {
        if (numberOfLeaf < 2) throw new IllegalArgumentException();

        this.elementCounter = numberOfLeaf - 1;
        this.size = numberOfLeaf * 2 - 1;
        this.inorder = new int[size];
        this.pos = new int[size];
        this.depth = new int[size];
        this.tree = new CCTNode[size];
        this.p = new int[numberOfLeaf];
        this.log2Array = new int[size + 1];
        this.pow2Array = new int[(int) (Math.log(size) / Math.log(2)) + 1];
        this.M = new int[size][];

        for (int i = 0; i < numberOfLeaf; i++) {
            tree[i] = new CCTNode(i, i);
            p[i] = i;
        }

        for (int i = 0; i < size; i++) {
            this.M[i] = new int[(int) (Math.log(size - i) / Math.log(2)) + 1];
        }

        for (int i = 1; i < size + 1; i++) {
            this.log2Array[i] = (int) (Math.log(i) / Math.log(2)); // ! log 0 = -inf
        }

        for (int i = 0; i <= Math.log(size) / Math.log(2); i++) {
            this.pow2Array[i] = (int) Math.pow(2, i);
        }

        this.uf = new UF(size);
    }

    /**
     * Update the CCTree when a new edge of the MST
     *
     * @param v      first node of the edge
     * @param w      second node of the edge
     * @param newEdge the new edge.
     */
    public void updateCCTree(int v, int w, Edge newEdge) { // TOD : simplify
        int rv = uf.find(v);
        int rw = uf.find(w);
        uf.union(rv, rw);

        elementCounter++;

        CCTNode left = tree[p[rv]];
        CCTNode right = tree[p[rw]];
        CCTNode newNode = new CCTNode(left, right, newEdge, elementCounter);

        p[rv] = elementCounter;
        p[rw] = elementCounter; // rv or rw will be the new canonical element
        tree[elementCounter] = newNode;
    }

    /**
     * Travel the CCTree to setup inorder, pos and height vectors.
     */
    public void travelCCTree() {
        globalCounter = 0;
        CCTNode root = tree[elementCounter]; // last added value is the root of the tree
        inorderTravelling(root, size);
    }

    /**
     * Inorder visit of the ccTree.
     *
     * @param node node to visit.
     * @param h    height of the node.
     */
    private void inorderTravelling(CCTNode node, int h) {
        if (node.left != null) inorderTravelling(node.left, h + 1);
        inorder[globalCounter] = node.index;
        depth[globalCounter] = h;
        pos[node.index] = globalCounter;
        globalCounter++;
        if (node.right != null) inorderTravelling(node.right, h + 1);

    }

    /**
     * Precompute the RMQ.
     */
    public void precomputeRMQ() {
        for (int i = 0; i < size; i++) {
            M[i][0] = i;
        }

        for (int j = 1; j <= log2Array[size]; j++) {
            for (int i = 0; i <= size - pow2Array[j]; i++) {
                int minL = M[i][j - 1];
                int minR = M[i + pow2Array[j - 1]][j - 1];
                M[i][j] = depth[minL] <= depth[minR] ? minL : minR;
            }
        }
    }

    /**
     * Compute the range min query with the precomputed RMQ.
     *
     * @param i lower bound
     * @param j upper bound
     * @return the index of the RMQ inorder.
     */
    public int RMQ(int i, int j) {
        int logWidth = log2Array[j - i + 1];
        int minL = M[i][logWidth];
        int minR = M[j - pow2Array[logWidth] + 1][logWidth];
        return depth[minL] <= depth[minR] ? minL : minR;
    }

    /**
     * Return the last common ancestor between i and j in the ccTree.
     *
     * @param i lower bound of the interval.
     * @param j upper bound of the interval.
     * @return index of the lca node in the ccTree.
     */
    public int LCA(int i, int j) {
        int posi = pos[i];
        int posj = pos[j];
        return posi <= posj ? inorder[RMQ(posi, posj)] : inorder[RMQ(posj, posi)];
    }

}

class CCTNode {
    public boolean isLeaf;
    public CCTNode left;
    public CCTNode right;
    public CCTNode parent;
    public Edge edge; // if not a leaf
    public int vertice; // if leaf
    public int index; // position of the node in de tree array

    /**
     * Leaf constructor
     *
     * @param vetice    node id
     * @param index     index in tree[]
     */
    public CCTNode(int vetice, int index) {
        this.isLeaf = true;
        this.vertice = vetice;
        this.index = index;
    }

    /**
     * Intermediate node constructor
     */
    public CCTNode(CCTNode left, CCTNode right, Edge edge, int index) {
        this.left = left;
        this.right = right;
        this.edge = edge;
        this.index = index;

        right.parent = this;
        left.parent = this;
    }
}

class Edge implements Comparable<Edge> {
    int v;
    int w;
    double weight;

    public Edge(int v, int w, double weight) {
        this.v = v;
        this.w = w;
        this.weight = weight;
    }

    @Override
    public int compareTo(Edge o) {
        if (this.weight < o.weight) return -1;
        else if (this.weight > o.weight) return +1;
        else return 0;
    }

    @Override
    public boolean equals(Object o) {

        // If the object is compared with itself then return true
        if (o == this) {
            return true;
        }

        /* Check if o is an instance of Complex or not
          "null instanceof [type]" also returns false */
        if (!(o instanceof Edge)) {
            return false;
        }

        // typecast o to Complex so that we can compare data members
        Edge e = (Edge) o;

        // Compare the data members and return accordingly
        return (Integer.compare(v, e.v) == 0 && Integer.compare(w, e.w) == 0) || (Integer.compare(v, e.w) == 0 && Integer.compare(w, e.v) == 0);
    }
}