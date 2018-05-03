#include "ligra.h"
#include <stdio.h>
#include <stdlib.h>

long num_vertices;

/*  Args: 
      vertices: a list of names of vertices
      vertices_size: length of the list vertices
      n: number of vertices in the graph
    Returns:
      List of length n, where result[i] = 1 if i is in vertices, 0 otherwise 
    Description:
      	Converts between frontier representations, 
        from vertices[i] = 1 indicating that vertex 1 is in the frontier to 
        vertices[i] = 1 indicating that vertex i is in the frontier 
*/
intE *listToIndicator(uintE *vertices, uintE vertices_size, long n) {
    int *indicator = newA(int, n);
    parallel_for(int i = 0; i < vertices_size; i++) {
        indicator[vertices[i]] = 1;
    } 
    return indicator;
}


/*  Args: 
      graph G: a graph G with n vertices
    Returns:
      List of length n*n, where result[i * n + j] represents the similarity of vertices i and j
    Description: 
      Computes the similarity for all vertex pairs (u, v) and stores them 
      in the returned similarity matrix, where similarity between u and v is 
      the number of neighbors they have in common
*/  
template <class vertex>
int *Similarity(graph<vertex>& G) {
    long n = G.n;
    int *similarity = newA(int, n * n);
    parallel_for(long i = 0; i < n; i++) {
        vertex v = G.V[i];
        uintE v_degree = v.getOutDegree();
	uintE *v_neighbors = (uintE*) v.getOutNeighbors();
        parallel_for(long j = 0; j < v_degree; j ++) {
	    vertex u = G.V[v_neighbors[j]];
            uintE u_degree = u.getOutDegree();
	    uintE *u_neighbors = (uintE*) u.getOutNeighbors();	    

            int *v_ind = listToIndicator(v_neighbors, v_degree, n);
	    int *u_ind = listToIndicator(u_neighbors, u_degree, n);            
            int *intersection = newA(int, n);

            //Multiplication here is being used like the "and" operator.
            parallel_for(long k = 0; k < n; k++) {
                 intersection[k] = v_ind[k] * u_ind[k];
            }
            similarity[n*i + v_neighbors[j]] = sequence::plusReduce(intersection, n); 
        }
    }
    return similarity;
}

/*  Args: 
      n: the length of the permutation to generate
    Returns:
      List of length n representing a random permutation of 0..n-1
    Description:
      Generates a random permutation of the numbers 0..n-1. 
*/
long *randomPermutation(long n) {
    long *array = newA(long, n);
    parallel_for(long i = 0; i < n; i++) {
        array[i] = i;
    }

    for(long i = 0; i < n; i++) {
        int r = (rand() % (n - i)) + i; //rand int in [i, n) 
        long temp = array[i];
        array[i] = array[r];
        array[r] = temp;
    }
    return array;
}

/*  Args:
      C: a list of length n where C[i] = s indicates that 
         vertex i is a member of the cluster that s represents 
      n: the number of vertices
    Returns:
      The number of unique clusters with members in C
    Description:
      Computes the number of clusters that C has members in.   
*/
int num_clusters(long* C, long n) {
    int *found = newA(int, n);
    parallel_for(int i = 0; i < n; i++) {
        found[C[i]] = 1;
    }
    return sequence::plusReduce(found, n);
}

/*  Args:
      C: a list of length n where C[i] = s indicates that 
         vertex i is a member of the cluster that s represents
      n: the number of vertices 
    Returns:
      List of length n where result[i] = 1 indicates that i is the
      representative element for some cluster with members in C.
*/
long *get_cluster_names(long *C, long n) {
    long *names = newA(long, n);
    parallel_for(int i = 0; i < n; i++) {
        names[C[i]] = 1;
    }
    return names;
}

/*  Args:
      similarity: an n*n similarity matrix where similarity[n * i + j] represents the 
        the similarity of vertices i and j
      C: a list of length n where C[i] = c indicates that 
         vertex i is a member of the cluster that c represents
      s: the representative of the cluster that we'd like to compute the closest_cluster for
      n: the number of vertices
    Returns:
      The representative element of the closest cluster to s.
    Description: 
      Computes the closest cluster to s by returning the cluster with the highest average 
      similarity to s.
*/
long closest_cluster(int *similarity, long *C, long s, long n) {
    long* names = get_cluster_names(C, n);
    double max_sim = 0;
    long arg_max_sim = s;

    long s_count = 0;
    long t_count;
    for (int i = 0; i < n; i++) {
        if (names[i] == s) { s_count += 1; }
    }

    for(int t = 0; t < n; t++) {
        t_count = 0;
        if (names[t] == 1) {
            long sim_sum = 0;
            double avg_sim = 0;
 
        for (int i = 0; i < n; i++) {
            if (names[i] == t) { t_count += 1; }
        }

            for (int i = 0; i < n; i++) {
                if (C[i] == s) {
                    for (int j = 0; j < n; j++) {
                        if (C[j] == t) {
                            sim_sum += similarity[i * n + j];
                        } 
                    }
                }
            }

            avg_sim = sim_sum / ((double)(s_count * t_count));

            if (avg_sim > max_sim) {
                max_sim = avg_sim;
                arg_max_sim = t;
            }
        }
    }
    return arg_max_sim;
}

/*  Args:
      a: label that tracks the representative elements of clusters that a vertex belonged to
         where a is a list of length n and a[i] = s indicates that this vertex was a part of 
         cluster s on the (n-i)th iteration of the algorithm
      b: see a 
    Returns:
      true if a < b, false otherwise.
    Description:
      Used as a custom comparator to help sort labels of vertices. 
*/
bool labelCompare(long *a, long *b) {
    long n = num_vertices;
    for (int i = n - 1; i >= 0; i--) {
        if (a[i] > b[i]) {
            return false;
        }
        else if (a[i] < b[i]) {
            return true;
        }
    }
}

/*  Args:
      similarity: a list of length n*n that 
      n: the number of vertices
    Returns:
      A list of length n that defines an ordering of vertices where result[i] holds the 
      name of the i'th vertex in the ordering.
    Description:
      Computes an ordering of the vertices using the affinity ordering algorithm.     
*/
long *affinityOrdering(int *similarity, long n) {
    long **labels = newA(long*, n);
    for(int i = 0; i < n; i++) { labels[i] = newA(long, n); }
    parallel_for(int i = 0; i < n; i++) { labels[i][0] = i; }    

    long *C = newA(long, n);
    parallel_for(int i = 0; i < n; i++) {C[i] = i;}
    
    long *newC = newA(long, n);

    int iter = 0;
    bool converged = false;

    while(!converged) {

        //compute closest cluster for every vertex, save into newC
        long *names = get_cluster_names(C, n);
        for(long s = 0; s < n; s++) {
            if (names[s] == 1) {
                long t = closest_cluster(similarity, C, s, n);
                parallel_for(int i = 0; i < n; i++) {
                    if (C[i] == s || C[i] == t || newC[i] == s || newC[i] == t) {
                        newC[i] = min(s,t); 
                    }
                }
            }
        }

 	iter++;
        parallel_for(int i = 0; i < n; i++) {labels[i][iter] = newC[i];}
        converged = (num_clusters(C, n) == num_clusters(newC, n));
	parallel_for(int i = 0; i < n; i++) {C[i] = newC[i];}
    }

    sort(labels, labels + n, labelCompare);
    
    long *result = newA(long, n);
    parallel_for(int i = 0; i < n; i++) { result[i] = labels[i][0];  }

    return result;
}

/*  Args: 
      perm: a list of length n representing a permutation of 0..n-1
      n: the length of the permutation
    Returns:
      A list of length n where result[i] = j indicates that perm[j] = 1
    Description:
      Returns the inversion of a permutation.
*/
long *inv_perm(long *perm, long n) {
  long *inv_perm = newA(long, n);
  parallel_for(long i = 0; i < n; i++) {
      uintE v = perm[i];
      inv_perm[v] = i;
  }
  return inv_perm;
}

/*  Args:
      G: a graph G with n vertices
      perm: a list of length n representing a permutation of the vertices in G
      k: the number of partitions to create
    Returns:
      A double representing the fraction of edges in G that are cut by the partition
      created by interpreting perm to have the first n/k vertices as the first 
      partition, the second n/k vertices as the second partition, and so on.
    Description:
      Given a graph and a permutation of its vertices, interprets the permutation
      as a k-way partition and returns the fraction of edges cut by this partition. 
*/
template <class vertex>
double countCutEdges(graph<vertex>& G, long *perm, int k) {
    long n = G.n;    
    long *inverse = inv_perm(perm, n);
    
    long num_cut = 0;
    long part_size = n/k + !!(n%k);
    
    for (long i = 0; i < n; i++) {
       	long v = perm[i];
	vertex vert = G.V[v];
        uintE degree = vert.getOutDegree();
	uintE* neighbors = (uintE*) vert.getOutNeighbors();
        for (int j = 0; j < degree; j++) {
            int p1 = i / part_size;
            int p2 = inverse[neighbors[j]] / part_size;
            if (p1 != p2) { num_cut++; }
        }
    }

    return num_cut/((double)(G.m * 2));
}

/*  Args:
      G: a graph with n vertices
      perm: a list representing an ordering of the vertices in G
      k: the number of partitions to create
   Returns:
      A list representing a permutation of the vertices in G
   Description:
      Given an ordering of the vertices of G, where the first n/k vertices
      are interpreted to be in the first partition and so on, for n*k rounds
      randomly swaps vertices i and j if the swap reduces the number of cut edges. 
*/
template <class vertex>
void randomSwap(graph<vertex>& G, long *perm, int k) {
    long n = G.n;
    long part_size = n/k + !!(n%k);
    long *inverse = inv_perm(perm, n);
 
    long p = rand() % n;
    long q = rand() % n;
    
    long i = perm[p];
    long j = perm[q];

    int i_part = p / part_size;
    int j_part = q / part_size;
    vertex i_vert = G.V[i];
    vertex j_vert = G.V[j];
    uintE i_degree = i_vert.getOutDegree();
    uintE j_degree = j_vert.getOutDegree(); 
    uintE *i_neighbors = (uintE*) i_vert.getOutNeighbors();
    uintE *j_neighbors = (uintE*) j_vert.getOutNeighbors();
    
    long orig_cut = 0;

    for (int x = 0; x < i_degree; x++) {
        uintE neigh = i_neighbors[x];
        int part = inverse[neigh] / part_size; 
        if (part != i_part) { orig_cut ++; }
    }
    for (int y = 0; y < j_degree; y++) {
        uintE neigh = j_neighbors[y];
        if (neigh != i) {
            int part = inverse[neigh] / part_size;
            if (part != j_part) { orig_cut ++; }
        }
    }
   
    long new_cut = 0;

    for (int x = 0; x < i_degree; x++) {
        uintE neigh = i_neighbors[x];
        int part = inverse[neigh] / part_size;
        if (part != j_part) { new_cut ++; }
    }
    
    for (int y = 0; y < j_degree; y++) {
        uintE neigh = j_neighbors[y];
        if (neigh != i) {
            int part = inverse[neigh] / part_size;
            if (part != i_part) { new_cut ++; }
        }
    }

    if (new_cut < orig_cut) {
        perm[p] = j;
        perm[q] = i;
    }
    
    return;
}

/*  Args: 
      G: a graph G with n vertices
      P: command line arguments
    Returns:
      void
    Description:
      Computes and compares the number of cut edges of a k-way partition generated randomly
      to one generated via the affinity ordering algorithm, as well as the improvement that
      making random swaps to these algorithms can provide.
*/ 
template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {

    int k = 5;

    long n = G.n;
    num_vertices = n;

    printf("Number of Vertices: %li, Number of Edges: %li\n", G.n, G.m);

    int *similarity = Similarity(G);

    long *random_perm = randomPermutation(n);
    long *affinity_perm = affinityOrdering(similarity, n);

    printf("Fraction of cut edges under random permutation: %f\n",
           countCutEdges(G, random_perm, k)); 

    printf("Fraction of cut edges under affinity permutation: %f\n",
           countCutEdges(G, affinity_perm, k));

    printf("\nAfter random swaps:\n\n");
    
    for(long i = 0; i < k*n; i++) { 
        randomSwap(G, random_perm, k); 
        randomSwap(G, affinity_perm, k);
    }
    
    printf("Fraction of cut edges under random permutation: %f\n",
           countCutEdges(G, random_perm, k));

    printf("Fraction of cut edges under affinity permutation: %f\n",
           countCutEdges(G, affinity_perm, k));

}
