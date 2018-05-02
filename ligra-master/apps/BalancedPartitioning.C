#include "ligra.h"
#include <stdio.h>
#include <stdlib.h>

long num_vertices;
bool random_nonlocal;

intE *listToIndicator(uintE *vertices, uintE vertices_size, long n) {
    int *indicator = newA(int, n);
    parallel_for(int i = 0; i < vertices_size; i++) {
        indicator[vertices[i]] = 1;
    } 
    return indicator;
}

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

            parallel_for(long k = 0; k < n; k++) {
                 intersection[k] = v_ind[k] * u_ind[k];
            }
            similarity[n*i + v_neighbors[j]] = sequence::plusReduce(intersection, n); 
        }
    }
    return similarity;
}

long *randomPermutation(long n) {
    long *array = newA(long, n);
    parallel_for(long i = 0; i < n; i++) {
        array[i] = i;
    }

    for(long i = 0; i < n; i++) {
        int r = (rand() % (n - i)) + i;
        long temp = array[i];
        array[i] = array[r];
        array[r] = temp;
    }
    return array;
}

int num_clusters(long* C, long n) {
    int *found = newA(int, n);
    parallel_for(int i = 0; i < n; i++) {
        found[C[i]] = 1;
    }
    return sequence::plusReduce(found, n);
}


long *get_cluster_names(long *C, long n) {
    long *names = newA(long, n);
    parallel_for(int i = 0; i < n; i++) {
        names[C[i]] = 1;
    }
    return names;
}


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

long *inv_perm(long *perm, long n) {
  long *inv_perm = newA(long, n);
  parallel_for(long i = 0; i < n; i++) {
      uintE v = perm[i];
      inv_perm[v] = i;
  }
  return inv_perm;
}


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


template <class vertex>
void randomSwap(graph<vertex>& G, long *perm, int k) {
    long n = G.n;
    long part_size = n/k + !!(n%k);
    long *inverse = inv_perm(perm, n);
 
    long p = rand() % n;
    long q;
    if (rand() % 2) {
        q = (p + (rand() % part_size)) % n;
        if (q < p) { q = n - 1; }
    } else {
        q = (p - (rand() % part_size)) % n;
        if (q > p) { q = 0; } 
    }
    if (random_nonlocal) { q = rand() % n; }
    
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


template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {

    int k = 5;
    random_nonlocal = true;  // random swaps instead of local swaps

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
