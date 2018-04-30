#include "ligra.h"
#include <stdio.h>
#include <stdlib.h>

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

long closest_clusters(int *similarity, long *C, long s, long n) {
    return (long)0;
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

        long *names = get_cluster_names(C, n);
        parallel_for(long s = 0; s < n; s++) {
            if (names[s] == 1) {
                long t = closest_clusters(similarity, C, s, n);
                parallel_for(int i = 0; i < n; i++) {
                    if (C[i] == s) {newC[i] = t;}
                }
            }
        }

        //compute "closest" cluster for every cluster (helper fn) 
        //invert fn p, find which clusters to be combined
        //create new C, making cluster names the lowest vertex in each cluster 
        

 	iter++;
        parallel_for(int i = 0; i < n; i++) {labels[i][iter] = newC[i];}
        converged = (num_clusters(C, n) == num_clusters(newC, n));
	parallel_for(int i = 0; i < n; i++) {C[i] = newC[i];}
    }

    //sort according to labels 
}

template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {

    long n = G.n;

    int *similarity = Similarity(G);
    long *permutation = randomPermutation(n); // replace with AffinityOrdering

    // semi local swaps 
}
