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

template <class vertex>
long *affinityOrdering(graph<vertex>& G, int *similarity) {
    long *C = newA(long, G.n);
    long *labels = malloc(sizeof(labels));
}

template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {

    long n = G.n;

    int *similarity = Similarity(G);
    long *permutation = randomPermutation(n); // replace with AffinityOrdering

    // semi local swaps 
}
