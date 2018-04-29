#include "ligra.h"
#include <stdio.h>

intE *listToIndicator(intE *vertices, intE vertices_size, long n) {
    int *indicator = newA(int, n);
    parallel_for(long i = 0; i < vertices_size; i++) {
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
        intE v_degree = v.getOutDegree();
	intE *v_neighbors = v.getOutNeighbors();
        parallel_for(long j = 0; j < v_degree; j ++) {
	    vertex u = v_neighbors[j];
            intE u_degree = u.getOutDegree();
	    intE *u_neighbors = u.getOutNeighbors();	    

            intE *v_ind = listToIndicator(v_neighbors, v_degree, n);
	    intE *u_ind = listToIndicator(u_neighbors, u_degree, n);            
            intE *intersection = newA(int, n);

            parallel_for(long k = 0; k < n; k++) {
                 intersection[k] = v_ind[k] * u_ind[k];
            }

            similarity[n*v + u] = sequence::plusReduce(intersection, n); 

        }
    }
    return similarity;
}


template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {

    //int *similarity = Similarity(G);
    printf("Completed successfully.\n");
}
