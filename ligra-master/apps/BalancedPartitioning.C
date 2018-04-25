#include "ligra.h"
#include <stdio.h>

template <class vertex>
void Compute(graph<vertex>& G, commandLine P) {
	long n = G.n;
	long m = G.m;
	printf("Number of vertices in G: %li\n", n);
	printf("Number of edges in G: %li\n", m);	
}
