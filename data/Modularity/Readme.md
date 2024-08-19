# Readme

Modularity clustering [1] is a special case of clique partitioning where the edge weights are computed from the degrees 
of an undirected graph. The dataset contains six instances of the modularity clustering problem that are also part
of the openGM benchmark dataset [2]. The sources of the underlying graphs are: 
- adjnoun [3],
- dolphins [4],
- football [5],
- karate [6],
- lesmis [7],
- polbooks [8].

These graphs can be obtained from [here](http://www-personal.umich.edu/~mejn/netdata/).

Note that typically, the edge costs are normalized by the number of edges in order to archive optimal solution values 
between 0 and 1. The edge costs in this dataset are not normalized in oder to preserve integer costs.
The cost for a pair of nodes u, v is computed as 
    
    2 * m * 1_{uv in E} - deg(u) * deg(v)

where m = |E| is the number of edges, 1_{uv in E} indicates whether uv is an edge in the graph and deg(v) is the degree 
of the node v in the graph.  


# References

[1] Brandes, U., Delling, D., Gaertler, M., Goerke, R., Hoefer, M., Nikoloski, Z., Wagner, D. On modularity clustering. 
IEEE Transactions on Knowledge and Data Engineering 20(2), 172-188 (2008).

[2] Jörg H. Kappes, Bjoern Andres, Fred A. Hamprecht, Christoph Schnörr, Sebastian Nowozin, Dhruv Batra, Sungwoong Kim, 
Thorben Kroeger, Bernhard X. Kausler, Jan Lellmann, Bogdan Savchynskyy, Nikos Komodakis, Carsten Rother. 
A Comparative Study of Modern Inference Techniques for Discrete Energy Minimization Problems.
International Journal of Computer Vision (2015).

[3] Newman, Mark EJ. Finding community structure in networks using the eigenvectors of matrices. 
Physical review E 74.3 (2006).

[4] D. Lusseau, K. Schneider, O. J. Boisseau, P. Haase, E. Slooten, and S. M. 
Dawson, Behavioral Ecology and Sociobiology 54, 396-405 (2003).

[5] Girvan, Michelle, and Mark EJ Newman. Community structure in social and biological networks. 
Proceedings of the national academy of sciences 99.12 (2002).

[6] W. W. Zachary. An information flow model for conflict and fission in small groups. Journal of Anthropological 
Research 33, 452-473 (1977).

[7] D. E. Knuth, The Stanford GraphBase: A Platform for Combinatorial Computing, Addison-Wesley, Reading, MA (1993).

[8] V. Krebs. http://www.orgnet.com/.
