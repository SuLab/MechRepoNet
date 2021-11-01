# edge\_semtypes.csv
Edges were mapped to a simple hierarchy of edge types.  All general, undirected associations were labeled "associated_with". Directed edges were divided into relationships that reflect "membership" (involved in, part of) and relationships that reflect "general regulation" (affects, regulates).  Regulatory relationships were further subdivided into those with a "positive regulatory effect" (activates) and those with a "negative regulatory effect" (inhibits).  All edges with either a positive or a negative regulatory effect were duplicated with a "general regulation" edge.

* `rel_dir` - indicates whether or not an edge has a direction. {-1, 0 ,1}. 1 indicates a positive directed edge, -1 indicates a negative directed edge, and 0 indicates either positive or negative directed edge (generalized)
* `directed` - indicates whether an edge has a general regulatory effect
* `parent_rel` - a mapping of edges from a directed case to a non-directed case.

These edges are used in identifying mechanistic paths from metapaths in MechRepoNet. This is done following the function highlighted https://github.com/mmayers12/data\_tools/blob/master/data\_tools/graphs/\_metapaths.py#L257 
