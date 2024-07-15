#!/usr/bin/env python
# coding: utf-8

# chu-liu code
def reverse_graph(G):
    '''Return the reversed graph where g[dst][src]=G[src][dst]'''
    g = {}
    for src in G.keys():
        for dst in G[src].keys():
            if dst not in g.keys():
                g[dst] = {}
            g[dst][src] = G[src][dst]
    return g


def build_max(rg, root):
    '''Find the max in-edge for every node except for the root.'''
    mg = {}
    for dst in rg.keys():
        if dst == root:
            continue
        max_ind = -100
        max_value = -100
        for src in rg[dst].keys():
            if rg[dst][src] >= max_value:
                max_ind = src
                max_value = rg[dst][src]
        mg[dst] = {max_ind: max_value}
    return mg


def find_circle(mg):
    '''Return the first circle if found, otherwise return None'''

    for start in mg.keys():
        visited = []
        stack = [start]
        while stack:
            n = stack.pop()
            if n in visited:
                C = []
                while n not in C:
                    C.append(n)
                    n = list(mg[n].keys())[0]
                return C
            visited.append(n)
            if n in mg.keys():
                stack.extend(list(mg[n].keys()))
    return None


def chu_liu_edmond(G, root):
    ''' G: dict of dict of weights
            G[i][j] = w means the edge from node i to node j has weight w.
            Assume the graph is connected and there is at least one spanning tree existing in G.
        root: the root node, has outgoing edges only.
    '''
    # reversed graph rg[dst][src] = G[src][dst]
    rg = reverse_graph(G)
    # root only has out edge
    rg[root] = {}
    # the maximum edge for each node other than root
    mg = build_max(rg, root)

    # check if mg is a tree (contains a circle)
    C = find_circle(mg)
    # if there is no circle, it means mg is what we want
    if not C:
        return reverse_graph(mg)
    # Now consider the nodes in the circle C as one new node vc
    all_nodes = G.keys()
    vc = max(all_nodes) + 1

    # The new graph G_prime with V_prime=V\C+{vc}
    V_prime = list(set(all_nodes) - set(C)) + [vc]
    G_prime = {}
    vc_in_idx = {}
    vc_out_idx = {}
    # Now add the edges to G_prime
    for u in all_nodes:
        for v in G[u].keys():
            # First case: if the source is not in the circle, and the dest is in the circle, i.e. in-edges for C
            # Then we only keep one edge from each node that is not in C to the new node vc with the largest difference (G[u][v] - list(mg[v].values())[0])
            # To specify, for each node u in V\C, there is an edge between u and vc if and only if there is an edge between u and any node v in C,
            # And the weight of edge u->vc = max_{v in C} (G[u][v] - mg[v].values) The second term represents the weight of max in-edge of v.
            # Then we record that the edge u->vc is originally the edge u->v with v=argmax_{v in C} (G[u][v] - mg[v].values)

            if (u not in C) and (v in C):
                if u not in G_prime.keys():
                    G_prime[u] = {}
                w = G[u][v] - list(mg[v].values())[0]
                if (vc not in G_prime[u]) or (vc in G_prime[u] and w > G_prime[u][vc]):
                    G_prime[u][vc] = w
                    vc_in_idx[u] = v
            # Second case: if the source is in the circle, but the dest is not in the circle, i.e out-edge for C
            # Then we only keep one edge from the new node vc to each node that is not in C
            # To specify, for each node v in V\C, there is an edge between vc and v iff there is an edge between any edge u in C and v.
            # And the weight of edge vc->v = max_{u in C} G[u][v]
            # Then we record that the edge vc->v originally the edge u->v with u=argmax_{u in C} G[u][v]
            elif (u in C) and (v not in C):
                if vc not in G_prime.keys():
                    G_prime[vc] = {}
                w = G[u][v]
                if (v not in G_prime[vc]) or (v in G_prime[vc] and w > G_prime[vc][v]):
                    G_prime[vc][v] = w
                    vc_out_idx[v] = u
            # Third case: if the source and dest are all not in the circle, then just add the edge to the new graph.
            elif (u not in C) and (v not in C):
                if u not in G_prime.keys():
                    G_prime[u] = {}
                G_prime[u][v] = G[u][v]
    # Recursively run the algorithm on the new graph G_prime
    # The result A should be a tree with nodes V\C+vc, then we just need to break the circle C and plug the subtree into A
    # To break the circle, we need to use the in-edge of vc, say u->vc to replace the original selected edge u->v,
    # where v was the original edge we recorded in the first case above.
    # Then if vc has out-edges, we also need to replace them with the original edges, recorded in the second case above.
    A = chu_liu_edmond(G_prime, root)
    print(A)
    all_nodes_A = list(A.keys())
    for src in all_nodes_A:
        # The number of out-edges varies, could be 0 or any number <=|V\C|
        if src == vc:
            for node_in in A[src].keys():
                orig_out = vc_out_idx[node_in]
                if orig_out not in A.keys():
                    A[orig_out] = {}
                A[orig_out][node_in] = G[orig_out][node_in]
        else:
            for dst in A[src]:
                # There must be only one in-edge to vc.
                if dst == vc:
                    orig_in = vc_in_idx[src]
                    A[src][orig_in] = G[src][orig_in]
                    del A[src][dst]
    del A[vc]

    # Now add the edges from the circle to the result.
    # Remember not to include the one with new in-edge
    for node in C:
        if node != orig_in:
            src = list(mg[node].keys())[0]
            if src not in A.keys():
                A[src] = {}
            A[src][node] = mg[node][src]
    return A
