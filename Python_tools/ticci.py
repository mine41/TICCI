#!/usr/bin/env python
# coding: utf-8

import numpy as np
import networkx as nx
from matplotlib import cm
from scipy.sparse import csr_matrix
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl
from Python_tools.chuliu import chu_liu_edmond


class ticci:
    def __init__(self):
        pass

    def addCCI(adata, cci_df, cci_k):
        # Operate on connectivities
        arr = adata.obsp['connectivities'].toarray()  # Cell adjacency matrix
        myshape = arr.shape
        adatadf = adata.to_df()

        # cci weights
        print("Multiplier:", cci_k)
        for i in cci_df.iterrows():
            if i[0] % 10000 == 0:
                print(i[0], "/", len(cci_df))
            ligand_cell = i[1]['sender_cells']
            receptor_cell = i[1]['receiver_cells']
            prob = i[1]['prob']
            ligandLoc = adatadf.index.get_loc(ligand_cell)
            receptorLoc = adatadf.index.get_loc(receptor_cell)

            if arr[ligandLoc][receptorLoc] != 0:
                arr[ligandLoc][receptorLoc] += (prob * cci_k)

            if arr[receptorLoc][ligandLoc] != 0:
                arr[receptorLoc][ligandLoc] += (prob * cci_k)

        # New connectivities
        row, col = np.nonzero(arr)
        values = arr[row, col]
        adjacency = csr_matrix((values, (row, col)), shape=myshape)
        adata.obsp['connectivities'] = adjacency

    def getLineage(adata, stable_entropy_filepath, lineage_filepath):
        # Examine single-cell entropy
        sceList = []
        louvainList = []
        sce_threshold = 0.25

        for paga in adata.obs['louvain'].cat.categories:
            temp_scE = 0
            for i in range(0, int(len(
                    adata.obs['entropy'][adata.obs['louvain'] == paga]) * sce_threshold)):  # Take the top 25%
                temp_scE += adata.obs['entropy'][adata.obs['louvain'] == paga].sort_values()[i]
            louvainList.append(paga)
            sceList.append(temp_scE / int(len(adata.obs['entropy'][adata.obs['louvain'] == paga]) * sce_threshold))
        pagaSCE = pd.Series(sceList, index=louvainList)
        adata.uns['louvain_entropy'] = sceList

        # Normalize single-cell entropy
        maxE = pagaSCE.sort_values(ascending=False)[0]
        minE = pagaSCE.sort_values(ascending=True)[0]
        normalizationSCE = []
        for i in pagaSCE:
            normalizationSCE.append((i - minE) / (maxE - minE))

        adata.obs['stable_entropy'] = adata.obs['entropy']
        for i in pagaSCE.index:
            adata.obs['stable_entropy'][adata.obs['louvain'] == i] = pagaSCE[i]

        sc.pl.paga(adata, color=['stable_entropy'], save=stable_entropy_filepath)

        # Visualize normalized entropy
        nodes = []
        nodes_color_list = []
        for louvain in adata.obs['louvain'].cat.categories:
            node_color = cm.OrRd(normalizationSCE[int(louvain)])
            nodes.append(int(louvain))
            nodes_color_list.append(node_color)

        position = {}
        for i in range(0, len(adata.uns['paga']['pos'])):
            position[i] = adata.uns['paga']['pos'][i]

        # Process as directed
        arr = adata.uns['paga']['connectivities'].toarray()
        size = arr.shape[0]
        for row in range(0, size):
            for col in range(0, size):
                if pagaSCE[row] < pagaSCE[col] and arr[row][col] > 0:  # Delete rows where entropy is less than columns
                    arr[row][col] = 0

        # Generate dictionary graph structure G
        g = {}
        for row in range(0, size):
            nonzero = np.nonzero(arr[row])[0]
            if len(nonzero) > 0:
                weight_dict = {}
                for i in nonzero:
                    weight_dict[i] = arr[row][i]
                g[row] = weight_dict

        # Specify starting point
        g_root = int(pagaSCE.sort_values(ascending=False).index[0])

        # Run chu-liu
        h = chu_liu_edmond(g, g_root)

        # Convert dictionary tree result to array, then to sparse matrix, save to adata.uns['paga']['directed_tree']
        tree_arr = np.zeros((size, size))
        for i in h:
            for j in h[i]:
                tree_arr[i][j] = h[i][j]

        row, col = np.nonzero(tree_arr)
        values = tree_arr[row, col]
        csr_tree_arr = csr_matrix((values, (row, col)), shape=tree_arr.shape)
        adata.uns['paga']['directed_tree'] = csr_tree_arr

        # Process edges
        edges = []
        tree_arr = adata.uns['paga']['directed_tree'].toarray()
        row, col = np.nonzero(tree_arr)
        widths = tree_arr[row, col] * 5
        for i in range(0, len(row)):
            edge = (0, 0)
            if pagaSCE[row[i]] > pagaSCE[col[i]]:
                edge = (row[i], col[i])
            else:
                edge = (col[i], row[i])
            edges.append(edge)

        # Draw the graph
        G = nx.DiGraph()
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)

        # Draw network graph
        pos = nx.spring_layout(G)  # Support layout for displaying node positions

        node_colors = np.array(nodes_color_list)
        edge_widths = np.array(widths)
        vmin, vmax = pagaSCE.values.min(), pagaSCE.values.max()

        nx.draw(G,
                with_labels=True,
                pos=position,
                node_color=nodes_color_list,
                width=widths,
                cmap=cm.OrRd)

        # Add color bar
        sm = pl.cm.ScalarMappable(cmap=cm.get_cmap("OrRd"), norm=pl.Normalize(vmin=vmin, vmax=vmax))
        sm._A = []
        cbar = pl.colorbar(sm, label='scEntropy')

        # Save as SVG file
        pl.savefig(lineage_filepath, format='svg')
