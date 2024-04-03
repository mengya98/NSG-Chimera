import muon as mu
import numpy as np
import scanpy as sc
import scirpy as ir
from cycler import cycler
from matplotlib import cm as mpl_cm
from matplotlib import pyplot as plt


#Importing data
adata = sc.read_h5ad("/home/ug1379/TCR/Tcellsub2.h5ad")
adata_tcr = ir.io.read_10x_vdj("/home/ug1379/TCR/TCRcombind_filtered_contig_annotations_inmeta.csv")
mdata = mu.MuData({"gex": adata, "airr": adata_tcr})
#TCR Quality Control
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
mu.pp.filter_obs(
    mdata, "airr:chain_pairing", lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))
#Define clonotypes and clonotype clusters
ir.pp.ir_dist(mdata)
ir.tl.define_clonotypes(mdata, receptor_arms="all", dual_ir="primary_only",same_v_gene = False)
ir.tl.clonotype_network(mdata, min_cells=3)
mdata.obs.groupby("gex:tissue", dropna=False).size()
ax = ir.pl.clonotype_network(
    mdata, color="gex:Rename", base_size=5, label_fontsize=5, panel_size=(20, 20)
)
plt.show(ax)


ir.pp.ir_dist(
    mdata,
    metric="alignment",
    sequence="aa",
    cutoff=15)#sequences lower than cutoff will be connected


ir.tl.define_clonotype_clusters(
    mdata, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any")
ir.tl.clonotype_network(mdata, min_cells=3, min_nodes=6, sequence="aa", metric="alignment")
ax = ir.pl.clonotype_network(
    mdata, color="gex:sample", label_fontsize=9, panel_size=(20, 20), base_size=5)
plt.show(ax)

ax = ir.pl.clonotype_network(
    mdata, color="gex:Rename", label_fontsize=9, panel_size=(30, 30), base_size=6
)
plt.show(ax)

#Clonotype modularity
ir.tl.clonotype_modularity(mdata, target_col="airr:cc_aa_alignment")
clonotypes_top_modularity = list(
    mdata.obs.set_index("airr:cc_aa_alignment")["airr:clonotype_modularity"]
    .sort_values(ascending=False)
    .index.unique()
    .values[:2]
)
test_ad = mu.pl.embedding(
    mdata,
    basis="gex:umap",
    color="airr:cc_aa_alignment",
    groups=clonotypes_top_modularity,
    palette=cycler(color=mpl_cm.Dark2_r.colors),
)

