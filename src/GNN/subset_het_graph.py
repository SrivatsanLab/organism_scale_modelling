import scipy
import torch
import torch_geometric as tg
from genomic_gnn.construct_het_graph import load_het_graph

def subset_het_graph(
    het_graph_path: str,
    genome_preserve_fraction: float, 
    gene_preserve_fraction: float = 1.0,   
    save_path: str | None = None,
    return_graph: bool = True,
):
    """
    Subset a heterogeneous graph to include only specified genomes and their associated genes.

    Parameters:
        het_graph_path: Path to the full heterogeneous graph file.
        genome_preserve_fraction: Fraction of genome nodes to retain in the subset graph.
        gene_preserve_fraction: Fraction of gene nodes to retain in the subset graph.
        save_path: Optional path to save the subset graph.
        return_graph: If True, return the in-memory subset graph object.
    Returns:
        torch_geometric.data.HeteroData | None: The subset heterogeneous graph if `return_graph` is True.
    """

    het_graph = load_het_graph(het_graph_path)

    num_genomes = het_graph['genome'].num_nodes
    num_genomes_to_preserve = int(num_genomes * genome_preserve_fraction)
    preserved_genome_indices = torch.randperm(num_genomes)[:num_genomes_to_preserve]

    num_genes = het_graph['gene'].num_nodes
    num_genes_to_preserve = int(num_genes * gene_preserve_fraction)
    preserved_gene_indices = torch.randperm(num_genes)[:num_genes_to_preserve]

    subset_dict = {
        'genome': preserved_genome_indices,
        'gene': preserved_gene_indices
    }
    if type(het_graph['genome'].x) == scipy.sparse._csr.csr_matrix:
        het_graph['genome'].x = het_graph['genome'].x.todense()
    subset_het_graph = het_graph.subgraph(subset_dict)

    # Remove gene nodes with no connections to preserved genomes
    genome_src, gene_dest = subset_het_graph['genome', 'present', 'gene'].edge_index
    preserved_gene_indices = torch.unique(gene_dest[torch.isin(genome_src, preserved_genome_indices)])
    subset_dict = {
        'gene': preserved_gene_indices
    }
    subset_het_graph = subset_het_graph.subgraph(subset_dict)

    subset_het_graph.validate(raise_on_error=True)
    if save_path is not None:
        save_het_graph(subset_het_graph, save_path)
    if return_graph:
        return subset_het_graph
