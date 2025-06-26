# -*- coding: utf-8 -*-
import os
import pymnet
import matplotlib.pyplot as plt

# own modules
from .Reader import Reader
from .inference import inferring_intralayer_edges, inferring_interlayer_edges


class MultiDGene:
    def __init__(
        self,
        genes_of_interest: list,
        tissues_of_interest: list = [],
        p_threshold: float = 0.05,
        data_path: str = None,
        mode: str = 'equal_exp',
        verbose: bool = False
    ):
        ###
        # mode:str can be defined as:
        # - 'equal_exp' (default) [edges refer to equal expressed genes between tissues (layers)]
        # - 'diff_exp' [edges refer to differential expressed genes between tissues (layers)]
        # - 'diff_co_exp' [edges refer to differential co-expressed pairs of genes between tissues (layers)]
        ###     

        if not genes_of_interest:
            raise ValueError("Parameter 'genes_of_interest' cannot be empty.")

        self.verbose = verbose
        self.data_path = data_path or os.path.dirname(os.path.realpath(__file__)) + "/"
        self.genes_of_interest = genes_of_interest
        self.tissues_of_interest = tissues_of_interest
        self.p_threshold = p_threshold
        self.mode = mode
        self.s_threshold = [0.5, 0.9]

        self._print_banner()
        self._initialize_data()

    def _print_banner(self):
        if not self.verbose:
            return

        banner = """
***************************************************************
 --------------------------------------------------------------
 >        multi-Dimensional Gene Network (multi-DGeNe)        <
 --------------------------------------------------------------
 Description:\tMulti-Dimensional Gene Network (multi-DGeNe)
              \tinference based on Weighted Gene Co-expression
              \tNetwork Analysis (WGCNA) across human tissues.
 --------------------------------------------------------------
 Version:\t0.1
 GitHub:\thttps://github.com/pietrocinaglia/multidgene
 Contact:\tPietro Cinaglia (cinaglia@unicz.it)
 --------------------------------------------------------------
***************************************************************
"""
        print(banner)
        print("[INFO] Initialization...")

    def _initialize_data(self):
        self.GTExData = Reader(data_path=self.data_path, verbose=self.verbose)

        self.gene_symbol_ids_of_interest, self.gene_code_ids_of_interest = self.GTExData.load_metadata_from_file(
            self.genes_of_interest
        )

        self.gene_expression_read_counts_by_tissue = self.GTExData.load_gene_expression_read_counts_from_file(
            gene_code_ids_of_interest=self.gene_code_ids_of_interest,
            tissues_of_interest=self.tissues_of_interest
        )

        self.gene_expression_tpm_by_tissue, self.tissues_of_interest = self.GTExData.retrieve_gene_expression_tpm_from_api(
            gene_code_ids_of_interest=self.gene_code_ids_of_interest,
            tissues_of_interest=self.tissues_of_interest
        )

        if self.verbose:
            print("[OK]\n")

    def construct_multidimensional_network(self) -> tuple[pymnet.MultilayerNetwork, set]:
        if self.verbose:
            print("[INFO] Constructing Multi-Dimensional Network...")

        net = pymnet.MultilayerNetwork(aspects=1, directed=False, fullyInterconnected=True)

        net, edge_colors = inferring_intralayer_edges(
            gene_symbol_ids_of_interest=self.gene_symbol_ids_of_interest,
            gene_exp_by_tissue=self.gene_expression_tpm_by_tissue,
            tissues_of_interest=self.tissues_of_interest
        )

        interlayer_input = self.gene_expression_read_counts_by_tissue if self.mode in ['diff_exp', 'equal_exp'] else self.gene_expression_tpm_by_tissue

        net, edge_colors = inferring_interlayer_edges(
            gene_symbol_ids_of_interest=self.gene_symbol_ids_of_interest,
            gene_exp_by_tissue=interlayer_input,
            tissues_of_interest=self.tissues_of_interest,
            mode=self.mode,
            p_threshold=self.p_threshold,
            multidimensional_network=net,
            edge_colors=edge_colors
        )

        if self.verbose:
            print("\n[INFO] Multi-Dimensional Network Successfully Modelled\n")

        return net, edge_colors

    def plot(self, multidimensional_network: pymnet.MultilayerNetwork, edge_colors: set = {}, output_path: str = None):
        fig = pymnet.draw(
            multidimensional_network,
            layershape='rectangle',
            layerPadding=0.2,
            defaultLayerColor="#BDBBBB",
            defaultLayerLabelLoc=(1.1, 1.1),
            defaultLayerLabelColor="#6E6E6E",
            defaultNodeLabelColor="#000000",
            defaultNodeColor="#000000",
            nodeSizeRule={'rule': 'scaled', 'scalecoeff': 0.2},
            defaultEdgeColor='black',
            edgeColorDict=edge_colors,
            defaultNodeLabelSize=6,
            alignedNodes=True,
            figsize=(7, 4),
            autoscale=True,
        )

        if output_path:
            filepath = os.path.join(output_path, "multidimensional_network.pdf")
            fig.savefig(filepath)
            print(f"- Plot saved as file [OK] -> {filepath}")

        plt.show()

    def save_as_files(self, multidimensional_network: pymnet.MultilayerNetwork, output_path: str):
        if not output_path:
            raise ValueError("You must provide an output path for storing data.")

        filepath = os.path.join(output_path, "multidimensional_network.txt")

        # Placeholder for future implementation
        print(f"[WARNING] Saving to file is not implemented yet. Expected path: {filepath}")
