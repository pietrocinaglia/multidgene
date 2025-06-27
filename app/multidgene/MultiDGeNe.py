# -*- coding: utf-8 -*-

import os

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference, bridging Co-Expression and Expression Dynamics across tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

# Own modules
from .MultiDNetwork import MultiDNetwork
from .Reader import Reader
from .inference import inferring_intralayer_edges, inferring_interlayer_edges

class MultiDGeNe:
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
        # Default
        self.s_threshold = [0.5, 0.9]

        if self.verbose:
            self._show_banner()

        self._initialization()

    def _show_banner(self):
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

    def _initialization(self):
        if self.verbose:
            print("[INFO] Initialization...")

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

    def construct_multidimensional_network(self) -> tuple[MultiDNetwork, set]:
        if self.verbose:
            print("[INFO] Constructing Multi-Dimensional Network...")

        multidimensional_network = MultiDNetwork()
        multidimensional_network.add_layers(self.tissues_of_interest) #add to all layers
        multidimensional_network.add_nodes(self.gene_symbol_ids_of_interest) #add to all layers

        multidimensional_network, edge_colors = inferring_intralayer_edges(
            gene_symbol_ids_of_interest=self.gene_symbol_ids_of_interest,
            gene_exp_by_tissue=self.gene_expression_tpm_by_tissue,
            tissues_of_interest=self.tissues_of_interest,
            multidimensional_network=multidimensional_network,
            verbose=self.verbose
        )

        interlayer_input = self.gene_expression_read_counts_by_tissue if self.mode in ['diff_exp', 'equal_exp'] else self.gene_expression_tpm_by_tissue

        multidimensional_network, edge_colors = inferring_interlayer_edges(
            gene_symbol_ids_of_interest=self.gene_symbol_ids_of_interest,
            gene_exp_by_tissue=interlayer_input,
            tissues_of_interest=self.tissues_of_interest,
            mode=self.mode,
            p_threshold=self.p_threshold,
            multidimensional_network=multidimensional_network,
            edge_colors=edge_colors,
            verbose=self.verbose
        )

        if self.verbose:
            print("\n[INFO] Multi-Dimensional Network Successfully Modelled\n")

        return multidimensional_network, edge_colors
    