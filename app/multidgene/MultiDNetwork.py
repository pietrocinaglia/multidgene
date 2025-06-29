# -*- coding: utf-8 -*-

import networkx as nx

import csv
from collections import defaultdict

import warnings

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference, bridging Co-Expression and Expression Dynamics across tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

# Own modules
from .VisualDesigner import VisualDesigner

class MultiDNetwork:
    layers = {}
    graph = None

    def __init__(self):
        self.layers = set()
        self.graph = nx.MultiDiGraph()

    def add_layer(self, layer_id: str):
        self.layers.add(layer_id)
    
    def add_layers(self, layer_ids: list):
        self.layers.update(layer_ids)

    def add_node(self, node_id: str, layer_id: str):
        if layer_id is None and self.layers is None:
            warnings.warn("You can not add nodes without define layers!")
            return
        
        self.add_layer(layer_id)
        self.graph.add_node((node_id, layer_id))

    def add_nodes(self, nodes: set, layer_id: str = None):
        if layer_id is None and self.layers is None:
            warnings.warn("You can not add nodes without define layers!")
            return

        if layer_id is None:
            for layer in self.layers:
                for node in nodes:
                    self.add_node(node, layer)
        else:
            for node in nodes:
                self.add_node(node, layer_id)

    def add_edge(self, node1_id: str, node2_id: str, layer1: str = None, layer2: str = None, weight: int = 1) -> bool:
        if layer2 is None:
            layer2 = layer1
        
        if (not self.has_node(node1_id,layer1)) or (not self.has_node(node2_id,layer2)):
            return False
        
        self.graph.add_edge((node1_id, layer1), (node2_id, layer2), weight=weight)
        return True

    def has_node(self, node_id:str, layer_id: str = None) -> bool:
        for (node, layer) in self.graph.nodes:
            if node_id == node and (layer_id is None or layer_id == layer):
                return True
        return False
    
    def get_nodes(self, layer_id: str = None) -> dict:
        result = defaultdict(set)
        for (node, layer) in self.graph.nodes:
            if layer_id is None or layer_id == layer:
                result[layer].add(node)
        return dict(result)

    def get_edges(self, layer_id: str = None) -> list:
        edges = []
        for (u, v, data) in self.graph.edges(data=True):
            (n1, l1) = u
            (n2, l2) = v
            if l1 == layer_id and l2 == layer_id:
                edges.append((n1, n2))
            elif layer_id is None:
                edges.append((n1, n2, l1, l2))
            else:
                continue
        return edges

    def get_layers(self) -> list:
        return list(self.layers)

    def write_edge_list(self, output_path, data: bool = False):
        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            for u, v, d in self.graph.edges(data=True):
                row = [u[0], v[0], u[1], v[1]]
                if data:
                    row.append(d.get('weight', 1))
                writer.writerow(row)

    def read_edge_list(self, input_path, data: bool = False):
        with open(input_path, 'r', newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                node1, node2, layer1, layer2 = row[:4]
                weight = int(row[4]) if data and len(row) > 4 else 1
                self.add_edge(node1, node2, layer1, layer2, weight)

    def write_adjacency_list(self, output_path):
        with open(output_path, 'w') as f:
            for node in self.graph.nodes:
                neighbors = list(self.graph.successors(node))
                line = f"{node[0]}:{node[1]}"
                for nbr in neighbors:
                    line += f" {nbr[0]}:{nbr[1]}"
                f.write(line + '\n')

    def read_adjacency_list(self, input_path):
        with open(input_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                node1_id, layer1 = parts[0].split(':')
                for neighbor in parts[1:]:
                    node2_id, layer2 = neighbor.split(':')
                    self.add_edge(node1_id, node2_id, layer1, layer2)

    def get_adjacency_matrix(self, layer_id: str):
        nodes = [n for n in self.graph.nodes if n[1] == layer_id]
        G_layer = self.graph.subgraph(nodes).to_directed()
        return nx.adjacency_matrix(G_layer, nodelist=nodes)
    
    def plot(self, output_path:str=None, format='eps', dpi='figure'):
        print("Plotting (2D)...")
        visual_designer = VisualDesigner(self, output_path, format, dpi='figure')
        visual_designer.plot()
    
    def plot3d(self, output_path:str=None, format='eps'):
        print("Plotting (3D)...")
        visual_designer = VisualDesigner(self, output_path, format)
        visual_designer.plot3d()

    def plot3d_as_html(self):
        print("Generating 3D Plot (HTML)...")
        visual_designer = VisualDesigner(self)
        return visual_designer.generate_3d_plotly()