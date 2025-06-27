# -*- coding: utf-8 -*-

import networkx as nx
import csv
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.lines as mlines
import warnings

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference, bridging Co-Expression and Expression Dynamics across tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

class MultiDNetwork:
    layers = {}

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
        pos = {}
        layer_positions = {}
        color_map = list(mcolors.TABLEAU_COLORS.values())
        layer_list = sorted(self.layers)
        spacing = 5

        fig, ax = plt.subplots()

        for i, layer in enumerate(layer_list):
            y_offset = i * spacing
            layer_nodes = [n for n in self.graph.nodes if n[1] == layer]
            layer_pos = nx.spring_layout(self.graph.subgraph(layer_nodes), seed=42)
            for node, coord in layer_pos.items():
                pos[node] = (coord[0], coord[1] + y_offset)

            layer_positions[layer] = color_map[i % len(color_map)]
            
            # Draw layer boundary
            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            if x_vals and y_vals:
                min_x, max_x = min(x_vals), max(x_vals)
                min_y, max_y = min(y_vals), max(y_vals)
                width = (max_x - min_x + 1) + 0.0
                height = (max_y - min_y + 1) + 1.0
                ellipse = patches.Ellipse(((min_x + max_x) / 2, (min_y + max_y) / 2), width, height,
                                        color=layer_positions[layer], alpha=0.2)
                ax.add_patch(ellipse)
                ax.text((min_x + max_x) / 2, max_y + 1.2, layer, ha='center', fontsize=8, weight='bold')

        # Draw nodes (all in black)
        nx.draw_networkx_nodes(self.graph, pos, node_color='grey', node_size=250, ax=ax)

        # Draw intralayer edges
        for layer in layer_list:
            layer_edges = [(u, v) for u, v, _ in self.graph.edges(data=True) if u[1] == layer and v[1] == layer]
            nx.draw_networkx_edges(self.graph, pos, edgelist=layer_edges,
                                edge_color=layer_positions[layer], ax=ax, arrows=False)

        # Interlayer edges with unique colors
        layer_pairs = set((u[1], v[1]) for u, v, _ in self.graph.edges(data=True) if u[1] != v[1])
        
        #interlayer_colors = list(mcolors.CSS4_COLORS.values())
        interlayer_colors = [c for c in mcolors.CSS4_COLORS.values() if mcolors.rgb_to_hsv(mcolors.to_rgb(c))[2] < 0.9]
        interlayer_colors = {pair: interlayer_colors[i % len(interlayer_colors)] for i, pair in enumerate(layer_pairs)}

        for (u, v, _) in self.graph.edges(data=True):
            if u[1] != v[1]:
                pair = (u[1], v[1])
                color = interlayer_colors.get(pair, 'grey')
                nx.draw_networkx_edges(self.graph, pos, edgelist=[(u, v)], edge_color=color,
                                    style='dashed', ax=ax, arrows=False)

        # Draw node labels
        nx.draw_networkx_labels(self.graph, pos, font_size=4, font_color="black",
                                labels={n: n[0] for n in self.graph.nodes}, ax=ax)


        # Create legend handles
        legend_handles = []

        # Legend for intralayer edges
        for layer in layer_list:
            color = layer_positions[layer]
            handle = mlines.Line2D([], [], color=color, label=f"{layer}", linewidth=2)
            legend_handles.append(handle)

        # Legend for interlayer edges
        for (layer1, layer2), color in interlayer_colors.items():
            handle = mlines.Line2D([], [], color=color, linestyle='dashed', label=f"{layer1} ↔ {layer2}", linewidth=2)
            legend_handles.append(handle)
    
        # Interlayer handle (shared style)
        interlayer_handle = mlines.Line2D([], [], color='dimgray', linestyle='dashed', label='Interlayer', linewidth=2)
        legend_handles.append(interlayer_handle)

        # Draw the legend and move it more to the right
        ax.legend(
            handles=legend_handles,
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            borderaxespad=0.,
            frameon=True,
            ncol=1,
            handlelength=2.5,
            fontsize=6
        )

        plt.title("MultiDNetwork Visualization")
        plt.axis('off')
        plt.tight_layout()

        if output_path:
            plt.savefig(output_path+'.'+format, bbox_inches='tight', format=format, dpi=dpi)
        else:
            plt.show()
