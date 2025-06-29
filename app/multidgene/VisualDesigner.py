# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import plotly.graph_objects as go

class VisualDesigner:

    def __init__(self, multidimensional_network, output_path: str = None, format='eps', dpi='figure'):
        self.multidimensional_network = multidimensional_network
        self.output_path = output_path
        self.format = format
        self.dpi = dpi
        self.title = "MultiDNetwork Visualization"
        self.color_map = list(mcolors.TABLEAU_COLORS.values())

    def _compute_positions(self, layer_nodes, offset=(0, 0), seed=42):
        layout = nx.spring_layout(self.multidimensional_network.graph.subgraph(layer_nodes), seed=seed)
        return {node: (coord[0] + offset[0], coord[1] + offset[1]) for node, coord in layout.items()}

    def _draw_layer_boundary(self, ax, layer_nodes, pos, color, label):
        if not layer_nodes:
            return

        x_vals = [pos[n][0] for n in layer_nodes]
        y_vals = [pos[n][1] for n in layer_nodes]
        min_x, max_x = min(x_vals), max(x_vals)
        min_y, max_y = min(y_vals), max(y_vals)
        width, height = max_x - min_x + 1.0, max_y - min_y + 2.0

        ellipse = patches.Ellipse(((min_x + max_x) / 2, (min_y + max_y) / 2), width, height, color=color, alpha=0.2)
        ax.add_patch(ellipse)
        ax.text((min_x + max_x) / 2, max_y + 1.2, label, ha='center', fontsize=12, weight='bold')

    def _get_edges_by_layer(self, layer):
        return [(u, v) for u, v, _ in self.multidimensional_network.graph.edges(data=True) if u[1] == v[1] == layer]

    def _get_interlayer_edges(self):
        return [(u, v) for u, v, _ in self.multidimensional_network.graph.edges(data=True) if u[1] != v[1]]

    def _draw_edges(self, ax, pos, layer_positions):
        legend_handles = []

        # Draw intralayer edges
        for layer in sorted(self.multidimensional_network.layers):
            color = layer_positions[layer]
            layer_edges = self._get_edges_by_layer(layer)
            nx.draw_networkx_edges(self.multidimensional_network.graph, pos, edgelist=layer_edges, edge_color=color, ax=ax, arrows=False)
            legend_handles.append(mlines.Line2D([], [], color=color, label=layer, linewidth=4))

        # Draw interlayer edges
        layer_pairs = set((u[1], v[1]) for u, v in self._get_interlayer_edges())
        interlayer_colors = {pair: color for pair, color in zip(layer_pairs, mcolors.CSS4_COLORS.values())}

        for u, v in self._get_interlayer_edges():
            pair = (u[1], v[1])
            color = interlayer_colors.get(pair, 'grey')
            nx.draw_networkx_edges(self.multidimensional_network.graph, pos, edgelist=[(u, v)], edge_color=color, style='dashed', ax=ax, arrows=False)

        for (layer1, layer2), color in interlayer_colors.items():
            legend_handles.append(mlines.Line2D([], [], color=color, linestyle='dashed', label=f"{layer1} ↔ {layer2}", linewidth=4))

        return legend_handles

    def _draw_nodes_and_labels(self, ax, pos):
        nx.draw_networkx_nodes(self.multidimensional_network.graph, pos, node_color='grey', node_size=350, ax=ax)
        nx.draw_networkx_labels(self.multidimensional_network.graph, pos, font_size=8, font_color="black",
                                labels={n: n[0] for n in self.multidimensional_network.graph.nodes}, ax=ax)

    def _compute_layer_positions(self):
        pos, layer_positions = {}, {}
        for i, layer in enumerate(sorted(self.multidimensional_network.layers)):
            offset = (1.5 if i % 2 == 0 else 0, i * 5)
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            layer_pos = self._compute_positions(layer_nodes, offset)
            pos.update(layer_pos)
            layer_positions[layer] = self.color_map[i % len(self.color_map)]
        return pos, layer_positions
    
    def plot(self):
        fig, ax = plt.subplots()
        pos, layer_positions = self._compute_layer_positions()

        for layer in sorted(self.multidimensional_network.layers):
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            self._draw_layer_boundary(ax, layer_nodes, pos, layer_positions[layer], layer)

        self._draw_nodes_and_labels(ax, pos)
        
        legend_handles = self._draw_edges(ax, pos, layer_positions)

        ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0., frameon=True, fontsize=6)
        ax.set_title(self.title)
        ax.axis('off')
        plt.tight_layout()

        if self.output_path:
            plt.savefig(f"{self.output_path}.{self.format}", bbox_inches='tight', format=self.format, dpi=self.dpi)
        else:
            plt.show()


    def plot3d(self):
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        pos = {}
        layer_positions = {}
        color_map = list(mcolors.TABLEAU_COLORS.values())
        layer_list = sorted(self.multidimensional_network.layers)
        z_spacing = 10

        layer_bounds = {}
        x_offset_step = 2.5
        y_offset_step = 0

        for i, layer in enumerate(layer_list):
            z = i * z_spacing
            x_offset = x_offset_step if i % 2 == 0 else 0
            y_offset = i * y_offset_step

            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            layer_pos = nx.spring_layout(self.multidimensional_network.graph.subgraph(layer_nodes), seed=42)

            for node, (x, y) in layer_pos.items():
                pos[node] = (x + x_offset, y + y_offset, z)

            layer_positions[layer] = color_map[i % len(color_map)]

            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            if x_vals and y_vals:
                padding = 0.6
                layer_bounds[layer] = {
                    "x_min": min(x_vals) - padding,
                    "x_max": max(x_vals) + padding,
                    "y_min": min(y_vals) - padding,
                    "y_max": max(y_vals) + padding,
                    "z": z
                }

        # Draw nodes
        for node, (x, y, z) in pos.items():
            ax.scatter(x, y, z, color='black', s=55)
            ax.text(x, y, z + 0.5, node[0], fontsize=6)

        # Intralayer edges
        for layer in layer_list:
            color = layer_positions[layer]
            edges = self._get_edges_by_layer(layer)
            for u, v in edges:
                x_vals = [pos[u][0], pos[v][0]]
                y_vals = [pos[u][1], pos[v][1]]
                z_vals = [pos[u][2], pos[v][2]]
                ax.plot(x_vals, y_vals, z_vals, color=color, linewidth=1.5)

        # Interlayer edges
        interlayer_edges = self._get_interlayer_edges()
        interlayer_colors_raw = [
            c for c in mcolors.CSS4_COLORS.values()
            if mcolors.rgb_to_hsv(mcolors.to_rgb(c))[2] < 0.85
        ]
        interlayer_pairs = list(set((u[1], v[1]) for u, v in interlayer_edges))
        interlayer_colors = {
            pair: interlayer_colors_raw[i % len(interlayer_colors_raw)]
            for i, pair in enumerate(interlayer_pairs)
        }

        for u, v in interlayer_edges:
            pair = (u[1], v[1])
            color = interlayer_colors.get(pair, 'dimgray')
            x_vals = [pos[u][0], pos[v][0]]
            y_vals = [pos[u][1], pos[v][1]]
            z_vals = [pos[u][2], pos[v][2]]
            ax.plot(x_vals, y_vals, z_vals, color=color, linestyle='dashed', linewidth=1)

        # 3D layer boxes
        for layer, bounds in layer_bounds.items():
            x_min, x_max = bounds["x_min"], bounds["x_max"]
            y_min, y_max = bounds["y_min"], bounds["y_max"]
            z = bounds["z"]
            dz = 1.0

            verts = [
                [x_min, y_min, z], [x_max, y_min, z], [x_max, y_max, z], [x_min, y_max, z],
                [x_min, y_min, z + dz], [x_max, y_min, z + dz],
                [x_max, y_max, z + dz], [x_min, y_max, z + dz]
            ]
            faces = [
                [verts[0], verts[1], verts[2], verts[3]],
                [verts[4], verts[5], verts[6], verts[7]],
                [verts[0], verts[1], verts[5], verts[4]],
                [verts[1], verts[2], verts[6], verts[5]],
                [verts[2], verts[3], verts[7], verts[6]],
                [verts[3], verts[0], verts[4], verts[7]]
            ]
            box = Poly3DCollection(faces, alpha=0.15, facecolor=layer_positions[layer], edgecolor=layer_positions[layer])
            ax.add_collection3d(box)

            ax.text((x_min + x_max) / 2, y_max + 0.3, z + dz + 0.3, layer, ha='center', fontsize=10, weight='bold', color=layer_positions[layer])

        # Legend
        legend_handles = [mlines.Line2D([], [], color=layer_positions[layer], label=layer, linewidth=2) for layer in layer_list]
        legend_handles += [
            mlines.Line2D([], [], color=color, linestyle='dashed', label=f"{l1} ↔ {l2}", linewidth=2)
            for (l1, l2), color in interlayer_colors.items()
        ]
        ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=8)

        ax.set_title(self.title, fontsize=13, weight='bold')
        ax.set_axis_off()
        ax.view_init(elev=20, azim=135)

        plt.tight_layout()
        if self.output_path:
            plt.savefig(self.output_path, bbox_inches='tight', format=self.format)
        else:
            plt.show()


    def generate_3d_plotly(self):
        pos = {}
        color_map = list(mcolors.TABLEAU_COLORS.values())
        layer_list = sorted(self.multidimensional_network.layers)
        z_spacing = 10
        layer_positions = {}

        x_offset_step = 3
        y_offset_step = 0

        # Layout
        for i, layer in enumerate(layer_list):
            z = i * z_spacing
            x_offset = x_offset_step if i % 2 == 0 else 0
            y_offset = i * y_offset_step
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            layout = nx.spring_layout(self.multidimensional_network.graph.subgraph(layer_nodes), seed=42)

            for node, (x, y) in layout.items():
                pos[node] = (x + x_offset, y + y_offset, z)

            layer_positions[layer] = color_map[i % len(color_map)]

        # Nodes
        node_x, node_y, node_z, node_text = [], [], [], []
        for node, (x, y, z) in pos.items():
            node_x.append(x)
            node_y.append(y)
            node_z.append(z)
            node_text.append(node[0])

        node_trace = go.Scatter3d(
            x=node_x,
            y=node_y,
            z=node_z,
            mode='markers+text',
            text=node_text,
            textposition='top center',
            marker=dict(size=10, color='#7f7f7f'),
            hoverinfo='text',
            name='Nodes'
        )

        layer_names = []
        layer_x = []
        layer_y = []
        layer_z = []

        for layer in sorted(set(layer for (_, layer) in pos)):
            xs = [x for (n, l), (x, y, z) in pos.items() if l == layer]
            ys = [y for (n, l), (x, y, z) in pos.items() if l == layer]
            zs = [z for (n, l), (x, y, z) in pos.items() if l == layer]
            if xs and ys and zs:
                layer_x.append(min(xs) - 0.2)  # adjust position as needed
                layer_y.append(sum(ys) / len(ys))
                layer_z.append(sum(zs) / len(zs))
                layer_names.append(f"{layer}")

        layer_label_traces = []
        for i, layer in enumerate(layer_names):
            layer_label_traces.append(go.Scatter3d(
                x=[layer_x[i]], y=[layer_y[i]], z=[layer_z[i]],
                mode='text',
                text=[layer],
                textposition='middle right',
                textfont=dict(size=16, color=layer_positions[layer], family="Arial", weight='bold'),
                hoverinfo='none',
                showlegend=False
            ))

        # Edges
        interlayer_colors = {}
        interlayer_pairs = list(set((u[1], v[1]) for u, v in self.multidimensional_network.graph.edges() if u[1] != v[1]))
        interlayer_palette = [c for c in mcolors.CSS4_COLORS.values() if mcolors.rgb_to_hsv(mcolors.to_rgb(c))[2] < 0.85]
        for i, pair in enumerate(interlayer_pairs):
            interlayer_colors[pair] = interlayer_palette[i % len(interlayer_palette)]

        edge_traces = []
        for u, v in self.multidimensional_network.graph.edges():
            x_vals = [pos[u][0], pos[v][0], None]
            y_vals = [pos[u][1], pos[v][1], None]
            z_vals = [pos[u][2], pos[v][2], None]
            if u[1] == v[1]:  # Intralayer
                color = layer_positions[u[1]]
                style = 'solid'
            else:  # Interlayer
                pair = (u[1], v[1])
                rev_pair = (v[1], u[1])
                color = interlayer_colors.get(pair, interlayer_colors.get(rev_pair, 'dimgray'))
                style = 'dash'
            edge_traces.append(go.Scatter3d(
                x=x_vals, y=y_vals, z=z_vals,
                mode='lines',
                line=dict(color=color, width=4, dash=style),
                hoverinfo='none',
                showlegend=False
            ))

        # Layers
        layer_boxes = []
        dz = 0.5  # thickness of the layer box
        extra_space = 3
        for i, layer in enumerate(layer_list):
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            if not layer_nodes:
                continue

            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            z = i * z_spacing

            x_min, x_max = min(x_vals) - 0.5, max(x_vals) + 0.5
            y_min, y_max = min(y_vals) - 0.5, max(y_vals) + 0.5

            color = layer_positions[layer]

            # Define vertices for bottom and top faces
            vertices_x = [x_min, x_max, x_max, x_min, x_min, x_max, x_max, x_min]
            vertices_y = [y_min, y_min, y_max, y_max, y_min, y_min, y_max, y_max]
            vertices_z = [z-dz] * 4 + [z+dz] * 4

            # Define the triangles composing the 6 faces of the box (12 triangles)
            i_faces = [
                0, 1, 2,  0, 2, 3,  # bottom
                4, 5, 6,  4, 6, 7,  # top
                0, 1, 5,  0, 5, 4,  # front
                1, 2, 6,  1, 6, 5,  # right
                2, 3, 7,  2, 7, 6,  # back
                3, 0, 4,  3, 4, 7   # left
            ]

            layer_boxes.append(go.Mesh3d(
                x=vertices_x,
                y=vertices_y,
                z=vertices_z,
                color=color,
                opacity=0.15,
                i=i_faces[0::3],
                j=i_faces[1::3],
                k=i_faces[2::3],
                hoverinfo='skip',
                showscale=False,
                name=f"Layer {layer}",
                showlegend=False
            ))
        
        # Layer legend markers
        legend_layer_traces = []
        for layer in layer_list:
            legend_layer_traces.append(go.Scatter3d(
                x=[None], y=[None], z=[None],  # no points displayed
                mode='markers',
                marker=dict(size=12, color=layer_positions[layer]),
                name=f"Layer {layer}",
                showlegend=True
            ))

        # Interlayer edges legend traces
        legend_interlayer_traces = []
        for (l1, l2), color in interlayer_colors.items():
            legend_interlayer_traces.append(go.Scatter3d(
                x=[None, None], y=[None, None], z=[None, None],
                mode='lines',
                line=dict(color=color, width=4, dash='dash'),
                name=f"Interlayer {l1} ↔ {l2}",
                showlegend=True
            ))

        # Combine all traces
        data = layer_boxes + layer_label_traces + edge_traces + [node_trace] + legend_layer_traces + legend_interlayer_traces

        fig = go.Figure(data=data)

        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
            ),
            title=self.title,
            margin=dict(l=0, r=0, b=0, t=40),
        )

        return fig
