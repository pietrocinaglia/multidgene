# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

###
# @project: multi-DGeNe
# @author: Pietro Cinaglia
# @mail: cinaglia@unicz.it
# @description: Multi-Dimensional Gene Network (multi-DGeNe) inference, bridging Co-Expression and Expression Dynamics across tissues.
# @url: https://github.com/pietrocinaglia/multidgene
###

# Own modules

class VisualDesigner:

    multidimensional_network = None

    def __init__(self, multidimensional_network, output_path:str=None, format='eps', dpi='figure'):
        self.multidimensional_network = multidimensional_network
        self.output_path = output_path
        self.format = format
        self.dpi = dpi
        self.title = "MultiDNetwork Visualization"

    def plot(self):
        pos = {}
        layer_positions = {}
        color_map = list(mcolors.TABLEAU_COLORS.values())
        layer_list = sorted(self.multidimensional_network.layers)
        spacing = 5
        # Create legend handles
        legend_handles = []

        fig, ax = plt.subplots()

        x_offset_step = 1.5  # horizontal offset between layers
        y_offset_step = 5    # vertical offset between layers

        for i, layer in enumerate(layer_list):
            x_offset = 0
            y_offset = i * y_offset_step
            if i %2 == 0:
                x_offset = x_offset_step
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            layer_pos = nx.spring_layout(self.multidimensional_network.graph.subgraph(layer_nodes), seed=42)
            for node, coord in layer_pos.items():
                # shift x and y coordinates
                pos[node] = (coord[0] + x_offset, coord[1] + y_offset)

            layer_positions[layer] = color_map[i % len(color_map)]
            
            # Draw layer boundary
            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            if x_vals and y_vals:
                min_x, max_x = min(x_vals), max(x_vals)
                min_y, max_y = min(y_vals), max(y_vals)
                width = (max_x - min_x + 1) + 0.0
                height = (max_y - min_y + 1) + 1.0
                
                ellipse = patches.Ellipse(((min_x + max_x) / 2, (min_y + max_y) / 2), width, height, color=layer_positions[layer], alpha=0.2)
                ax.add_patch(ellipse)
                #rectangle = patches.Rectangle((min_x - 0.5, min_y - 0.5), width, height, linewidth=0.8, edgecolor=layer_positions[layer], facecolor=layer_positions[layer], alpha=0.2)
                #ax.add_patch(rectangle)

                ax.text((min_x + max_x) / 2, max_y + 1.2, layer, ha='center', fontsize=8, weight='bold')

        # Draw nodes (all in black)
        nx.draw_networkx_nodes(self.multidimensional_network.graph, pos, node_color='grey', node_size=250, ax=ax)

        # Intralayer edges
        for layer in layer_list:
           # Draw intralayer edges
            layer_edges = [(u, v) for u, v, _ in self.multidimensional_network.graph.edges(data=True) if u[1] == layer and v[1] == layer]
            nx.draw_networkx_edges(self.multidimensional_network.graph, pos, edgelist=layer_edges,
                                edge_color=layer_positions[layer], ax=ax, arrows=False)
            # Legend for intralayer edges
            color = layer_positions[layer]
            handle = mlines.Line2D([], [], color=color, label=f"{layer}", linewidth=2)
            legend_handles.append(handle)

        # Interlayer edges with unique colors
        layer_pairs = set((u[1], v[1]) for u, v, _ in self.multidimensional_network.graph.edges(data=True) if u[1] != v[1])
        
        #interlayer_colors = list(mcolors.CSS4_COLORS.values())
        interlayer_colors = [c for c in mcolors.CSS4_COLORS.values() if mcolors.rgb_to_hsv(mcolors.to_rgb(c))[2] < 0.9]
        interlayer_colors = {pair: interlayer_colors[i % len(interlayer_colors)] for i, pair in enumerate(layer_pairs)}

        for (u, v, _) in self.multidimensional_network.graph.edges(data=True):
            if u[1] != v[1]:
                pair = (u[1], v[1])
                color = interlayer_colors.get(pair, 'grey')
                nx.draw_networkx_edges(self.multidimensional_network.graph, pos, edgelist=[(u, v)], edge_color=color,
                                    style='dashed', ax=ax, arrows=False)

        # Draw node labels
        nx.draw_networkx_labels(self.multidimensional_network.graph, pos, font_size=4, font_color="black",
                                labels={n: n[0] for n in self.multidimensional_network.graph.nodes}, ax=ax)

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

        plt.title(self.title)
        plt.axis('off')
        plt.tight_layout()

        if self.output_path:
            plt.savefig(self.output_path+'.'+format, bbox_inches='tight', format=self.format, dpi=self.dpi)
        else:
            plt.show()

    def plot3d(self,):

        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111, projection='3d')

        pos = {}
        layer_positions = {}
        color_map = list(mcolors.TABLEAU_COLORS.values())
        layer_list = sorted(self.multidimensional_network.layers)

        z_spacing = 10  # space between layers

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

            # Store bounding box coordinates
            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            if x_vals and y_vals:
                min_x, max_x = min(x_vals), max(x_vals)
                min_y, max_y = min(y_vals), max(y_vals)
                padding = 0.6
                layer_bounds[layer] = {
                    "x_min": min_x - padding, "x_max": max_x + padding,
                    "y_min": min_y - padding, "y_max": max_y + padding,
                    "z": z
                }

        # Draw nodes
        for node, (x, y, z) in pos.items():
            ax.scatter(x, y, z, color='grey', s=45)

        # Node labels
        for node, (x, y, z) in pos.items():
            ax.text(x, y, z + 0.5, node[0], fontsize=6)

        # Intralayer edges
        for layer in layer_list:
            color = layer_positions[layer]
            for u, v, _ in self.multidimensional_network.graph.edges(data=True):
                if u[1] == v[1] == layer:
                    x_vals = [pos[u][0], pos[v][0]]
                    y_vals = [pos[u][1], pos[v][1]]
                    z_vals = [pos[u][2], pos[v][2]]
                    ax.plot(x_vals, y_vals, z_vals, color=color, linewidth=1.5)

        # Interlayer edges
        interlayer_colors_raw = [c for c in mcolors.CSS4_COLORS.values() if mcolors.rgb_to_hsv(mcolors.to_rgb(c))[2] < 0.85]
        interlayer_pairs = list(set((u[1], v[1]) for u, v, _ in self.multidimensional_network.graph.edges(data=True) if u[1] != v[1]))
        interlayer_colors = {pair: interlayer_colors_raw[i % len(interlayer_colors_raw)] for i, pair in enumerate(interlayer_pairs)}

        for u, v, _ in self.multidimensional_network.graph.edges(data=True):
            if u[1] != v[1]:
                pair = (u[1], v[1])
                color = interlayer_colors.get(pair, 'dimgray')
                x_vals = [pos[u][0], pos[v][0]]
                y_vals = [pos[u][1], pos[v][1]]
                z_vals = [pos[u][2], pos[v][2]]
                ax.plot(x_vals, y_vals, z_vals, color=color, linestyle='dashed', linewidth=1)

        # Draw 3D colored boxes around each layer

        # Parallelepiped
        for layer, bounds in layer_bounds.items():
            x_min = bounds["x_min"]
            x_max = bounds["x_max"]
            y_min = bounds["y_min"]
            y_max = bounds["y_max"]
            z = bounds["z"]
            dz = 1.0  # height of the box

            # Vertices of the bottom and top face
            verts = [
                [x_min, y_min, z],
                [x_max, y_min, z],
                [x_max, y_max, z],
                [x_min, y_max, z],
                [x_min, y_min, z + dz],
                [x_max, y_min, z + dz],
                [x_max, y_max, z + dz],
                [x_min, y_max, z + dz]
            ]

            faces = [
                [verts[0], verts[1], verts[2], verts[3]],  # bottom
                [verts[4], verts[5], verts[6], verts[7]],  # top
                [verts[0], verts[1], verts[5], verts[4]],  # front
                [verts[1], verts[2], verts[6], verts[5]],  # right
                [verts[2], verts[3], verts[7], verts[6]],  # back
                [verts[3], verts[0], verts[4], verts[7]]   # left
            ]
        
            box = Poly3DCollection(faces, alpha=0.15, facecolor=layer_positions[layer], edgecolor=layer_positions[layer])
            ax.add_collection3d(box)

            # Layer label
            ax.text((x_min + x_max) / 2, y_max + 0.5, z + dz + 0.3, layer, ha='center', fontsize=10, weight='bold', color=layer_positions[layer])

        '''
        for layer, bounds in layer_bounds.items():
            x_min, x_max = bounds["x_min"], bounds["x_max"]
            y_min, y_max = bounds["y_min"], bounds["y_max"]
            z = bounds["z"]
            
            # Calculate center and radii
            x_center = (x_min + x_max) / 2
            y_center = (y_min + y_max) / 2
            z_center = z + 0.5  # center in z, raise half the dz of your previous box height
            
            rx = (x_max - x_min) / 2  # radius x-axis
            ry = (y_max - y_min) / 2  # radius y-axis
            rz = 0.5                  # radius z-axis (height/2)
            
            # Create a grid of angles
            u = np.linspace(0, 2 * np.pi, 40)
            v = np.linspace(0, np.pi, 20)
            
            # Parametric equations for the ellipsoid surface
            x = x_center + rx * np.outer(np.cos(u), np.sin(v))
            y = y_center + ry * np.outer(np.sin(u), np.sin(v))
            z_ellipsoid = z_center + rz * np.outer(np.ones_like(u), np.cos(v))
            
            # Plot the surface
            ax.plot_surface(x, y, z_ellipsoid, color=layer_positions[layer], alpha=0.15, edgecolor='k', linewidth=0.3)
            
            # Layer label above the ellipsoid
            ax.text(x_center, y_center, z_center + rz + 0.3, layer, ha='center', fontsize=10, weight='bold', color=layer_positions[layer])
        '''
        # Legend
        legend_handles = []
        for layer in layer_list:
            handle = mlines.Line2D([], [], color=layer_positions[layer], label=layer, linewidth=2)
            legend_handles.append(handle)

        for (l1, l2), color in interlayer_colors.items():
            handle = mlines.Line2D([], [], color=color, linestyle='dashed', label=f"{l1} ↔ {l2}", linewidth=2)
            legend_handles.append(handle)

        ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=8)

        ax.set_title(self.title, fontsize=13, weight='bold')
        ax.set_axis_off()
        ax.view_init(elev=20, azim=135)

        plt.tight_layout()
        if self.output_path:
            plt.savefig(self.output_path, bbox_inches='tight', format=self.format)
        else:
            plt.show()

