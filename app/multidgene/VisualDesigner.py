# -*- coding: utf-8 -*-

import networkx as nx
import matplotlib.colors as mcolors

import plotly.graph_objects as go
import plotly.io as pio

import uuid


class VisualDesigner:

    def __init__(self, multidimensional_network, scroll_height_px:int=800, output_path: str = None, format='eps', dpi='figure'):
        self.multidimensional_network = multidimensional_network
        self.scroll_height_px = scroll_height_px
        self.output_path = output_path
        self.format = format
        self.dpi = dpi
        self.title = "MultiDNetwork Visualization"
        self.color_map = list(mcolors.TABLEAU_COLORS.values())


    def generate_clickable_html(self, fig):
        div_id = f"plot_{uuid.uuid4().hex}"

        plot_html = pio.to_html(fig, full_html=False, include_plotlyjs='cdn', div_id=div_id)

        custom_js = f"""
        <script>
        document.addEventListener('DOMContentLoaded', function () {{
            var gd = document.getElementById('{div_id}');
            gd.on('plotly_click', function (data) {{
                var clicked_point = data.points[0];
                if (!clicked_point || !clicked_point.text) return;

                var clicked_node = clicked_point.text;
                alert(clicked_node);

                //var all_data = gd.data;
                //console.log(all_data);
            }});
        }});
        </script>
        """

        return f"""
        <div style="height:{self.scroll_height_px}px; overflow-y:auto; border:1px solid #ccc;">
            {plot_html}
            {custom_js}
        </div>
        """


    def generate_3d_plotly(self):
        """
        Generate a 3D Plotly visualization of the multidimensional network.

        Args:
            scroll_height_px (int): Height of the scrollable HTML container in pixels.

        Returns:
            str: HTML string containing the Plotly plot.
        """
        layer_list = sorted(self.multidimensional_network.layers)

        layer_thickness = 0.1
        x_offset_step = 3
        y_offset_step = 0
        z_offset_step = 10

        # Layout
        layer_positions = {}
        pos = {}
        for i, layer in enumerate(layer_list):
            x_offset = x_offset_step if i % 2 == 0 else 0
            y_offset = i * y_offset_step
            z = i * z_offset_step

            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            layout = nx.spring_layout(self.multidimensional_network.graph.subgraph(layer_nodes), seed=42)

            for node, (x, y) in layout.items():
                pos[node] = (x + x_offset, y + y_offset, z)

            layer_positions[layer] = self.color_map[i % len(self.color_map)]

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
                x=x_vals,
                y=y_vals,
                z=z_vals,
                mode='lines',
                line=dict(color=color, width=4, dash=style),
                hoverinfo='none',
                showlegend=False
            ))

        # Layers
        layer_boxes = []
        for i, layer in enumerate(layer_list):
            layer_nodes = [n for n in self.multidimensional_network.graph.nodes if n[1] == layer]
            if not layer_nodes:
                continue

            x_vals = [pos[n][0] for n in layer_nodes]
            y_vals = [pos[n][1] for n in layer_nodes]
            z = i * z_offset_step

            x_min, x_max = min(x_vals) - 0.5, max(x_vals) + 0.5
            y_min, y_max = min(y_vals) - 0.5, max(y_vals) + 0.5

            color = layer_positions[layer]

            # Define vertices for bottom and top faces
            vertices_x = [x_min, x_max, x_max, x_min, x_min, x_max, x_max, x_min]
            vertices_y = [y_min, y_min, y_max, y_max, y_min, y_min, y_max, y_max]
            vertices_z = [z-layer_thickness] * 4 + [z+layer_thickness] * 4

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
                x=[None], y=[None], z=[None],
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

        '''
        plot_html = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')
        scrollable_html = f"""
        <div style="height:{scroll_height_px}px; overflow-y:auto; border:1px solid #ccc;">
            {plot_html}
        </div>
        """

        if self.output_path:
            with open(self.output_path, "w", encoding="utf-8") as f:
                f.write(scrollable_html)
        '''

        scrollable_html = self.generate_clickable_html(fig)

        if self.output_path:
            with open(self.output_path, "w", encoding="utf-8") as f:
                f.write(scrollable_html)

        return scrollable_html
        