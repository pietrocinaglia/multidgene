from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import os
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import uuid
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from your_module import run_multidgene_analysis  # Replace with actual import

app = FastAPI()

# CORS for frontend communication
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

@app.post("/api/run")
async def run_analysis(request: Request):
    data = await request.json()
    genes = data.get("genes", [])
    tissues = data.get("tissues", [])
    mode = data.get("mode", "equal_exp")

    # Run your core function, get NetworkX or edge list
    edge_list = run_multidgene_analysis(genes, tissues, mode)

    # Generate unique ID
    session_id = str(uuid.uuid4())

    # Save CSV
    csv_path = os.path.join(OUTPUT_DIR, f"{session_id}.csv")
    pd.DataFrame(edge_list, columns=["Source", "Target", "Layer", "Weight"]).to_csv(csv_path, index=False)

    # Save PNG visualization
    png_path = os.path.join(OUTPUT_DIR, f"{session_id}.png")

    G = nx.Graph()
    tissue_colors = {}
    nodes_by_color = {}

    # Extract unique tissues (layers)
    unique_layers = list({layer for _, _, layer, _ in edge_list})
    colormap = cm.get_cmap("tab10", len(unique_layers))  # Use tab10 for good color variety

    # Assign color to each layer
    for idx, layer in enumerate(unique_layers):
        hex_color = mcolors.to_hex(colormap(idx))
        tissue_colors[layer] = hex_color
        nodes_by_color[hex_color] = []

    for source, target, layer, weight in edge_list:
        n1 = f"{source} ({layer})"
        n2 = f"{target} ({layer})"
        G.add_edge(n1, n2, weight=weight)
        nodes_by_color[tissue_colors[layer]].extend([n1, n2])

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=42)  # consistent layout

    for color, nodes in nodes_by_color.items():
        nx.draw_networkx_nodes(G, pos, nodelist=list(set(nodes)), node_color=color, label=color)

    nx.draw_networkx_edges(G, pos, alpha=0.5)
    nx.draw_networkx_labels(G, pos, font_size=8)

    # Add legend with tissue names
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=tissue_colors[t], label=t) for t in tissue_colors]
    plt.legend(handles=legend_elements, title="Tissues", loc="lower left", bbox_to_anchor=(0, 0))

    plt.axis("off")
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.clf()


    return JSONResponse(content={"csv": f"/api/download/{session_id}.csv", "png": f"/api/download/{session_id}.png"})

@app.get("/api/download/{filename}")
async def download_file(filename: str):
    file_path = os.path.join(OUTPUT_DIR, filename)
    if os.path.exists(file_path):
        return FileResponse(path=file_path, filename=filename)
    return JSONResponse(status_code=404, content={"message": "File not found."})

@app.get("/api/network")
def get_network_json():
    json_edges = []
    for source, target, layer, weight in edge_list:
        json_edges.append({
            "source": f"{source} ({layer})",
            "target": f"{target} ({layer})",
            "layer": layer,
            "weight": weight
        })
    return JSONResponse(content=json_edges)
