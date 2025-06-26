from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
import os
import json

# Own modules
from multidgene import MultiDGene
from multidgene import MultiDNetwork

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
app = FastAPI()

# CORS for frontend communication
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
async def root():
    return {"message": "Multi-Dimensional Gene Network (multi-DGeNe) inference. bridging Co-Expression and Expression Dynamics across tissues."}

'''
Example of body: {"genes":["BRCA1","MZT2B","TP53","DDB2","USP24","ANAPC15","PSMD14","COP1"],"tissues":["Brain_Cortex","Brain_Hippocampus","Liver","Thyroid"]}
'''
@app.post("/construct")
async def construct_multidimensional_network(request: Request):
    data = await request.json()
    genes = data.get("genes", [])
    tissues = data.get("tissues", [])
    mode = data.get("mode", "equal_exp")

    # Run your core function, get NetworkX or edge list
    multidgene = MultiDGene(genes_of_interest=genes, tissues_of_interest=tissues, mode=mode, p_threshold=0.05, data_path=WORKSPACE+'../data/')
    multidimensional_network, _ = multidgene.construct_multidimensional_network()

    edge_list = multidimensional_network.get_edges()

    output = json.dumps(edge_list)

    return output