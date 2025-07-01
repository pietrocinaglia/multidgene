from flask import Flask, render_template, request, jsonify
import csv
import os

# Own packages
from multidgene import MultiDGeNe

WORKSPACE = os.path.dirname(os.path.realpath(__file__)) + "/"
app = Flask(__name__)

def load_tissues_from_csv(filepath=WORKSPACE+'data/tissues.csv'):
    tissues = []
    with open(filepath, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) >= 2:
                tissues.append({'label': row[0].strip(), 'value': row[1].strip(), 'group': row[2].strip})
    return tissues

@app.route('/')
def index():
    tissues = load_tissues_from_csv()
    return render_template('index.html', tissues=tissues)

@app.route('/generate-network', methods=['POST'])
def generate_network():
    data = request.json
    genes_input = data['genes']
    tissues_selected = data['tissues']
    p_threshold = float(data['p_threshold'])
    mode = data['mode']

    # Sanitize whitespace and comma-based to \n
    genes_input = genes_input.replace(" ", "")
    genes_input = genes_input.replace(",", "\n")

    genes = [g.strip() for g in genes_input.splitlines() if g.strip()]

    multidgene = MultiDGeNe(
        genes_of_interest=genes,
        tissues_of_interest=tissues_selected,
        p_threshold=p_threshold,
        mode=mode,
        data_path=WORKSPACE+'data/',
        verbose=False
    )

    network, _ = multidgene.construct_multidimensional_network()

    html_plot = network.generate_3d_plotly()
    
    return jsonify({'html_plot': html_plot})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=int(os.environ.get("PORT", 5000)), debug=True)
