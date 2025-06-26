# Multi-Dimensional Gene Network (multi-DGeNe)

**multi-DGeNe** is a tool for inferring Multi-Dimensional Gene Network (multi-DGeNe) inference based on Co-Expression and Expression Dynamics Across Tissues.

## Additional Information

Genotype-Tissue Expression (GTEx) data was retrieved from GTEx Portal (https://gtexportal.org) via its own API.

The human tissues are registered by GTEx Portal and must be recalled using the labels indicated by them as identification.

## Table of Contents

- [Getting Started](#getting-started)
- [Usage](#usage)


## Getting Started <a name="getting-started"></a>

To begin using **multi-DGeNe** you have to follow these steps:

1. Install dependecies

```
pip install -r requirements.txt
```

2. Copy **multidgene** and import it into your code:

- 2.1. Copy "multidgene" folder into your workspace.

- 2.2. Import multi-DGeNe in your code:

```
from multidgene import MultiDGene
```

3. Create an instance of **multi-DGeNe**:

```
multidgene = MultiDGene(genes_of_interest: list, tissues_of_interest: list = [], p_threshold: float = 0.05, data_path: str = None, mode: str = 'equal_exp', verbose: bool = False)
```

'mode: str' concerns the method applied for infering interlayer edges, and it can be defined as:
- 'equal_exp' (default): interlayer edges refer to equal expressed genes between tissues (layers)
- 'diff_exp': interlayer edges refer to differential expressed genes between tissues (layers)
- 'diff_co_exp': interlayer edges refer to differential co-expressed pairs of genes between tissues (layers)


4. Use the following methods:

- construct_multidimensional_network(self) -> tuple[pymnet.MultilayerNetwork, set]:
- plot(self, multidimensional_network: pymnet.MultilayerNetwork, edge_colors: set = {}, output_path: str = None):
- save_as_files(self, multidimensional_network: pymnet.MultilayerNetwork, output_path: str):


## Usage <a name="usage"></a>

See **test.py** for illustrative purpose, it implements our solution for a real-world use case.


### License

MIT License - Copyright (c) Pietro Cinaglia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.