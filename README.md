# DILIPredictor Online tools and Website - Drug-Induced Liver Injury Prediction 

## Overview
DILI is a comprehensive repository aimed at enhancing the early detection of Drug-Induced Liver Injury (DILI) through the integration of predicted in vivo and in vitro data. This project utilizes advanced machine learning models and chemical informatics to predict the likelihood of DILI for various compounds.
For code see: https://github.com/srijitseal/DILI

### Using PyPI 
You can install the DILI Predictor using pip: `pip install dilipred`. ***Please use Python <3.12, >=3.9***

### Building from Source
You can also build from source using python-poetry:
1. Clone the repository: `git clone https://github.com/Manas02/dili-pip.git`
2. Navigate to the project directory: `cd dili-pip/`
3. Install dependencies: `poetry install`
4. Activate the virtual environment: `poetry shell`
5. Build the project: `poetry build`

## Usage

### Running DILIPredictor as CLI
To get started with the CLI, use: `dili -h`

### Predicting DILI for a Single Molecule
Select from the sidebar to predict DILI for a single molecule.

### Running DILIPredictor as Library
Here's a basic example of how to use DILIPredictor as a Python library:

```python
from dilipred import DILIPRedictor

if __name__ == '__main__':
    dp = DILIPRedictor()
    smiles = "CCCCCCCO"
    result = dp.predict(smiles)
    print(result)
```

## Using local implementation
Download key files from https://github.com/srijitseal/DILI/raw/main/local_implementation.zip and run locally!

## Run Online on Server #1 v4.0.0 [Recommended]
If you prefer to use the predictor online via Uppsala University SciLifeLab Serve: https://dili.serve.scilifelab.se/

## Run Online on Server #2 v4.0.0
If you prefer to use the predictor online via streamlit: https://dilipredictor.streamlit.app/

## Citation
If you use DILI Predictor in your work, please cite:

Improved Early Detection of Drug-Induced Liver Injury by Integrating Predicted in vivo and in vitro Data; Srijit Seal, Dominic P. Williams, Layla Hosseini-Gerami, Manas Mahale, Anne E. Carpenter, Ola Spjuth, Andreas Bender bioRxiv 2024.01.10.575128; doi: https://doi.org/10.1101/2024.01.10.575128

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements
Developed and maintained by Srijit Seal and contributors.

## Contact
For any questions or issues, please open an issue on the GitHub repository.
