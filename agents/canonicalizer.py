try:
    import pubchempy as pcp
except ImportError:
    pcp = None

try:
    from rdkit import Chem
except ImportError:
    Chem = None

from typing import Dict, Any, Optional

class CanonicalizerAgent:
    """
    Agent 2: Canonicalizer (Name -> Structure)
    Normalizes names to canonical identifiers (SMILES, InChI).
    """
    def __init__(self):
        pass

    def canonicalize(self, chemical_name: str) -> Dict[str, Any]:
        """
        Converts a chemical name to canonical structure info.
        """
        result = {
            "name": chemical_name,
            "smiles": None,
            "inchi": None,
            "molecular_weight": None,
            "formula": None,
            "error": None
        }

        try:
            # 1. Try PubChem Resolution
            if pcp:
                compounds = pcp.get_compounds(chemical_name, 'name')
                
                if compounds:
                    compound = compounds[0]
                    result["smiles"] = compound.isomeric_smiles
                    result["inchi"] = compound.inchi
                    result["molecular_weight"] = compound.molecular_weight
                    result["formula"] = compound.molecular_formula
                    
                    # 2. Validate with RDKit
                    if Chem:
                        mol = Chem.MolFromSmiles(result["smiles"])
                        if mol:
                            result["canonical_smiles"] = Chem.MolToSmiles(mol, canonical=True)
                        else:
                            result["error"] = "PubChem returned invalid SMILES"
                    else:
                         result["canonical_smiles"] = result["smiles"] # Fallback if RDKit missing
                else:
                    result["error"] = "Compound not found in PubChem"
            else:
                result["error"] = "PubChemPy not installed"
                # Mock response for demo
                if chemical_name.lower() == "aspirin":
                     result["smiles"] = "CC(=O)OC1=CC=CC=C1C(=O)O"
                     result["canonical_smiles"] = "CC(=O)OC1=CC=CC=C1C(=O)O"
                     result["error"] = None

        except Exception as e:
            result["error"] = str(e)

        return result

if __name__ == "__main__":
    agent = CanonicalizerAgent()
    print(agent.canonicalize("Aspirin"))
    print(agent.canonicalize("Ibuprofen"))
