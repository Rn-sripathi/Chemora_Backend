try:
    import pubchempy as pcp
except ImportError:
    pcp = None

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    Chem = None
    Descriptors = None
    rdMolDescriptors = None

from typing import Dict, Any, Optional

try:
    import requests
except ImportError:
    requests = None

class CanonicalizerAgent:
    """
    Agent 2: Canonicalizer (Name -> Structure)
    Normalizes names to canonical identifiers (SMILES, InChI).
    Prioritizes OPSIN -> PubChem -> RDKit.
    """
    def __init__(self):
        pass

    def fetch_from_opsin(self, name: str) -> Optional[str]:
        """
        Fetches SMILES from OPSIN API.
        """
        if not requests:
            return None
        try:
            url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
            response = requests.get(url, timeout=5)
            if response.status_code == 200:
                data = response.json()
                return data.get("smiles")
        except Exception as e:
            print(f"OPSIN Error: {e}")
        return None

    def canonicalize(self, chemical_name: str) -> Dict[str, Any]:
        """
        Converts a chemical name to canonical structure info.
        """
        result = {
            "name": chemical_name,
            "smiles": None,
            "canonical_smiles": None,
            "inchi": None,
            "molecular_weight": None,
            "formula": None,
            "error": None
        }

        try:
            # 1. Try OPSIN (Primary Source)
            smiles = self.fetch_from_opsin(chemical_name)
            source = "OPSIN" if smiles else None

            # 2. Try PubChem (Fallback)
            if not smiles and pcp:
                try:
                    compounds = pcp.get_compounds(chemical_name, 'name')
                    if compounds:
                        compound = compounds[0]
                        smiles = compound.isomeric_smiles
                        # Populate metadata from PubChem if available
                        result["molecular_weight"] = compound.molecular_weight
                        result["formula"] = compound.molecular_formula
                        source = "PubChem"
                except Exception as e:
                    print(f"PubChem Error: {e}")

            # 3. Process with RDKit (Canonicalization)
            if smiles:
                result["smiles"] = smiles
                if Chem:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            # Standardize to Canonical SMILES
                            result["canonical_smiles"] = Chem.MolToSmiles(mol, isomericSmiles=True)
                            result["inchi"] = Chem.MolToInchi(mol)
                            
                            # Calculate MW/Formula if missing (e.g. from OPSIN)
                            if not result["molecular_weight"]:
                                result["molecular_weight"] = round(Chem.Descriptors.MolWt(mol), 2)
                            if not result["formula"]:
                                result["formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
                        else:
                            result["error"] = f"Invalid SMILES from {source}"
                    except Exception as e:
                         result["error"] = f"RDKit Error: {e}"
                         result["canonical_smiles"] = smiles # Fallback
                else:
                    result["canonical_smiles"] = smiles
                    result["error"] = "RDKit not installed"
            else:
                result["error"] = "Compound not found in OPSIN or PubChem"

        except Exception as e:
            result["error"] = str(e)

        return result

if __name__ == "__main__":
    agent = CanonicalizerAgent()
    print(agent.canonicalize("Aspirin"))
    print(agent.canonicalize("Ibuprofen"))
