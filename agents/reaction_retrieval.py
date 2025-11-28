try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
except ImportError:
    Chem = None
    AllChem = None
    DataStructs = None

from typing import List, Dict, Any

class ReactionRetrievalAgent:
    """
    Agent 4: Reaction Retrieval
    Queries reaction databases for precedent transformations matching target or substructures.
    """
    def __init__(self):
        # Mock Reaction Database
        # In reality, this would be a connection to a SQL or Vector DB
        self.reaction_db = [
            {"id": "rxn1", "smiles": "CC(=O)O.OC1=CC=CC=C1C(=O)O>>CC(=O)OC1=CC=CC=C1C(=O)O", "name": "Aspirin Synthesis"},
            {"id": "rxn2", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O>>CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "name": "Ibuprofen Synthesis (Dummy)"}, # Dummy reaction
        ]
        
        # Precompute fingerprints for the product of the reaction
        self.db_fingerprints = []
        if Chem and AllChem:
            for rxn in self.reaction_db:
                try:
                    # Extract product (right side of >>)
                    product_smiles = rxn["smiles"].split(">>")[-1]
                    mol = Chem.MolFromSmiles(product_smiles)
                    if mol:
                        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                        self.db_fingerprints.append((rxn, fp))
                except:
                    pass
        else:
            # Mock fingerprints if RDKit missing
            pass

    def search(self, target_smiles: str, similarity_threshold: float = 0.5) -> List[Dict[str, Any]]:
        """
        Finds reactions that produce similar molecules to the target.
        """
        if not (Chem and AllChem and DataStructs):
             # Mock search if RDKit missing
             return [r for r in self.reaction_db if target_smiles in r["smiles"]]

        target_mol = Chem.MolFromSmiles(target_smiles)
        if not target_mol:
            return [{"error": "Invalid target SMILES"}]

        target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, 2, nBits=2048)
        
        results = []
        for rxn, fp in self.db_fingerprints:
            similarity = DataStructs.TanimotoSimilarity(target_fp, fp)
            if similarity >= similarity_threshold:
                result = rxn.copy()
                result["similarity"] = similarity
                results.append(result)
        
        # Sort by similarity
        results.sort(key=lambda x: x["similarity"], reverse=True)
        return results

if __name__ == "__main__":
    agent = ReactionRetrievalAgent()
    # Search for Aspirin
    print(agent.search("CC(=O)OC1=CC=CC=C1C(=O)O"))
