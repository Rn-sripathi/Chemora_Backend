import json
import sys
import os

# Add project root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from agents.intent_parser import IntentParserAgent
from agents.canonicalizer import CanonicalizerAgent
from agents.literature_retrieval import LiteratureRetrievalAgent
from agents.reaction_retrieval import ReactionRetrievalAgent

def main():
    print("=== Chemical Agent System Demo ===")
    
    # 1. User Intent
    print("\n--- Agent 1: User Intent ---")
    intent_agent = IntentParserAgent()
    
    if len(sys.argv) > 1:
        user_query = " ".join(sys.argv[1:])
    else:
        print("Enter your query (e.g., 'I need the synthesis of Aspirin'):")
        user_query = input("> ")
        if not user_query:
            user_query = "I need the synthesis of Aspirin. It should be cheap and safe."
            print(f"No input provided. Using default: {user_query}")

    intent = intent_agent.parse(user_query)
    print(f"User Query: {user_query}")
    print(f"Parsed Intent: {json.dumps(intent, indent=2)}")
    
    target_name = intent.get("target")
    if not target_name:
        print("No target found.")
        return

    # 2. Canonicalizer
    print("\n--- Agent 2: Canonicalizer ---")
    canon_agent = CanonicalizerAgent()
    canon_result = canon_agent.canonicalize(target_name)
    print(f"Canonicalizing '{target_name}':")
    print(json.dumps(canon_result, indent=2))
    
    target_smiles = canon_result.get("canonical_smiles") or canon_result.get("smiles")
    if not target_smiles:
        print("Could not resolve SMILES.")
        # Fallback for demo if API fails
        target_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
        print(f"Using fallback SMILES: {target_smiles}")

    # 3. Literature Retrieval
    print("\n--- Agent 3: Literature Retrieval ---")
    lit_agent = LiteratureRetrievalAgent()
    lit_results = lit_agent.search(f"Synthesis of {target_name}")
    print(f"Found {len(lit_results)} documents:")
    for doc in lit_results:
        print(f"- {doc['text']} (Yield: {doc['yield']})")

    # 4. Reaction Retrieval
    print("\n--- Agent 4: Reaction Retrieval ---")
    rxn_agent = ReactionRetrievalAgent()
    rxn_results = rxn_agent.search(target_smiles)
    print(f"Found {len(rxn_results)} reactions:")
    for rxn in rxn_results:
        print(f"- {rxn['name']}: {rxn['smiles']} (Similarity: {rxn['similarity']:.2f})")

if __name__ == "__main__":
    main()
