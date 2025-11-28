from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import uvicorn
import os
import sys

# Add project root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from agents.intent_parser import IntentParserAgent
from agents.canonicalizer import CanonicalizerAgent
from agents.literature_retrieval import LiteratureRetrievalAgent
from agents.reaction_retrieval import ReactionRetrievalAgent

from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(title="Chemical Agent Platform API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # In production, replace with specific origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Initialize Agents
intent_agent = IntentParserAgent()
canonicalizer_agent = CanonicalizerAgent()
literature_agent = LiteratureRetrievalAgent()
reaction_agent = ReactionRetrievalAgent()

# Data Models
class UserQuery(BaseModel):
    query: str

class AgentResponse(BaseModel):
    intent: Dict[str, Any]
    canonical_data: Optional[Dict[str, Any]] = None
    literature_results: Optional[List[Dict[str, Any]]] = None
    reaction_results: Optional[List[Dict[str, Any]]] = None
    error: Optional[str] = None

@app.post("/api/process", response_model=AgentResponse)
async def process_query(user_query: UserQuery):
    """
    Orchestrates the full agent workflow.
    """
    response = AgentResponse(intent={})
    
    try:
        # 1. Parse Intent
        print(f"Processing query: {user_query.query}")
        intent = intent_agent.parse(user_query.query)
        response.intent = intent
        
        target_name = intent.get("target")
        if not target_name:
            return response

        # 2. Canonicalize
        print(f"Canonicalizing: {target_name}")
        canon_result = canonicalizer_agent.canonicalize(target_name)
        response.canonical_data = canon_result
        
        target_smiles = canon_result.get("canonical_smiles") or canon_result.get("smiles")
        
        # 3. Literature Search
        print(f"Searching literature for: {target_name}")
        response.literature_results = literature_agent.search(f"Synthesis of {target_name}")

        # 4. Reaction Search
        if target_smiles:
            print(f"Searching reactions for SMILES: {target_smiles}")
            response.reaction_results = reaction_agent.search(target_smiles)
        
    except Exception as e:
        response.error = str(e)
        print(f"Error: {e}")

    return response

@app.get("/health")
async def health_check():
    return {"status": "ok"}

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
