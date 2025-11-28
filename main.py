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
from agents.retrosynthesis import RetrosynthesisAgent
from agents.forward_prediction import ForwardPredictionAgent
from agents.route_scorer import RouteScorerAgent
from agents.safety_checker import SafetyCheckerAgent
from agents.protocol_generator import ProtocolGeneratorAgent
from agents.procurement import ProcurementAgent
from agents.data_curation import DataCurationAgent
from agents.feedback import FeedbackAgent

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
retro_agent = RetrosynthesisAgent()
prediction_agent = ForwardPredictionAgent()
scorer_agent = RouteScorerAgent()
safety_agent = SafetyCheckerAgent()
protocol_agent = ProtocolGeneratorAgent()
procurement_agent = ProcurementAgent()
curation_agent = DataCurationAgent()
feedback_agent = FeedbackAgent()

# Data Models
class UserQuery(BaseModel):
    query: str

class AgentResponse(BaseModel):
    intent: Dict[str, Any]
    molecule_info: Dict[str, Any]
    literature: List[Dict[str, Any]]
    reactions: List[Dict[str, Any]]
    routes: Optional[List[Dict[str, Any]]] = None
    protocol: Optional[str] = None
    error: Optional[str] = None

@app.post("/api/process", response_model=AgentResponse)
async def process_query(user_query: UserQuery):
    """
    Orchestrates the full agent workflow.
    """
    try:
        # 1. Parse Intent
        print(f"Processing query: {user_query.query}")
        intent = intent_agent.parse(user_query.query)
        
        target_name = intent.get("target_molecule")
        molecule_info = {}
        literature = []
        reactions = []
        routes = []
        protocol = "No suitable route found."

        if target_name:
            # 2. Canonicalize
            print(f"Canonicalizing: {target_name}")
            molecule_info = canonicalizer_agent.canonicalize(target_name)
            
            target_smiles = molecule_info.get("canonical_smiles") or molecule_info.get("smiles")
            
            # 3. Literature Search
            print(f"Searching literature for: {target_name}")
            literature = literature_agent.search(f"Synthesis of {target_name}")

            # 4. Reaction Search
            if target_smiles:
                print(f"Searching reactions for SMILES: {target_smiles}")
                reactions = reaction_agent.search(target_smiles)
            
            # 5. Retrosynthesis Planning
            print(f"Planning routes for: {target_name}")
            routes = retro_agent.plan(target_name)

            # 6-8. Process Routes (Prediction, Safety, Procurement, Scoring)
            for route in routes:
                # Forward Prediction for each step
                for step in route.get("steps", []):
                    prediction = prediction_agent.predict(step)
                    step["prediction"] = prediction
                
                # Safety Check
                safety_info = safety_agent.check(route)
                route["safety"] = safety_info
                
                # Procurement Estimation
                cost_info = procurement_agent.estimate(route)
                route["cost"] = cost_info

            # Rank Routes
            ranked_routes = scorer_agent.score(routes)
            routes = ranked_routes # Update routes with ranked version
            
            # 9. Generate Protocol for Top Route
            if ranked_routes:
                top_route = ranked_routes[0]
                protocol = protocol_agent.generate(top_route)

            # 11. Log Provenance
            curation_agent.log("process_query", {
                "query": user_query.query,
                "intent": intent,
                "top_route_id": ranked_routes[0].get("id") if ranked_routes else None
            })

        return AgentResponse(
            intent=intent,
            molecule_info=molecule_info,
            literature=literature,
            reactions=reactions,
            routes=routes,
            protocol=protocol
        )
        
    except Exception as e:
        print(f"Error: {e}")
        return AgentResponse(
            intent={},
            molecule_info={},
            literature=[],
            reactions=[],
            error=str(e)
        )

@app.get("/health")
async def health_check():
    return {"status": "ok"}

@app.post("/api/feedback")
async def submit_feedback(route_id: str, rating: int, comments: str):
    return feedback_agent.process_feedback(route_id, rating, comments)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
