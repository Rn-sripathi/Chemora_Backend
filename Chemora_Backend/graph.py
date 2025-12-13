from typing import TypedDict, Annotated, List, Dict, Any, Optional
from langgraph.graph import StateGraph, END
import asyncio

# Import Agents
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

# Define State
class AgentState(TypedDict):
    user_query: str
    image_path: Optional[str]
    intent: Dict[str, Any]
    target_molecule: Optional[str]
    molecule_data: Dict[str, Any]
    literature: List[Dict[str, Any]]
    reactions: List[Dict[str, Any]]
    routes: List[Dict[str, Any]]
    protocol: Optional[str]
    routes: List[Dict[str, Any]]
    protocol: Optional[str]
    errors: List[str]
    visited_agents: List[str]

# Define Nodes
async def node_intent_parser(state: AgentState):
    print("--- Intent Parser ---")
    intent = intent_agent.parse(state["user_query"], state.get("image_path"))
    return {
        "intent": intent,
        "target_molecule": intent.get("target"),
        "visited_agents": ["Intent Parser"]
    }

async def node_canonicalizer(state: AgentState):
    print("--- Canonicalizer ---")
    target = state["target_molecule"]
    if not target:
        return {"errors": ["No target molecule identified"]}
    
    data = canonicalizer_agent.canonicalize(target)
    return {
        "molecule_data": data,
        "visited_agents": ["Canonicalizer"]
    }

async def node_retrieval(state: AgentState):
    print("--- Retrieval (Lit & Rxn) ---")
    target = state["target_molecule"]
    smiles = state["molecule_data"].get("canonical_smiles") or state["molecule_data"].get("smiles")
    molecule_name = state["molecule_data"].get("name", "")
    
    # Run in parallel
    async def run_lit():
        # Call retrieve() method instead of search()
        return literature_agent.retrieve(f"Synthesis of {target}", molecule_name)
    
    async def run_rxn():
        if smiles:
            # Call find_similar() method instead of search()
            return reaction_agent.find_similar(smiles, molecule_name)
        return []

    lit, rxn = await asyncio.gather(run_lit(), run_rxn())
    return {
        "literature": lit, 
        "reactions": rxn,
        "visited_agents": ["Literature Retrieval", "Reaction Retrieval"]
    }

async def node_planning(state: AgentState):
    print("--- Retrosynthesis Planning ---")
    target = state["target_molecule"]
    molecule_name = state["molecule_data"].get("name", "")
    # Call plan() with target and molecule_name
    routes = retro_agent.plan(target, molecule_name)
    return {
        "routes": routes,
        "visited_agents": ["Retrosynthesis Planner"]
    }

async def node_analysis(state: AgentState):
    print("--- Route Analysis (Pred, Safety, Cost) ---")
    routes = state["routes"]
    molecule_name = state["molecule_data"].get("name", "")
    
    async def process_route(route):
        async def run_prediction():
            for step in route.get("steps", []):
                step["prediction"] = prediction_agent.predict(step)
        
        async def run_safety():
            route["safety"] = safety_agent.check(route)
        
        async def run_procurement():
            route["cost"] = procurement_agent.estimate(route)

        await asyncio.gather(run_prediction(), run_safety(), run_procurement())
        return route

    if routes:
        routes = await asyncio.gather(*(process_route(route) for route in routes))
    
    return {"routes": routes}

async def node_ranking(state: AgentState):
    print("--- Route Ranking ---")
    routes = state["routes"]
    ranked = scorer_agent.score(routes)
    return {
        "routes": ranked,
        "visited_agents": ["Route Scorer"]
    }

async def node_finalization(state: AgentState):
    print("--- Finalization (Protocol & Curation) ---")
    top_route = state["routes"][0] if state["routes"] else None
    molecule_name = state["molecule_data"].get("name", "")
    
    if top_route:
        protocol = protocol_agent.generate(top_route, molecule_name)
    else:
        protocol = f"No viable route found for {molecule_name}. Unable to generate protocol."
    
    # Log to curation
    curation_agent.log_provenance(
        event_type="synthesis_workflow",
        user_id="system",
        data={
            "molecule": molecule_name,
            "routes_count": len(state["routes"]),
            "protocol_generated": bool(protocol)
        }
    )
    
    return {
        "protocol": protocol,
        "visited_agents": ["Protocol Generator", "Data Curation"]
    }

# Define Standalone Nodes
async def node_safety_check(state: AgentState):
    print("--- Standalone Safety Check ---")
    target = state["target_molecule"]
    if not target:
        return {"errors": ["No target molecule"]}
    
    result = safety_agent.check_molecule(target)
    # Wrap in a mock route structure for frontend compatibility or just return raw
    return {
        "routes": [{"id": "safety_info", "safety": result}],
        "visited_agents": ["Safety Checker"]
    }

async def node_procurement_check(state: AgentState):
    print("--- Standalone Procurement Check ---")
    target = state["target_molecule"]
    if not target:
        return {"errors": ["No target molecule"]}
    
    result = procurement_agent.check_price(target)
    return {
        "routes": [{"id": "price_info", "cost": result}],
        "visited_agents": ["Procurement Agent"]
    }

# Define Conditional Logic
def router(state: AgentState):
    # Stop if there are errors (e.g., no target identified)
    if state.get("errors"):
        print("--- Router: Errors detected, stopping ---")
        return "end"

    intent_type = state["intent"].get("intent_type", "synthesis")
    print(f"--- Router: {intent_type} ---")
    
    if intent_type == "synthesis":
        return "synthesis"
    elif intent_type == "safety_check":
        return "safety"
    elif intent_type == "procurement":
        return "procurement"
    else:
        return "synthesis" # Default

# Build Graph
workflow = StateGraph(AgentState)

# Add Nodes
workflow.add_node("intent_parser", node_intent_parser)
workflow.add_node("canonicalizer", node_canonicalizer)

# Synthesis Pipeline Nodes
workflow.add_node("retrieval", node_retrieval)
workflow.add_node("planning", node_planning)
workflow.add_node("analysis", node_analysis)
workflow.add_node("ranking", node_ranking)
workflow.add_node("finalization", node_finalization)  # Updated from node_protocol

# Standalone Nodes
workflow.add_node("safety_check", node_safety_check)
workflow.add_node("procurement_check", node_procurement_check)

# Add Edges
workflow.set_entry_point("intent_parser")
workflow.add_edge("intent_parser", "canonicalizer")

# Router after Canonicalizer (so we have the clean name)
workflow.add_conditional_edges(
    "canonicalizer",
    router,
    {
        "synthesis": "retrieval",
        "safety": "safety_check",
        "procurement": "procurement_check",
        "end": END
    }
)

# Synthesis Flow
workflow.add_edge("retrieval", "planning")
workflow.add_edge("planning", "analysis")
workflow.add_edge("analysis", "ranking")
workflow.add_edge("ranking", "finalization")
workflow.add_edge("finalization", END)

# Standalone Flows
workflow.add_edge("safety_check", END)
workflow.add_edge("procurement_check", END)

# Compile
app_graph = workflow.compile()
