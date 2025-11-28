from fastapi import FastAPI, HTTPException, File, UploadFile
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import uvicorn
import os
import sys
import json
import asyncio
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

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
from agents.chat_agent import ChatAgent

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
chat_agent = ChatAgent()

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
    visited_agents: List[str] = []
    error: Optional[str] = None

from graph import app_graph

from fastapi import File, UploadFile, Form
import shutil

@app.post("/api/stream")
async def stream_process(
    query: str = Form(...),
    file: Optional[UploadFile] = File(None)
):
    """
    Streams real-time agent execution updates using Server-Sent Events (SSE).
    """
    async def event_generator():
        temp_file_path = None
        try:
            # Handle file upload
            if file:
                temp_file_path = f"temp_{file.filename}"
                with open(temp_file_path, "wb") as buffer:
                    shutil.copyfileobj(file.file, buffer)
           
            # Initial state
            initial_state = {
                "user_query": query,
                "image_path": temp_file_path,
                "intent": {},
                "target_molecule": None,
                "molecule_data": {},
                "literature": [],
                "reactions": [],
                "routes": [],
                "protocol": None,
                "errors": [],
                "visited_agents": []
            }
            
            # Emit start event
            yield f"data: {json.dumps({'type': 'start', 'message': 'Starting agent pipeline'})}\n\n"
            await asyncio.sleep(0.1)
            
            # Execute graph with streaming
            agent_names = [
                "Intent Parser", "Canonicalizer", "Literature Retrieval",
                "Reaction Retrieval", "Retrosynthesis", "Forward Prediction",
                "Safety Checker", "Procurement", "Route Scorer",
                "Protocol Generator", "Data Curation"
            ]
            
            current_agent_idx = 0
            
            # Stream each step
            async for chunk in app_graph.astream(initial_state):
                # Emit agent progress
                if current_agent_idx < len(agent_names):
                    agent_data = {
                        'type': 'agent',
                        'name': agent_names[current_agent_idx],
                        'status': 'executing',
                        'index': current_agent_idx
                    }
                    yield f"data: {json.dumps(agent_data)}\n\n"
                    await asyncio.sleep(0.1)
                    
                    # Mark as complete
                    agent_data['status'] = 'complete'
                    yield f"data: {json.dumps(agent_data)}\n\n"
                    
                    current_agent_idx += 1
            
            # Get final state
            final_state = await app_graph.ainvoke(initial_state)
            
            # Clean up temp file
            if temp_file_path and os.path.exists(temp_file_path):
                os.remove(temp_file_path)
            
            # Emit final result
            result_data = {
                'type': 'result',
                'data': {
                    'intent': final_state.get("intent", {}),
                    'canonical_data': final_state.get("molecule_data", {}),
                    'literature_results': final_state.get("literature", []),
                    'reaction_results': final_state.get("reactions", []),
                    'routes': final_state.get("routes", []),
                    'protocol': final_state.get("protocol"),
                    'visited_agents': final_state.get("visited_agents", [])
                }
            }
            yield f"data: {json.dumps(result_data)}\n\n"
            
            # Emit done event
            yield f"data: {json.dumps({'type': 'done'})}\n\n"
            
        except Exception as e:
            print(f"Streaming error: {e}")
            error_data = {
                'type': 'error',
                'message': str(e)
            }
            yield f"data: {json.dumps(error_data)}\n\n"
            
            # Clean up temp file on error
            if temp_file_path and os.path.exists(temp_file_path):
                os.remove(temp_file_path)
    
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
            "Connection": "keep-alive"
        }
    )

@app.post("/api/process", response_model=AgentResponse)
async def process_query(
    query: str = Form(...),
    file: Optional[UploadFile] = File(None)
):
    """
    Orchestrates the full agent workflow using LangGraph.
    Accepts multipart/form-data to handle optional image uploads.
    """
    temp_file_path = None
    try:
        print(f"Processing query via Graph: {query}")
        
        # Handle File Upload
        if file:
            temp_file_path = f"temp_{file.filename}"
            with open(temp_file_path, "wb") as buffer:
                shutil.copyfileobj(file.file, buffer)
            print(f"Saved uploaded file to {temp_file_path}")

        # Initial State
        initial_state = {
            "user_query": query,
            "image_path": temp_file_path,
            "intent": {},
            "target_molecule": None,
            "molecule_data": {},
            "literature": [],
            "reactions": [],
            "routes": [],
            "protocol": None,
            "errors": [],
            "visited_agents": []
        }

        # Invoke Graph
        final_state = await app_graph.ainvoke(initial_state)
        
        # Clean up temp file
        if temp_file_path and os.path.exists(temp_file_path):
            os.remove(temp_file_path)
            print(f"Removed temp file {temp_file_path}")
        
        # Extract Results
        return AgentResponse(
            intent=final_state.get("intent", {}),
            molecule_info=final_state.get("molecule_data", {}),
            literature=final_state.get("literature", []),
            reactions=final_state.get("reactions", []),
            routes=final_state.get("routes", []),
            protocol=final_state.get("protocol"),
            visited_agents=final_state.get("visited_agents", []),
            error=str(final_state.get("errors")) if final_state.get("errors") else None
        )
        
    except Exception as e:
        print(f"Error: {e}")
        # Clean up temp file on error
        if temp_file_path and os.path.exists(temp_file_path):
            os.remove(temp_file_path)
            
        return AgentResponse(
            intent={},
            molecule_info={},
            literature=[],
            reactions=[],
            error=str(e),
            visited_agents=[]
        )

@app.get("/health")
async def health_check():
    return {"status": "ok"}

@app.post("/api/feedback")
async def submit_feedback(route_id: str, rating: int, comments: str):
    return feedback_agent.process_feedback(route_id, rating, comments)

class ChatMessage(BaseModel):
    message: str
    history: Optional[List[Dict[str, str]]] = []

@app.post("/api/chat")
async def chat(chat_message: ChatMessage):
    """
    Conversational AI endpoint for natural language chemistry questions.
    """
    try:
        response = chat_agent.chat(chat_message.message, chat_message.history)
        return response
    except Exception as e:
        print(f"Chat error: {e}")
        return {
            "response": "I apologize, but I encountered an error. Please try rephrasing your question.",
            "history": chat_message.history,
            "context": {},
            "suggestions": ["Ask about synthesis", "Ask about safety", "Ask about mechanisms"]
        }

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
