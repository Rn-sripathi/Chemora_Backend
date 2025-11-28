"""
Conversational AI Chat Agent for Chemora
Handles natural language chemistry questions with context awareness
Integrates with all backend agents for professional-grade responses
"""

from typing import List, Dict, Any, Optional
import json
import os

# Try to import LLM libraries (OpenAI or Gemini)
try:
    import openai
    openai_available = True
except ImportError:
    openai_available = False

try:
    import google.generativeai as genai
    genai_available = True
except ImportError:
    genai_available = False

# Check which LLM is configured
llm_type = None
if os.getenv("GEMINI_API_KEY"):
    llm_type = "gemini"
    if genai_available:
        genai.configure(api_key=os.getenv("GEMINI_API_KEY"))
elif os.getenv("OPENAI_API_KEY"):
    llm_type = "openai"
    if openai_available:
        openai.api_key = os.getenv("OPENAI_API_KEY")

print(f"Chat Agent LLM: {llm_type or 'None (agent-only mode)'}")

# Import all Chemora agents for professional responses
try:
    from agents.retrosynthesis import RetrosynthesisAgent
    from agents.literature_retrieval import LiteratureRetrievalAgent
    from agents.reaction_retrieval import ReactionRetrievalAgent
    from agents.safety_checker import SafetyCheckerAgent
    from agents.protocol_generator import ProtocolGeneratorAgent
    from agents.canonicalizer import CanonicalizerAgent
    agents_available = True
except ImportError:
    agents_available = False
    print("Warning: Chemora agents not available for chat integration")

class ChatAgent:
    """
    Conversational AI agent that:
    - Maintains conversation context
    - Routes questions to specialized agents
    - Provides natural language responses
    - Handles follow-up questions
    - Integrates with full Chemora agent pipeline
    """
    
    def __init__(self):
        self.conversation_history = []
        self.context = {}
        
        # Initialize all Chemora agents
        if agents_available:
            self.retro_agent = RetrosynthesisAgent()
            self.literature_agent = LiteratureRetrievalAgent()
            self.reaction_agent = ReactionRetrievalAgent()
            self.safety_agent = SafetyCheckerAgent()
            self.protocol_agent = ProtocolGeneratorAgent()
            self.canonicalizer = CanonicalizerAgent()
            print("Chat Agent: Initialized with full agent pipeline")
        else:
            self.retro_agent = None
            self.literature_agent = None
            self.reaction_agent = None
            self.safety_agent = None
            self.protocol_agent = None
            self.canonicalizer = None
            print("Chat Agent: Running in standalone mode (agents not available)")
        
        # System prompt for chemistry assistant
        self.system_prompt = """You are Chemora, an expert AI chemistry assistant - essentially "ChatGPT for Chemistry". You are:

**Core Capabilities:**
- An expert in organic, inorganic, physical, analytical, and computational chemistry
- Knowledgeable about synthesis, mechanisms, spectroscopy, biochemistry, and materials science
- Able to explain complex concepts clearly at any level (high school to PhD)
- Skilled at problem-solving, troubleshooting, and providing practical lab advice

**Special Powers (Agent Integration):
- You have access to specialized agents for:
  * Retrosynthetic analysis (ML-powered route planning)
  * Safety assessment (GHS codes, PPE, hazard analysis)
  * Literature search (precedent databases)
  * Protocol generation (with RDKit stoichiometry)
  * Molecule canonicalization

**Communication Style:**
- Conversational and friendly, like ChatGPT
- Scientifically rigorous and accurate
- Use proper chemical notation (H₂O, H₂SO₄, etc.)
- Explain reasoning step-by-step when helpful
- Ask clarifying questions when needed
- Admit uncertainty rather than making up information

**Response Format:**
- Use markdown formatting (headers, bold, lists, code blocks)
- Include relevant equations, mechanisms, or structures when helpful
- Provide safety warnings for hazardous materials
- Cite sources or databases when available
- Offer follow-up suggestions

**Examples of Questions You Can Answer:**
- "Explain why benzene is aromatic"
- "How do I run an NMR experiment?"
- "What's the mechanism of the aldol condensation?"
- "Help me troubleshoot my recrystallization"
- "Calculate the pH of a buffer solution"
- "What are the best solvents for nucleophilic substitution?"
- "Design a synthesis of ibuprofen"

You can discuss chemistry broadly - from basic concepts to advanced research topics. Be helpful, accurate, and comprehensive."""
    
    def chat(self, user_message: str, history: List[Dict[str, str]] = None) -> Dict[str, Any]:
        """
        Process a conversational message and return AI response.
        Uses LLM (ChatGPT-style) with agent augmentation for specialized tasks.
        
        Args:
            user_message: User's question/message
            history: Previous conversation history
            
        Returns:
            Dict with 'response', 'context', and 'suggestions'
        """
        if history:
            self.conversation_history = history
        
        # Add user message to history
        self.conversation_history.append({
            "role": "user",
            "content": user_message
        })
        
        # Detect if this requires specialized agent processing
        intent = self._classify_intent(user_message)
        should_use_agents = intent in ["synthesis", "safety", "literature"]
        
        # If it's a specialized chemistry task AND agents are available, use them
        if should_use_agents and agents_available:
            if intent == "synthesis":
                response = self._handle_synthesis_question(user_message)
            elif intent == "safety":
                response = self._handle_safety_question(user_message)
            elif intent == "literature":
                response = self._handle_literature_question(user_message)
            else:
                # Fallback to general LLM if intent is agent-related but not explicitly handled here
                response = self._generate_llm_response(user_message)
        else:
            # For all other questions, use LLM (ChatGPT-style)
            # This includes: mechanisms, concepts, troubleshooting, calculations, etc.
            response = self._generate_llm_response(user_message)
        
        # Add AI response to history
        self.conversation_history.append({
            "role": "assistant",
            "content": response
        })
        
        return {
            "response": response,
            "history": self.conversation_history,
            "context": self.context,
            "suggestions": self._generate_suggestions(user_message, response)
        }
    
    def _classify_intent(self, message: str) -> str:
        """Classify user intent from message"""
        message_lower = message.lower()
        
        synthesis_keywords = ["synthesis", "make", "prepare", "synthesize", "route", "procedure", "protocol"]
        safety_keywords = ["safe", "hazard", "toxic", "danger", "ppe", "risk"]
        mechanism_keywords = ["mechanism", "how does", "why", "reaction pathway"]
        literature_keywords = ["literature", "precedent", "paper", "reference", "research", "publication", "study"]
        
        if any(kw in message_lower for kw in synthesis_keywords):
            return "synthesis"
        elif any(kw in message_lower for kw in safety_keywords):
            return "safety"
        elif any(kw in message_lower for kw in mechanism_keywords):
            return "mechanism"
        elif any(kw in message_lower for kw in literature_keywords):
            return "literature"
        else:
            return "general"
    
    def _handle_synthesis_question(self, message: str) -> str:
        """Handle synthesis-related questions using actual retrosynthesis agent"""
        # Extract molecule name
        molecule = self._extract_molecule_name(message)
        
        if not molecule:
            return "I'd be happy to help with **synthesis planning**!\n\n" \
                   "Please tell me which molecule you're interested in synthesizing.\n" \
                   "For example: *'How do I synthesize aspirin?'* or *'I need to make paracetamol'*"
        
        # Store in context
        self.context["current_molecule"] = molecule
        
        # If agents not available, provide basic info
        if not agents_available or not self.retro_agent:
            return self._get_fallback_synthesis_info(molecule)
        
        try:
            # Step 1: Canonicalize the molecule
            canonical_result = self.canonicalizer.canonicalize(molecule)
            
            if canonical_result.get("error"):
                return f"I couldn't find structure information for **{molecule}**.\n\n" \
                       f"Could you provide:\n" \
                       f"- The IUPAC name\n" \
                       f"- A common name\n" \
                       f"- Or the SMILES notation\n\n" \
                       f"I can then provide a complete synthesis plan!"
            
            smiles = canonical_result.get("canonical_smiles", "")
            mol_name = canonical_result.get("name", molecule)
            formula = canonical_result.get("formula", "")
            mw = canonical_result.get("molecular_weight", "")
            
            # Store canonicalized data
            self.context["smiles"] = smiles
            self.context["formula"] = formula
            
            # Step 2: Get retrosynthesis routes
            routes = self.retro_agent.plan(smiles, mol_name)
            
            if not routes:
                return f"I analyzed **{mol_name}** but couldn't find viable synthesis routes.\n\n" \
                       f"This might be because:\n" \
                       f"- The molecule is very complex\n" \
                       f"- Limited precedent data available\n\n" \
                       f"Would you like me to search for literature precedents instead?"
            
            # Step 3: Get safety information
            safety = self.safety_agent.check_molecule(smiles)
            hazards = safety.get("hazards", [])
            
            # Step 4: Format response
            response = f"# Synthesis Planning for {mol_name}\n\n"
            response += f"**Molecular Formula:** {formula}\n"
            response += f"**Molecular Weight:** {mw} g/mol\n\n"
            
            # Show top routes
            response += f"## Retrosynthetic Routes\n\n"
            response += f"I found **{len(routes)} viable route(s)** for synthesizing {mol_name}:\n\n"
            
            for i, route in enumerate(routes[:3], 1):  # Show top 3
                response += f"### Route {i}"
                if route.get("confidence"):
                    response += f" (Confidence: {route['confidence']*100:.0f}%)"
                response += "\n\n"
                
                # Show steps
                for j, step in enumerate(route.get("steps", []), 1):
                    reaction = step.get("reaction", "Unknown reaction")
                    conditions = step.get("conditions", "")
                    response += f"{j}. **{reaction}**"
                    if conditions:
                        response += f"\n   - Conditions: {conditions}"
                    if step.get("prediction", {}).get("predicted_yield"):
                        response += f"\n   - Expected yield: {step['prediction']['predicted_yield']}%"
                    response += "\n\n"
            
            # Safety information
            if hazards:
                response += "## ⚠️ Safety Considerations\n\n"
                for hazard in hazards[:5]:
                    if isinstance(hazard, dict):
                        response += f"- {hazard.get('hazard', hazard.get('code', 'Unknown hazard'))}\n"
                    else:
                        response += f"- {hazard}\n"
                response += "\n"
            
            response += "*Would you like me to generate a detailed protocol for any specific route?*"
            
            return response
            
        except Exception as e:
            print(f"Synthesis planning error: {e}")
            return f"I encountered an issue analyzing **{molecule}**.\n\n" \
                   f"Let me try a different approach. Could you tell me:\n" \
                   f"- What starting materials do you have available?\n" \
                   f"- What scale are you working at?\n" \
                   f"- Any specific constraints (reagents to avoid, equipment limitations)?"
    
    def _get_fallback_synthesis_info(self, molecule: str) -> str:
        """Fallback synthesis info when agents not available"""
        molecule_lower = molecule.lower()
        
        # Common synthesis information (hardcoded for fallback)
        synthesis_info = {
            "aspirin": {
                "name": "Aspirin (Acetylsalicylic Acid)",
                "routes": [
                    "**Industrial Route:** Acetylation of salicylic acid with acetic anhydride in presence of acid catalyst (H₂SO₄ or H₃PO₄)",
                    "**Laboratory Route:** Salicylic acid + Acetic anhydride → Aspirin + Acetic acid (85-90% yield)"
                ],
                "conditions": "Temperature: 50-60°C, Time: 30-45 minutes",
                "safety": "Use fume hood, wear gloves. Acetic anhydride is corrosive and lachrymatory.",
                "reagents": "- Salicylic acid (2.0 g)\n- Acetic anhydride (5 mL)\n- Concentrated H₂SO₄ (3-5 drops as catalyst)\n- Ice bath for crystallization"
            },
            "paracetamol": {
                "name": "Paracetamol (Acetaminophen)",
                "routes": [
                    "**Industrial Route:** Acetylation of p-aminophenol with acetic anhydride",
                    "**Laboratory Route:** p-Aminophenol + Acetic anhydride → Paracetamol (75-85% yield)"
                ],
                "conditions": "Temperature: Room temperature to 60°C, Time: 30 minutes",
                "safety": "Use fume hood, handle acetic anhydride with care.",
                "reagents": "- p-Aminophenol (1.5 g)\n- Acetic anhydride (2 mL)\n- Sodium acetate (catalyst)\n- Water for recrystallization"
            },
            "ibuprofen": {
                "name": "Ibuprofen",
                "routes": [
                    "**Boots Process:** 6-step synthesis from isobutylbenzene via Friedel-Crafts acylation",
                    "**BHC (Hoechst) Process:** 3-step catalytic synthesis (greener alternative, 77% overall yield)"
                ],
                "conditions": "Multi-step process, requires specialized equipment",
                "safety": "Anhydrous conditions, inert atmosphere required for some steps.",
                "reagents": "Starting material: Isobutylbenzene or alternatives depending on route"
            }
        }
        
        for known_mol, info in synthesis_info.items():
            if known_mol in molecule_lower or molecule_lower in known_mol:
                return f"# Synthesis of {info['name']}\n\n" \
                       f"**Common Routes:**\n" + "\n".join(f"{i+1}. {route}" for i, route in enumerate(info['routes'])) + \
                       f"\n\n**Typical Conditions:**\n{info['conditions']}\n\n" \
                       f"**Safety Precautions:**\n{info['safety']}\n\n" \
                       f"**Reagents Needed:**\n{info['reagents']}\n\n" \
                       f"*Would you like a detailed step-by-step protocol for any specific route?*"
        
        return f"I can help you with the synthesis of **{molecule}**!\n\n" \
               f"For the best synthesis plan, I recommend:\n\n" \
               f"1. **Retrosynthetic Analysis:** Identify key functional groups and disconnections\n" \
               f"2. **Literature Search:** Check SciFinder, Reaxys, or organic synthesis databases\n" \
               f"3. **Route Selection:** Consider yield, cost, safety, and feasibility\n\n" \
               f"Would you like me to:\n" \
               f"- Search for literature precedents?\n" \
               f"- Suggest common reagents for similar transformations?\n" \
               f"- Provide safety information?"
    
    
    def _handle_safety_question(self, message: str) -> str:
        """Handle safety-related questions using actual safety checker"""
        # Check if we have a molecule in context
        if "smiles" in self.context and agents_available and self.safety_agent:
            smiles = self.context["smiles"]
            molecule = self.context.get("current_molecule", "this molecule")
            
            try:
                safety = self.safety_agent.check_molecule(smiles)
                
                response = f"# Safety Assessment for {molecule}\n\n"
                
                # Hazards
                hazards = safety.get("hazards", [])
                if hazards:
                    response += "## ⚠️ Identified Hazards\n\n"
                    for hazard in hazards:
                        if isinstance(hazard, dict):
                            response += f"- **{hazard.get('code', 'Unknown')}:** {hazard.get('description', hazard.get('hazard', ''))}\n"
                        else:
                            response += f"- {hazard}\n"
                    response += "\n"
                
                # PPE recommendations
                if safety.get("ppe_required"):
                    response += "## Required PPE\n\n"
                    for ppe in safety["ppe_required"]:
                        response += f"- {ppe}\n"
                    response += "\n"
                
                # Handling recommendations
                if safety.get("handling"):
                    response += "## Handling Recommendations\n\n"
                    response += f"{safety['handling']}\n\n"
                
                # Risk level
                if safety.get("risk_level"):
                    response += f"**Overall Risk Level:** {safety['risk_level']}\n\n"
                
                response += "*Always consult the Safety Data Sheet (SDS) before handling any chemical.*"
                
                return response
                
            except Exception as e:
                print(f"Safety check error: {e}")
        
        # Fallback response
        return "**General Safety Guidelines:**\n\n" \
               "For any chemical synthesis:\n" \
               "- **Always** consult Safety Data Sheets (SDS)\n" \
               "- Wear appropriate PPE: safety goggles, lab coat, nitrile gloves\n" \
               "- Work in a well-ventilated fume hood\n" \
               "- Have emergency equipment ready (eyewash, shower, fire extinguisher)\n" \
               "- Know the location of spill kits and first aid supplies\n" \
               "- Never work alone in the lab\n" \
               "- Properly label all containers\n\n" \
               "*What specific molecule or reaction do you need safety information for?*"
    
    def _handle_literature_question(self, message: str) -> str:
        """Handle literature search questions using actual literature retrieval agent"""
        molecule = self.context.get("current_molecule") or self._extract_molecule_name(message)
        
        if not molecule:
            return "I can help you search for **literature and precedents**!\n\n" \
                   "Please specify:\n" \
                   "- Which molecule or reaction you're researching\n" \
                   "- What type of information you need (synthesis procedures, characterization, properties)\n\n" \
                   "For example: *'Find literature on aspirin synthesis'*"
        
        if agents_available and self.literature_agent:
            try:
                # Search literature
                results = self.literature_agent.retrieve(f"Synthesis of {molecule}", molecule)
                
                if results:
                    response = f"# Literature Results for {molecule}\n\n"
                    response += f"I found **{len(results)} relevant reference(s)**:\n\n"
                    
                    for i, doc in enumerate(results[:5], 1):
                        response += f"### [{i}] {doc.get('text', 'Document')}\n\n"
                        if doc.get('source'):
                            response += f"- **Source:** {doc['source']}\n"
                        if doc.get('yield'):
                            response += f"- **Yield:** {doc['yield']}\n"
                        if doc.get('score'):
                            response += f"- **Relevance Score:** {doc['score']:.2f}\n"
                        response += "\n"
                    
                    response += "*Would you like me to provide synthesis protocols based on these precedents?*"
                    return response
                else:
                    return f"I searched for literature on **{molecule}** but couldn't find specific precedents in my database.\n\n" \
                           f"I recommend:\n" \
                           f"- Checking SciFinder or Reaxys for comprehensive coverage\n" \
                           f"- Searching PubMed for biological/medicinal chemistry aspects\n" \
                           f"- Looking at Organic Syntheses for verified procedures\n\n" \
                           f"Would you like me to suggest a synthesis route instead?"
            
            except Exception as e:
                print(f"Literature search error: {e}")
        
        return f"I can help search for literature on **{molecule}**!\n\n" \
               f"Unfortunately, my literature database is currently limited. For comprehensive searches, I recommend:\n" \
               f"- **SciFinder** - Most comprehensive chemistry database\n" \
               f"- **Reaxys** - Excellent for synthesis and reaction databases\n" \
               f"- **PubMed** - For biological and medicinal chemistry\n" \
               f"- **Organic Syntheses** - Verified procedures with detailed protocols\n\n" \
               f"Would you like me to suggest a synthesis route for {molecule} instead?"
    
    def _handle_mechanism_question(self, message: str) -> str:
        """Handle reaction mechanism questions"""
        return "I can explain reaction mechanisms! For the best explanation, please specify:\n" \
               "- The reaction type (e.g., nucleophilic substitution, elimination)\n" \
               "- The specific transformation\n" \
               "- Starting materials and products"
    
    def _handle_general_question(self, message: str) -> str:
        """Handle general chemistry questions"""
        # Check if this is a follow-up (yes/no responses)
        message_lower = message.lower().strip()
        
        # Handle affirmative responses
        if message_lower in ['yes', 'yes please', 'sure', 'ok', 'okay', 'yeah', 'yep', 'please']:
            if "current_molecule" in self.context:
                molecule = self.context["current_molecule"]
                return self._get_detailed_protocol(molecule)
            else:
                return "Great! What would you like to know about? I can help with:\n" \
                       "- **Synthesis planning** for specific molecules\n" \
                       "- **Safety information** for chemicals and reactions\n" \
                       "- **Reaction mechanisms** and explanations\n" \
                       "- **Literature searches** for precedents"
        
        # Handle negative responses
        if message_lower in ['no', 'no thanks', 'nope', 'not now']:
            return "No problem! Is there anything else I can help you with regarding chemistry?"
        
        # Try LLM if available, otherwise provide helpful response
        if openai_available and openai.api_key:
            return self._generate_llm_response(message)
        else:
            return "I'm here to help with chemistry questions! I specialize in:\n\n" \
                   "- **Synthesis Planning:** Route design, reagent selection, reaction conditions\n" \
                   "- **Safety Assessment:** Hazard identification, PPE recommendations, handling procedures\n" \
                   "- **Mechanism Explanations:** Reaction pathways, electron movement, intermediate identification\n" \
                   "- **Literature Searches:** Finding precedents and procedures\n\n" \
                   "What can I help you with today?"
    
    def _get_detailed_protocol(self, molecule: str) -> str:
        """Generate a detailed protocol using actual protocol generator agent"""
        
        # If agents available and we have SMILES in context
        if agents_available and self.protocol_agent and "smiles" in self.context:
            smiles = self.context["smiles"]
            
            try:
                # Get retrosynthesis routes first
                routes = self.retro_agent.plan(smiles, molecule) if self.retro_agent else []
                
                if routes:
                    # Use the top-ranked route
                    best_route = routes[0]
                    
                    # Generate protocol
                    protocol = self.protocol_agent.generate(
                        route=best_route,
                        target_smiles=smiles,
                        molecule_name=molecule
                    )
                    
                    if protocol:
                        return f"# Detailed Synthesis Protocol for {molecule}\n\n{protocol}\n\n" \
                               "*This protocol was generated using RDKit stoichiometry calculations, " \
                               "Jinja2 templates, and retrosynthetic analysis.*\n\n" \
                               "*Always verify with literature and perform risk assessment before starting!*"
                
            except Exception as e:
                print(f"Protocol generation error: {e}")
        
        # Fallback to hardcoded protocols for common molecules
        molecule_lower = molecule.lower()
        
        if "aspirin" in molecule_lower:
            return """# Detailed Aspirin Synthesis Protocol

**Objective:** Synthesize aspirin (acetylsalicylic acid) via acetylation of salicylic acid

**Materials:**
- Salicylic acid: 2.0 g (14.5 mmol)
- Acetic anhydride: 5.0 mL (53 mmol, 3.5 equiv)
- Concentrated H₂SO₄: 3-5 drops (catalyst)
- Ice bath
- 50 mL Erlenmeyer flask
- Buchner funnel and filter paper
- Vacuum filtration setup

**Procedure:**

**Step 1: Setup**
1. Set up a water bath at 50-60°C
2. Weigh 2.0 g salicylic acid into a dry 50 mL Erlenmeyer flask
3. Work in a fume hood for all steps involving acetic anhydride

**Step 2: Acetylation Reaction**
1. Add 5 mL acetic anhydride to the salicylic acid
2. Add 3-5 drops of concentrated H₂SO₄ as catalyst
3. Swirl gently to mix (exothermic - be careful!)
4. Place flask in 50-60°C water bath for 15 minutes
5. Swirl occasionally to ensure complete dissolution

**Step 3: Crystallization**
1. Remove flask from water bath
2. Carefully add 20 mL of cold water to decompose excess acetic anhydride
3. Cool in ice bath for 15-20 minutes
4. Scratch sides of flask to induce crystallization if needed

**Step 4: Isolation**
1. Collect crystals by vacuum filtration
2. Wash with cold water (2 × 10 mL)
3. Allow to air dry or use vacuum desiccator

**Step 5: Purification (Optional)**
1. Recrystallize from ethanol/water mixture
2. Expected yield: 1.7-2.0 g (85-95%)
3. Melting point should be 138-140°C

**Safety Notes:**
- ⚠️ Acetic anhydride is corrosive and causes severe burns
- ⚠️ H₂SO₄ is highly corrosive - use extreme caution
- ⚠️ Work in fume hood - acetic anhydride vapors irritate eyes and respiratory system
- Wear safety goggles, lab coat, and nitrile gloves

**Waste Disposal:**
- Neutralize acidic filtrate with sodium bicarbonate before disposal
- Dispose in appropriate aqueous waste container

*Would you like me to explain the reaction mechanism or discuss alternatives?*"""
        
        elif "paracetamol" in molecule_lower or "acetaminophen" in molecule_lower:
            return f"# Detailed Paracetamol Synthesis Protocol\n\n" \
                   f"**Starting Material:** p-Aminophenol\n" \
                   f"**Reaction:** Acetylation with acetic anhydride\n\n" \
                   f"**Procedure:**\n" \
                   f"1. Dissolve 1.5 g p-aminophenol in 25 mL water\n" \
                   f"2. Add 2 mL acetic anhydride with stirring\n" \
                   f"3. Add sodium acetate as catalyst\n" \
                   f"4. Stir at room temperature for 30 minutes\n" \
                   f"5. Cool in ice bath to crystallize\n" \
                   f"6. Filter and wash with cold water\n" \
                   f"7. Recrystallize from water/ethanol\n\n" \
                   f"**Expected Yield:** 75-85%\n" \
                   f"**Melting Point:** 169-171°C"
        
        else:
            return f"I can provide a detailed protocol for **{molecule}**!\n\n" \
                   f"However, I need more specific information:\n" \
                   f"- What is the starting material you have available?\n" \
                   f"- What scale are you working at (mg, g, kg)?\n" \
                   f"- Do you have any specific constraints (reagents to avoid, equipment limitations)?\n\n" \
                   f"Alternatively, try the **Synthesis Planner** tab for a comprehensive route analysis with automated protocol generation!"
    
    def _generate_llm_response(self, message: str) -> str:
        """Generate response using LLM (Gemini or GPT-4) - the core ChatGPT-style interface"""
        
        # Check if any LLM is available
        if llm_type == "gemini":
            return self._generate_gemini_response(message)
        elif llm_type == "openai":
            return self._generate_openai_response(message)
        else:
            return self._get_no_llm_fallback(message)
    
    def _generate_gemini_response(self, message: str) -> str:
        """Generate response using Google Gemini API"""
        if not genai_available:
            return "⚠️ Gemini library not installed. Run: pip install google-generativeai\n\n" + \
                   self._get_no_llm_fallback(message)
        
        try:
            # Initialize Gemini model
            # Fallback to the most stable model: gemini-1.0-pro
            try:
                model = genai.GenerativeModel(
                    model_name='gemini-1.0-pro',
                    system_instruction=self.system_prompt
                )
            except Exception as e:
                print(f"Gemini 1.0 Pro failed: {e}")
                print("Available models:")
                for m in genai.list_models():
                    print(f"- {m.name}")
                raise e
            
            # Build conversation history for Gemini
            chat_history = []
            for msg in self.conversation_history[-10:]:  # Last 10 for context
                role = "user" if msg["role"] == "user" else "model"
                chat_history.append({
                    "role": role,
                    "parts": [msg["content"]]
                })
            
            # Add context if available
            user_message = message
            if self.context.get("current_molecule"):
                user_message += f"\n\n[Context: Discussing {self.context['current_molecule']}"
                if self.context.get("smiles"):
                    user_message += f", SMILES: {self.context['smiles']}"
                user_message += "]"
            
            # Start chat session
            chat = model.start_chat(history=chat_history[:-1] if chat_history else [])
            
            # Generate response
            response = chat.send_message(user_message)
            
            return response.text
            
        except Exception as e:
            print(f"Gemini Error: {e}")
            return f"⚠️ Gemini API error: {str(e)}\n\n" + \
                   "Please check your GEMINI_API_KEY is set correctly.\n\n" + \
                   self._get_no_llm_fallback(message)
    
    def _generate_openai_response(self, message: str) -> str:
        """Generate response using OpenAI GPT-4"""
        if not openai_available:
            return self._get_no_llm_fallback(message)
        
        # Check if API key is set
        if not openai.api_key:
            return "I'd love to help, but I need an OpenAI API key to answer general chemistry questions.\n\n" \
                   "In the meantime, I can still help with:\n" \
                   "- Synthesis planning (using our agents)\n" \
                   "- Safety assessments\n" \
                   "- Literature searches\n\n" \
                   "What would you like to explore?"
        
        try:
            # Build context from conversation history
            messages = [{"role": "system", "content": self.system_prompt}]
            
            # Add recent conversation (last 10 messages for context)
            messages.extend(self.conversation_history[-10:])
            
            # Add agent context if available
            if self.context.get("current_molecule"):
                context_note = f"\n\n[Context: User is asking about {self.context['current_molecule']}"
                if self.context.get("smiles"):
                    context_note += f", SMILES: {self.context['smiles']}"
                context_note += "]"
                messages[-1]["content"] += context_note
            
            # Call GPT-4
            response = openai.ChatCompletion.create(
                model="gpt-4o-mini",  # Can upgrade to gpt-4 for better responses
                messages=messages,
                temperature=0.7,
                max_tokens=800,  # Increased for more detailed explanations
                presence_penalty=0.1,
                frequency_penalty=0.1
            )
            
            return response.choices[0].message.content
            
        except openai.error.AuthenticationError:
            return "⚠️ OpenAI API authentication failed. Please check your API key.\n\n" \
                   "I can still help with synthesis planning, safety, and literature using our specialized agents!"
        
        except openai.error.RateLimitError:
            return "⚠️ OpenAI API rate limit reached. Please try again in a moment.\n\n" \
                   "In the meantime, I can answer synthesis-related questions using our agents."
        
        except openai.error.InvalidRequestError as e:
            print(f"LLM Invalid Request: {e}")
            return "I had trouble processing that question. Could you rephrase it?\n\n" \
                   "For synthesis planning, safety, or literature questions, I can use our specialized agents instead!"
        
        except Exception as e:
            print(f"OpenAI Error: {e}")
            return self._get_no_llm_fallback(message)
    
    def _get_no_llm_fallback(self, message: str) -> str:
        """Fallback response when LLM is unavailable"""
        return "I'm a chemistry AI assistant, but I'm currently running without my general knowledge module.\n\n" \
               "**I can still help you with:**\n" \
               "- **Synthesis Planning:** Ask about making specific molecules\n" \
               "- **Safety Information:** Get hazard assessments and PPE recommendations\n" \
               "- **Literature Search:** Find precedents and references\n" \
               "- **Protocol Generation:** Get detailed lab procedures\n\n" \
               "**Examples:**\n" \
               "- 'How do I synthesize aspirin?'\n" \
               "- 'What are the safety concerns for working with acetic anhydride?'\n" \
               "- 'Find literature on Suzuki coupling reactions'\n\n" \
               "What would you like to know?"
    
    def _extract_molecule_name(self, text: str) -> Optional[str]:
        """Extract molecule name from user message"""
        # Simple extraction - look for common patterns
        import re
        
        # Pattern: "synthesis of [molecule]"
        match = re.search(r'(?:synthesis|make|prepare|synthesize)\s+(?:of\s+)?([A-Za-z0-9\-]+)', text, re.IGNORECASE)
        if match:
            return match.group(1).strip()
        
        # Pattern: "how to make [molecule]"
        match = re.search(r'how\s+to\s+(?:make|prepare|synthesize)\s+([A-Za-z0-9\-]+)', text, re.IGNORECASE)
        if match:
            return match.group(1).strip()
        
        return None
    
    def _generate_suggestions(self, user_message: str, ai_response: str) -> List[str]:
        """Generate follow-up suggestions"""
        intent = self._classify_intent(user_message)
        
        suggestions = {
            "synthesis": [
                "Show me a detailed protocol",
                "What are the safety concerns?",
                "Are there alternative routes?",
                "Calculate stoichiometry for 10g scale",
                "Find literature precedents"
            ],
            "safety": [
                "What PPE should I use?",
                "How do I dispose of waste?",
                "What are emergency procedures?",
                "Is there a safer alternative?"
            ],
            "mechanism": [
                "Show the electron-pushing mechanism",
                "What are the key intermediates?",
                "Why is this regioselective?",
                "How does temperature affect this?"
            ],
            "literature": [
                "Summarize the key findings",
                "Are there more recent studies?",
                "What's the best procedure?",
                "Compare yields from different methods"
            ],
            "general": [
                "Explain this concept in simple terms",
                "How do I synthesize aspirin?",
                "What's the difference between SN1 and SN2?",
                "Help me design an experiment",
                "Troubleshoot my reaction"
            ]
        }
        
        return suggestions.get(intent, suggestions["general"])
    
    def reset_conversation(self):
        """Clear conversation history and context"""
        self.conversation_history = []
        self.context = {}
