from typing import Dict, Any

class ProtocolGeneratorAgent:
    """
    Agent 9: Protocol Generator
    Generates lab-ready procedures.
    """
    def __init__(self):
        pass

    def generate(self, route: Dict[str, Any]) -> str:
        """
        Converts a route into a step-by-step text protocol.
        """
        protocol = f"### Protocol for Route {route.get('id', 'Unknown')}\n\n"
        
        for i, step in enumerate(route.get("steps", []), 1):
            reaction = step.get("reaction", "Unknown Reaction")
            conditions = step.get("conditions", "Standard Conditions")
            reagents = ", ".join(step.get("reagents", []))
            
            protocol += f"**Step {i}:**\n"
            protocol += f"- **Reaction:** {reaction}\n"
            protocol += f"- **Reagents:** {reagents}\n"
            protocol += f"- **Conditions:** {conditions}\n"
            protocol += f"- **Instructions:** Mix {reagents} under {conditions}. Monitor reaction progress. Quench and isolate product.\n\n"
            
        protocol += "**Safety Note:** Always wear appropriate PPE. Refer to SDS for all chemicals used."
        return protocol
