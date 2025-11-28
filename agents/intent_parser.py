import json
import re
from typing import Dict, Any, Optional

class IntentParserAgent:
    """
    Agent 1: User Intent / Parser
    Extracts user goal, constraints, delivery format, and priorities.
    """
    def __init__(self):
        # In a real scenario, we would load a model here.
        # self.nlp = spacy.load("en_core_web_sm")
        pass

    def parse(self, user_text: str, image_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Parses the user text and optional image to extract structured intent.
        """
        intent_data = {
            "target": None,
            "constraints": [],
            "mode": "automation", # default
            "priority": "normal",
            "original_text": user_text
        }

        # 1. Extract Target
        # Patterns to look for: "synthesis of X", "make X", "prepare X"
        target_patterns = [
            r"synthesis of\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|$)",
            r"make\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|$)",
            r"prepare\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|$)"
        ]
        
        for pattern in target_patterns:
            match = re.search(pattern, user_text, re.IGNORECASE)
            if match:
                target = match.group(1).strip()
                # Clean up trailing words that might have been caught
                stop_words = ["cheaply", "safely", "quickly", "it"]
                for word in stop_words:
                    if target.lower().endswith(f" {word}"):
                        target = target[:-len(word)-1]
                
                intent_data["target"] = target
                break
        else:
            # Fallback: assume the whole text might be the target if short, or look for keywords
            pass

        # 2. Extract Constraints (cost, safety, time)
        if "cheap" in user_text.lower() or "low cost" in user_text.lower():
            intent_data["constraints"].append("low_cost")
        if "safe" in user_text.lower() or "non-toxic" in user_text.lower():
            intent_data["constraints"].append("high_safety")
        if "fast" in user_text.lower() or "quick" in user_text.lower():
            intent_data["constraints"].append("fast_delivery")

        # 3. Extract Mode (bench vs automation)
        if "bench" in user_text.lower() or "manual" in user_text.lower():
            intent_data["mode"] = "bench"
        elif "automation" in user_text.lower() or "robot" in user_text.lower():
            intent_data["mode"] = "automation"

        # 4. Extract Priority
        if "urgent" in user_text.lower() or "asap" in user_text.lower():
            intent_data["priority"] = "high"

        # 5. Handle Image (Placeholder)
        if image_path:
            # In a real app, call OCR or Image2SMILES here
            intent_data["image_processed"] = True
            intent_data["image_path"] = image_path

        return intent_data

if __name__ == "__main__":
    agent = IntentParserAgent()
    sample_text = "I need the synthesis of Aspirin. It should be cheap and safe. I want to do this on the bench."
    print(json.dumps(agent.parse(sample_text), indent=2))
