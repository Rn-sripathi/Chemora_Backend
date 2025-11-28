import json
import re
import uuid
import time
import os
from typing import Dict, Any, Optional

try:
    import spacy
    from transformers import pipeline
    import pytesseract
    from PIL import Image
    # Optional: DECIMER for Image2SMILES
    try:
        from DECIMER import predict_SMILES
    except ImportError:
        predict_SMILES = None
except ImportError:
    spacy = None
    pipeline = None
    pytesseract = None
    Image = None
    predict_SMILES = None

class IntentParserAgent:
    """
    Agent 1: User Intent / Parser (Advanced)
    Extracts user goal, constraints, delivery format, and priorities using NLP and OCR.
    """
    def __init__(self):
        self.request_db = []
        
        # Initialize NLP models
        self.nlp = None
        self.classifier = None
        
        if spacy:
            try:
                self.nlp = spacy.load("en_core_web_sm")
            except Exception as e:
                print(f"Warning: Could not load spaCy model: {e}")

        if pipeline:
            try:
                # Zero-shot classification for intent
                # Using DistilBERT as requested
                self.classifier = pipeline("zero-shot-classification", model="typeform/distilbert-base-uncased-mnli")
            except Exception as e:
                print(f"Warning: Could not load Transformers pipeline: {e}")

    def parse(self, user_text: str, image_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Parses the user text and optional image to extract structured intent.
        """
        request_id = str(uuid.uuid4())
        user_id = "user_123" # Mocked user ID
        
        # 1. OCR & Image2SMILES Processing
        image_smiles = None
        if image_path:
            # Try OCR first
            if pytesseract and Image:
                try:
                    image_text = pytesseract.image_to_string(Image.open(image_path))
                    if image_text.strip():
                        user_text += f"\n[Image Content]: {image_text}"
                except Exception as e:
                    print(f"OCR Error: {e}")
            
            # Try Image2SMILES (DECIMER)
            if predict_SMILES:
                try:
                    image_smiles = predict_SMILES(image_path)
                    print(f"Image2SMILES detected: {image_smiles}")
                except Exception as e:
                    print(f"Image2SMILES Error: {e}")
            else:
                # Stub/Placeholder if DECIMER is not installed
                print("Image2SMILES (DECIMER) not installed. Skipping structure recognition.")

        intent_data = {
            "request_id": request_id,
            "user_id": user_id,
            "target": image_smiles, # Prioritize image-derived SMILES if found
            "constraints": [],
            "mode": "automation", # default
            "priority": "normal",
            "original_text": user_text,
            "timestamp": time.time(),
            "intent_type": "synthesis" # default
        }

        # 2. Intent Classification (DistilBERT)
        if self.classifier:
            labels = ["synthesis", "analysis", "procurement", "safety_check"]
            result = self.classifier(user_text, labels)
            intent_data["intent_type"] = result["labels"][0]

        # 3. Entity Extraction (spaCy + Regex)
        target = None
        
        # Regex (still useful for specific patterns)
        target_patterns = [
            r"synthesis (?:of|for)\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|on|in|$)",
            r"make\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|on|in|$)",
            r"prepare\s+([a-zA-Z0-9\-\[\]\(\)\s]+?)(?:\.|with|using|from|cheap|safe|on|in|$)"
        ]
        
        for pattern in target_patterns:
            match = re.search(pattern, user_text, re.IGNORECASE)
            if match:
                target = match.group(1).strip()
                break
        
        # spaCy NER fallback/refinement
        if not target and self.nlp:
            doc = self.nlp(user_text)
            # Look for entities that might be chemicals (often labeled as ORG, PRODUCT, or unknown in generic models)
            potential_targets = [ent.text for ent in doc.ents if ent.label_ in ["ORG", "PRODUCT", "WORK_OF_ART"]]
            if potential_targets:
                target = potential_targets[0]
            # Fallback: Look for the last NOUN in the sentence if no entities found
            elif not target:
                nouns = [chunk.text for chunk in doc.noun_chunks]
                # Filter out common non-chemical nouns
                stop_nouns = ["synthesis", "preparation", "method", "process", "way", "route", "it", "me", "analysis", "check"]
                candidates = [n for n in nouns if n.lower() not in stop_nouns]
                if candidates:
                    target = candidates[-1] # Assume the last noun is the target (e.g. "synthesis for [Paracetamol]")

        # Clean up target
        if target:
            stop_words = ["cheaply", "safely", "quickly", "it", "please", "me", "the", "a", "an"]
            for word in stop_words:
                if target.lower().endswith(f" {word}"):
                    target = target[:-len(word)-1]
            intent_data["target"] = target
            
        # Fallback for short queries (likely just the chemical name)
        if not intent_data.get("target"):
             words = user_text.split()
             if len(words) < 8: # Increased from 5 to 8 to catch "i need a safe synthesis for X"
                 # Heuristic: Take the last word(s) if it looks like a chemical name
                 # But for now, let's just trust the user input if it's short enough
                 # Remove common prefixes
                 clean_text = user_text.lower()
                 for prefix in ["i need", "synthesis of", "synthesis for", "make", "prepare", "how to make"]:
                     clean_text = clean_text.replace(prefix, "")
                 intent_data["target"] = clean_text.strip()

        # 4. Extract Constraints & Mode (Keyword/NER based)
        lower_text = user_text.lower()
        
        if "cheap" in lower_text or "low cost" in lower_text:
            intent_data["constraints"].append("low_cost")
        if "safe" in lower_text or "non-toxic" in lower_text:
            intent_data["constraints"].append("high_safety")
        if "fast" in lower_text or "quick" in lower_text:
            intent_data["constraints"].append("fast_delivery")

        if "bench" in lower_text or "manual" in lower_text:
            intent_data["mode"] = "bench"
        elif "automation" in lower_text or "robot" in lower_text:
            intent_data["mode"] = "automation"

        if "urgent" in lower_text or "asap" in lower_text:
            intent_data["priority"] = "high"

        # 5. LLM Refinement (GPT-4o-mini / LLaMA Stub)
        # This would be where the actual API call happens
        intent_data = self.refine_with_llm(user_text, intent_data)

        # 6. Validate and Save
        if self.validate_request(intent_data):
            self.save_request(intent_data)
        else:
            intent_data["error"] = "Could not identify target molecule."

        return intent_data

    def refine_with_llm(self, text: str, current_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Simulates a call to GPT-4o-mini or LLaMA 3.2 to refine slot filling and reasoning.
        """
        # Prompt Structure for LLM
        prompt = f"""
        You are an expert chemical assistant. Extract the following from the user query:
        - Target Molecule
        - Constraints (cost, safety, time)
        - Mode (bench vs automation)
        - Priority (normal vs high)
        
        User Query: "{text}"
        
        Current Extraction: {json.dumps(current_data, indent=2)}
        
        Return JSON only.
        """
        
        # Placeholder for API call
        # response = openai.ChatCompletion.create(model="gpt-4o-mini", messages=[...])
        # or
        # response = llama_pipeline(prompt)
        
        # For now, we return the data as is, but this structure is ready for the API.
        return current_data

    def validate_request(self, intent_data: Dict[str, Any]) -> bool:
        if not intent_data.get("target"):
            return False
        return True

    def save_request(self, intent_data: Dict[str, Any]):
        self.request_db.append(intent_data)

if __name__ == "__main__":
    agent = IntentParserAgent()
    sample_text = "I need the synthesis of Aspirin. It should be cheap and safe."
    print(json.dumps(agent.parse(sample_text), indent=2))
