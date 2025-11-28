from typing import Dict, Any
import random

class ForwardPredictionAgent:
    """
    Agent 6: Forward Reaction Predictor
    Predicts yield, selectivity, and feasibility.
    """
    def __init__(self):
        pass

    def predict(self, reaction_step: Dict[str, Any]) -> Dict[str, Any]:
        """
        Predicts outcome for a single reaction step.
        """
        # Mock logic: Random yield between 50% and 99%
        # In reality, this would use a ML model (e.g., YieldBERT)
        predicted_yield = round(random.uniform(50.0, 99.0), 1)
        
        # Feasibility check (mock)
        feasibility = "High"
        if predicted_yield < 60:
            feasibility = "Low"
        elif predicted_yield < 80:
            feasibility = "Medium"

        return {
            "predicted_yield": predicted_yield,
            "feasibility": feasibility,
            "side_products": ["Trace impurities"] if predicted_yield > 90 else ["Significant byproducts"]
        }
