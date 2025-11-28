from typing import List, Dict, Any

class RetrosynthesisAgent:
    """
    Agent 5: Retrosynthesis Planner
    Proposes synthetic routes from starting materials to target.
    """
    def __init__(self):
        # Mock database of routes
        self.mock_routes = {
            "aspirin": [
                {
                    "id": "route_1",
                    "steps": [
                        {"reaction": "Salicylic acid + Acetic anhydride -> Aspirin + Acetic acid", "reagents": ["Salicylic acid", "Acetic anhydride"], "conditions": "H2SO4 catalyst, 80C"}
                    ],
                    "confidence": 0.95
                },
                {
                    "id": "route_2",
                    "steps": [
                        {"reaction": "Salicylic acid + Acetyl chloride -> Aspirin + HCl", "reagents": ["Salicylic acid", "Acetyl chloride"], "conditions": "Pyridine, 0C"}
                    ],
                    "confidence": 0.85
                }
            ],
            "ibuprofen": [
                {
                    "id": "route_1",
                    "steps": [
                        {"reaction": "Isobutylbenzene + Acetic anhydride -> 4-Isobutylacetophenone", "reagents": ["Isobutylbenzene", "Acetic anhydride"], "conditions": "Friedel-Crafts"},
                        {"reaction": "4-Isobutylacetophenone -> Ibuprofen", "reagents": ["NaBH4", "CO", "Pd catalyst"], "conditions": "Reduction/Carbonylation"}
                    ],
                    "confidence": 0.90
                }
            ]
        }

    def plan(self, target_name: str, constraints: List[str] = None) -> List[Dict[str, Any]]:
        """
        Generates retrosynthetic routes.
        """
        key = target_name.lower()
        if key in self.mock_routes:
            return self.mock_routes[key]
        
        # Generic fallback route
        return [{
            "id": "route_generic",
            "steps": [
                {"reaction": f"Precursor A + Precursor B -> {target_name}", "reagents": ["Precursor A", "Precursor B"], "conditions": "Standard Conditions"}
            ],
            "confidence": 0.5
        }]
