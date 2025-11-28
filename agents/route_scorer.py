from typing import List, Dict, Any

class RouteScorerAgent:
    """
    Agent 7: Route Scorer / Ranker
    Ranks routes based on multi-criteria optimization.
    """
    def __init__(self):
        self.weights = {
            "yield": 0.5,
            "cost": 0.3, # Negative weight logic applied later
            "safety": 0.2
        }

    def score(self, routes: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Scores and ranks the provided routes.
        """
        scored_routes = []
        for route in routes:
            # Calculate aggregate metrics from steps
            total_yield = 1.0
            step_count = len(route.get("steps", []))
            
            # Mock metrics if not present
            for step in route.get("steps", []):
                step_yield = step.get("prediction", {}).get("predicted_yield", 90.0) / 100.0
                total_yield *= step_yield

            # Simple score: Yield / Steps (penalize long routes)
            # In reality, this would be complex
            raw_score = (total_yield * 100) - (step_count * 5)
            
            route_copy = route.copy()
            route_copy["score"] = round(raw_score, 2)
            scored_routes.append(route_copy)

        # Sort by score descending
        scored_routes.sort(key=lambda x: x["score"], reverse=True)
        return scored_routes
