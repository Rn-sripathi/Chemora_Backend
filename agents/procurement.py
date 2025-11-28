from typing import List, Dict, Any
import random

class ProcurementAgent:
    """
    Agent 10: Procurement & Cost Estimator
    Estimates reagent availability and cost.
    """
    def __init__(self):
        self.catalog = {
            "salicylic acid": {"price": 25.0, "unit": "500g", "vendor": "Sigma"},
            "acetic anhydride": {"price": 45.0, "unit": "1L", "vendor": "Fisher"},
            "isobutylbenzene": {"price": 30.0, "unit": "100g", "vendor": "TCI"},
            "sodium borohydride": {"price": 55.0, "unit": "100g", "vendor": "Sigma"}
        }

    def estimate(self, route: Dict[str, Any]) -> Dict[str, Any]:
        """
        Estimates total cost for a route.
        """
        total_cost = 0.0
        details = []
        
        for step in route.get("steps", []):
            for reagent in step.get("reagents", []):
                reagent_lower = reagent.lower()
                # Simple partial match
                found = False
                for item, info in self.catalog.items():
                    if item in reagent_lower:
                        cost = info["price"]
                        total_cost += cost
                        details.append({
                            "reagent": reagent,
                            "cost": cost,
                            "vendor": info["vendor"],
                            "status": "In Stock"
                        })
                        found = True
                        break
                
                if not found:
                    # Mock cost for unknown items
                    mock_cost = round(random.uniform(10.0, 100.0), 2)
                    total_cost += mock_cost
                    details.append({
                        "reagent": reagent,
                        "cost": mock_cost,
                        "vendor": "Generic Supplier",
                        "status": "Low Stock"
                    })

        return {
            "total_estimated_cost": round(total_cost, 2),
            "currency": "USD",
            "breakdown": details
        }
