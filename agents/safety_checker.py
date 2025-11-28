from typing import List, Dict, Any

class SafetyCheckerAgent:
    """
    Agent 8: Safety & Compliance Checker
    Evaluates hazards and regulatory constraints.
    """
    def __init__(self):
        self.hazard_db = {
            "benzene": "Carcinogen (Category 1A)",
            "sodium cyanide": "Acute Toxicity (Category 1)",
            "picric acid": "Explosive",
            "acetic anhydride": "Flammable Liquid, Corrosive",
            "sulfuric acid": "Skin Corrosion (Category 1A)",
            "acetyl chloride": "Reacts violently with water"
        }

    def check(self, route: Dict[str, Any]) -> Dict[str, Any]:
        """
        Checks a route for safety hazards.
        """
        hazards = []
        risk_level = "Low"

        for step in route.get("steps", []):
            reagents = step.get("reagents", [])
            for reagent in reagents:
                reagent_lower = reagent.lower()
                for haz_chem, haz_desc in self.hazard_db.items():
                    if haz_chem in reagent_lower:
                        hazards.append({
                            "reagent": reagent,
                            "hazard": haz_desc,
                            "step": step.get("reaction")
                        })
                        if "Explosive" in haz_desc or "Acute Toxicity" in haz_desc:
                            risk_level = "High"
                        elif risk_level != "High":
                            risk_level = "Medium"

        return {
            "hazards": hazards,
            "risk_level": risk_level,
            "compliant": risk_level != "High" # Simple rule
        }
