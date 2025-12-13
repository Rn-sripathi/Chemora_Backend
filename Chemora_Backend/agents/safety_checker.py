try:
    from rdkit import Chem
except ImportError:
    Chem = None

try:
    import pubchempy as pcp
except ImportError:
    pcp = None

# Optional: ChemBERTa for toxicity prediction
try:
    from transformers import AutoTokenizer, AutoModelForSequenceClassification
    chemberta_available = True
except ImportError:
    chemberta_available = False

import requests
from typing import List, Dict, Any, Optional
import re

class SafetyCheckerAgent:
    """
    Agent 8: Safety & Compliance Checker
    Evaluates hazards and regulatory constraints using:
    - GHS (Globally Harmonized System) hazard database
    - PubChem Safety API
    - ECHA (European Chemicals Agency) database
    - Python rule engine for compliance
    - ChemBERTa-Tox21 for ML-based toxicity prediction
    """
    def __init__(self):
        # Initialize ChemBERTa-Tox21 model
        self.chemberta_model = None
        self.chemberta_tokenizer = None
        
        if chemberta_available:
            try:
                # Could load "seyonec/ChemBERTa-zinc-base-v1" or "DeepChem/ChemBERTa-77M-MLM"
                # For Tox21: would use fine-tuned version
                # self.chemberta_tokenizer = AutoTokenizer.from_pretrained("...")
                # self.chemberta_model = AutoModelForSequenceClassification.from_pretrained("...")
                print("ChemBERTa-Tox21: Ready (Mock)")
            except Exception as e:
                print(f"ChemBERTa Init Error: {e}")
        
        # GHS Hazard Database (mock - would be actual database)
        self.ghs_hazards = {
            "H200": {"code": "H200", "hazard": "Unstable explosive", "category": "Physical"},
            "H225": {"code": "H225", "hazard": "Highly flammable liquid and vapour", "category": "Physical"},
            "H301": {"code": "H301", "hazard": "Toxic if swallowed", "category": "Health"},
            "H314": {"code": "H314", "hazard": "Causes severe skin burns and eye damage", "category": "Health"},
            "H315": {"code": "H315", "hazard": "Causes skin irritation", "category": "Health"},
            "H317": {"code": "H317", "hazard": "May cause allergic skin reaction", "category": "Health"},
            "H318": {"code": "H318", "hazard": "Causes serious eye damage", "category": "Health"},
            "H350": {"code": "H350", "hazard": "May cause cancer", "category": "Health"},
            "H400": {"code": "H400", "hazard": "Very toxic to aquatic life", "category": "Environmental"},
            "H410": {"code": "H410", "hazard": "Very toxic to aquatic life with long lasting effects", "category": "Environmental"}
        }
        
        # Python rule engine (simple rule-based system)
        self.safety_rules = self._initialize_safety_rules()

    def _initialize_safety_rules(self) -> List[Dict[str, Any]]:
        """
        Initializes Python-based rule engine for safety compliance.
        """
        return [
            {
                "rule_id": "RULE_001",
                "name": "Explosive compounds check",
                "condition": lambda chem: any(x in chem.lower() for x in ["nitro", "azide", "peroxide"]),
                "action": "FLAG_HIGH_RISK",
                "message": "Contains potentially explosive functional groups"
            },
            {
                "rule_id": "RULE_002",
                "name": "Carcinogenic compounds",
                "condition": lambda chem: any(x in chem.lower() for x in ["benzene", "asbestos", "formaldehyde"]),
                "action": "FLAG_CARCINOGEN",
                "message": "Known or suspected carcinogen"
            },
            {
                "rule_id": "RULE_003",
                "name": "Highly reactive metals",
                "condition": lambda chem: any(x in chem.lower() for x in ["sodium", "potassium", "lithium"]) and "metal" in chem.lower(),
                "action": "FLAG_REACTIVE",
                "message": "Highly reactive with water/air"
            },
            {
                "rule_id": "RULE_004",
                "name": "Toxic heavy metals",
                "condition": lambda chem: any(x in chem.lower() for x in ["mercury", "lead", "cadmium", "chromium"]),
                "action": "FLAG_TOXIC_METAL",
                "message": "Contains toxic heavy metal"
            }
        ]

    def query_pubchem_safety(self, compound_name: str) -> Optional[Dict[str, Any]]:
        """
        Queries PubChem Safety API for hazard information.
        """
        if not pcp:
            return None
        
        try:
            compounds = pcp.get_compounds(compound_name, 'name')
            if compounds:
                compound = compounds[0]
                
                # PubChem provides GHS classification
                return {
                    "cid": compound.cid,
                    "molecular_formula": compound.molecular_formula,
                    "ghs_hazards": self._extract_ghs_from_pubchem(compound),
                    "source": "PubChem Safety API"
                }
        except Exception as e:
            print(f"PubChem Safety API Error: {e}")
        
        return None

    def _extract_ghs_from_pubchem(self, compound) -> List[str]:
        """
        Extracts GHS hazard codes from PubChem compound data.
        """
        # In reality, would parse from compound.to_dict() or specific properties
        # Mock extraction for common compounds
        mock_ghs_data = {
            "aspirin": ["H302", "H315"],  # Harmful if swallowed, skin irritation
            "benzene": ["H225", "H350", "H340"],  # Flammable, carcinogen, mutagenic
            "acetic anhydride": ["H226", "H302", "H314"]  # Flammable, toxic, corrosive
        }
        
        compound_name = getattr(compound, 'iupac_name', '').lower()
        return mock_ghs_data.get(compound_name, [])

    def query_echa_database(self, compound_name: str) -> Optional[Dict[str, Any]]:
        """
        Queries ECHA (European Chemicals Agency) database for regulatory info.
        """
        # ECHA has a public API, would use actual endpoint
        # https://echa.europa.eu/information-on-chemicals
        
        try:
            # Mock ECHA query
            # In production: requests.get(f"https://echa.europa.eu/api/substance/{cas_number}")
            
            mock_echa_data = {
                "substance_name": compound_name,
                "reach_registered": True,
                "classification": ["Flam. Liq. 2", "Acute Tox. 4"],
                "hazard_statements": ["H225", "H302"],
                "regulatory_status": "Approved for use in EU",
                "restrictions": "None",
                "source": "ECHA Database"
            }
            
            return mock_echa_data
        except Exception as e:
            print(f"ECHA Database Error: {e}")
        
        return None

    def apply_rule_engine(self, compound_name: str, route_steps: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Applies Python rule engine to check safety compliance.
        """
        violations = []
        
        # Check compound name against rules
        for rule in self.safety_rules:
            if rule["condition"](compound_name):
                violations.append({
                    "rule_id": rule["rule_id"],
                    "rule_name": rule["name"],
                    "action": rule["action"],
                    "message": rule["message"],
                    "severity": "HIGH"
                })
        
        # Check all reagents in route steps
        for step in route_steps:
            reagents = step.get("reagents", [])
            for reagent in reagents:
                for rule in self.safety_rules:
                    if rule["condition"](reagent):
                        violations.append({
                            "rule_id": rule["rule_id"],
                            "rule_name": rule["name"],
                            "action": rule["action"],
                            "message": f"{rule['message']} (Reagent: {reagent})",
                            "severity": "MEDIUM",
                            "step": step.get("reaction", "Unknown")
                        })
        
        return violations

    def predict_toxicity_with_chemberta(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Uses ChemBERTa-Tox21 to predict toxicity endpoints.
        Tox21 includes 12 different toxicity assays.
        """
        if not self.chemberta_model or not chemberta_available:
            return None
        
        # Tokenize SMILES
        # inputs = self.chemberta_tokenizer(smiles, return_tensors="pt")
        # outputs = self.chemberta_model(**inputs)
        # predictions = torch.sigmoid(outputs.logits)
        
        # Mock prediction
        return {
            "nr_ar": 0.15,      # Nuclear receptor - Androgen receptor
            "nr_er": 0.12,      # Nuclear receptor - Estrogen receptor
            "sr_mmp": 0.08,     # Stress response - Mitochondrial membrane potential
            "overall_toxicity": "LOW",
            "confidence": 0.85,
            "model": "ChemBERTa-Tox21"
        }

    def check(self, route: Dict[str, Any]) -> Dict[str, Any]:
        """
        Comprehensive safety check using all available tools and datasets.
        Workflow:
        1. Query PubChem Safety API
        2. Query ECHA database
        3. Apply Python rule engine
        4. Predict toxicity with ChemBERTa
        5. Aggregate results with GHS hazard mapping
        """
        target_name = route.get("target_molecule", "Unknown")
        steps = route.get("steps", [])
        
        # 1. PubChem Safety API
        pubchem_data = self.query_pubchem_safety(target_name)
        
        # 2. ECHA Database
        echa_data = self.query_echa_database(target_name)
        
        # 3. Rule Engine
        rule_violations = self.apply_rule_engine(target_name, steps)
        
        # 4. ChemBERTa Toxicity Prediction (if SMILES available)
        toxicity_prediction = None
        target_smiles = route.get("target_smiles")
        if target_smiles:
            toxicity_prediction = self.predict_toxicity_with_chemberta(target_smiles)
        
        # 5. Map GHS hazards
        ghs_codes = []
        if pubchem_data:
            ghs_codes.extend(pubchem_data.get("ghs_hazards", []))
        if echa_data:
            ghs_codes.extend(echa_data.get("hazard_statements", []))
        
        ghs_details = [self.ghs_hazards.get(code, {"code": code, "hazard": "Unknown"}) for code in ghs_codes]
        
        # Determine overall risk level
        overall_risk = "low"
        if rule_violations or any("H3" in code for code in ghs_codes):  # H3xx are health hazards
            overall_risk = "high"
        elif ghs_codes:
            overall_risk = "medium"
        
        return {
            "overall_risk": overall_risk,
            "ghs_hazards": ghs_details,
            "pubchem_data": pubchem_data,
            "echa_data": echa_data,
            "rule_violations": rule_violations,
            "toxicity_prediction": toxicity_prediction,
            "recommendations": self._generate_recommendations(overall_risk, rule_violations),
            "data_sources": ["GHS Database", "PubChem Safety API", "ECHA", "ChemBERTa-Tox21"],
            "compliance_status": "PASS" if overall_risk == "low" else "REVIEW_REQUIRED"
        }

    def _generate_recommendations(self, risk_level: str, violations: List[Dict[str, Any]]) -> List[str]:
        """
        Generates safety recommendations based on assessment.
        """
        recommendations = []
        
        if risk_level == "high":
            recommendations.append("Conduct risk assessment before proceeding")
            recommendations.append("Ensure proper PPE (Personal Protective Equipment)")
            recommendations.append("Use fume hood for all operations")
        
        if violations:
            recommendations.append("Review regulatory compliance for flagged compounds")
        
        if risk_level == "medium":
            recommendations.append("Follow standard laboratory safety protocols")
        
        recommendations.append("Maintain Material Safety Data Sheets (MSDS) on file")
        
        return recommendations
