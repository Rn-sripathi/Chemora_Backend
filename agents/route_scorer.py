try:
    import lightgbm as lgb
    lightgbm_available = True
except ImportError:
    lightgbm_available = False

try:
    import shap
    shap_available = True
except ImportError:
    shap_available = False

from typing import List, Dict, Any, Optional
import json
import os

class RouteScorerAgent:
    """
    Agent 7: Route Scorer / Ranker
    Ranks synthetic routes using multi-criteria optimization:
    - Rule-based scoring (weighted formula)
    - LightGBM ranker (learns from chemist feedback)
    - SHAP explanations for interpretability
    Dataset: Chemist feedback (grows with use)
    """
    def __init__(self, config_path: str = "scoring_config.json"):
        self.config_path = config_path
        
        # Load or initialize scoring weights (stored in DB/file)
        self.scoring_weights = self._load_scoring_config()
        
        # Initialize LightGBM ranker
        self.lgb_ranker = None
        if lightgbm_available:
            try:
                # In production, would load pre-trained model
                # self.lgb_ranker = lgb.Booster(model_file='route_ranker.txt')
                print("LightGBM Ranker: Ready (Mock)")
            except Exception as e:
                print(f"LightGBM Init Error: {e}")
        
        # Initialize SHAP explainer
        self.shap_explainer = None
        if shap_available and self.lgb_ranker:
            try:
                # self.shap_explainer = shap.TreeExplainer(self.lgb_ranker)
                print("SHAP Explainer: Ready (Mock)")
            except Exception as e:
                print(f"SHAP Init Error: {e}")
        
        # Chemist feedback dataset (grows with use)
        self.feedback_dataset = []
        self._load_feedback_data()

    def _load_scoring_config(self) -> Dict[str, float]:
        """
        Loads weighted scoring formula from configuration file/DB.
        """
        default_weights = {
            "yield": 0.25,          # Higher yield is better
            "cost": 0.20,           # Lower cost is better
            "safety": 0.20,         # Higher safety score is better
            "complexity": 0.15,     # Lower complexity is better
            "green_chemistry": 0.10, # Higher sustainability is better
            "availability": 0.10    # Higher reagent availability is better
        }
        
        try:
            if os.path.exists(self.config_path):
                with open(self.config_path, 'r') as f:
                    weights = json.load(f)
                    return weights
        except Exception as e:
            print(f"Config load error: {e}, using defaults")
        
        return default_weights

    def _save_scoring_config(self):
        """
        Saves current scoring weights to configuration file.
        """
        try:
            with open(self.config_path, 'w') as f:
                json.dump(self.scoring_weights, f, indent=2)
        except Exception as e:
            print(f"Config save error: {e}")

    def _load_feedback_data(self):
        """
        Loads chemist feedback dataset (grows with use).
        """
        try:
            if os.path.exists("chemist_feedback.json"):
                with open("chemist_feedback.json", 'r') as f:
                    self.feedback_dataset = json.load(f)
                    print(f"Loaded {len(self.feedback_dataset)} feedback entries")
        except Exception as e:
            print(f"Feedback load error: {e}")

    def _save_feedback_data(self):
        """
        Saves feedback dataset.
        """
        try:
            with open("chemist_feedback.json", 'w') as f:
                json.dump(self.feedback_dataset, f, indent=2)
        except Exception as e:
            print(f"Feedback save error: {e}")

    def calculate_rule_based_score(self, route: Dict[str, Any]) -> float:
        """
        Rule-based scoring using weighted multi-criteria optimization.
        Formula: Score = Σ(weight_i × normalized_criterion_i)
        """
        # Extract route metrics
        avg_yield = route.get("confidence", 0.5) * 100  # Use confidence as proxy
        cost_score = 1.0 - (route.get("cost", {}).get("total_estimated_cost", 50) / 100.0)
        safety_score = route.get("safety", {}).get("overall_risk", "medium")
        
        # Normalize safety score
        safety_map = {"low": 1.0, "medium": 0.6, "high": 0.2}
        safety_normalized = safety_map.get(safety_score, 0.5)
        
        # Calculate complexity (inverse of number of steps)
        num_steps = len(route.get("steps", []))
        complexity_score = 1.0 / max(num_steps, 1) if num_steps > 0 else 0.5
        
        # Green chemistry score (placeholder)
        green_score = 0.7
        
        # Availability score (placeholder)
        availability_score = 0.8
        
        # Apply weighted formula
        score = (
            self.scoring_weights["yield"] * (avg_yield / 100.0) +
            self.scoring_weights["cost"] * cost_score +
            self.scoring_weights["safety"] * safety_normalized +
            self.scoring_weights["complexity"] * complexity_score +
            self.scoring_weights["green_chemistry"] * green_score +
            self.scoring_weights["availability"] * availability_score
        )
        
        return round(score * 100, 2)  # Return as percentage

    def rank_with_lightgbm(self, routes: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Uses LightGBM ranker to learn from chemist feedback and rank routes.
        """
        if not self.lgb_ranker or not lightgbm_available:
            # Fallback to rule-based
            return routes
        
        # Extract features from routes
        # features = self._extract_features(routes)
        # predictions = self.lgb_ranker.predict(features)
        
        # For now, return as-is
        return routes

    def explain_with_shap(self, route: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Generates SHAP explanations for route ranking.
        Shows which factors contributed most to the score.
        """
        if not self.shap_explainer or not shap_available:
            return None
        
        # Extract features
        # features = self._extract_features([route])
        # shap_values = self.shap_explainer.shap_values(features)
        
        # Mock explanation
        return {
            "yield_contribution": 0.25,
            "cost_contribution": 0.20,
            "safety_contribution": 0.18,
            "complexity_contribution": 0.15,
            "green_chemistry_contribution": 0.12,
            "availability_contribution": 0.10,
            "note": "SHAP analysis shows yield and cost are the primary drivers of this route's score."
        }

    def add_feedback(self, route_id: str, chemist_rating: float, comments: str = ""):
        """
        Adds chemist feedback to the growing dataset.
        This data is used to retrain the LightGBM ranker.
        """
        feedback_entry = {
            "route_id": route_id,
            "chemist_rating": chemist_rating,  # 1-5 scale
            "comments": comments,
            "timestamp": "2024-01-01"  # Would use actual timestamp
        }
        
        self.feedback_dataset.append(feedback_entry)
        self._save_feedback_data()
        
        print(f"Feedback added. Dataset now has {len(self.feedback_dataset)} entries.")
        
        # Trigger retraining if enough data
        if len(self.feedback_dataset) >= 100:
            self._retrain_lgb_ranker()

    def _retrain_lgb_ranker(self):
        """
        Retrains LightGBM ranker with accumulated chemist feedback.
        """
        if not lightgbm_available:
            return
        
        print("Retraining LightGBM ranker with chemist feedback...")
        # Would extract features and labels from feedback_dataset
        # train_data = lgb.Dataset(features, labels)
        # params = {'objective': 'lambdarank', 'metric': 'ndcg'}
        # self.lgb_ranker = lgb.train(params, train_data)

    def rank(self, routes: List[Dict[str, Any]], user_preferences: Dict[str, float] = None) -> List[Dict[str, Any]]:
        """
        Ranks routes using multi-criteria optimization.
        Combines rule-based scoring and LightGBM ranker.
        
        Args:
            routes: List of synthetic routes
            user_preferences: Optional weight adjustments from user
        
        Returns:
            Ranked routes with scores and explanations
        """
        if not routes:
            return []
        
        # Apply user preference weights if provided
        if user_preferences:
            original_weights = self.scoring_weights.copy()
            self.scoring_weights.update(user_preferences)
        
        # Calculate rule-based scores
        for route in routes:
            rule_score = self.calculate_rule_based_score(route)
            route["score"] = rule_score
            route["scoring_method"] = "Rule-based (weighted)"
            route["weights_used"] = self.scoring_weights.copy()
        
        # Apply LightGBM ranking if available
        if self.lgb_ranker and len(self.feedback_dataset) >= 10:
            routes = self.rank_with_lightgbm(routes)
            for route in routes:
                route["scoring_method"] = "Hybrid (Rule-based + LightGBM)"
        
        # Sort by score
        ranked_routes = sorted(routes, key=lambda x: x.get("score", 0), reverse=True)
        
        # Add SHAP explanations to top routes
        for i, route in enumerate(ranked_routes[:3]):
            if self.shap_explainer:
                route["shap_explanation"] = self.explain_with_shap(route)
        
        # Restore original weights if modified
        if user_preferences:
            self.scoring_weights = original_weights
        
        return ranked_routes

    def update_weights(self, new_weights: Dict[str, float]):
        """
        Updates scoring weights based on user/chemist preferences.
        """
        self.scoring_weights.update(new_weights)
        self._save_scoring_config()
        print(f"Scoring weights updated: {self.scoring_weights}")

    def score(self, routes: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Backwards-compatible method for graph.py integration.
        Scores and ranks routes using the new rank() method.
        
        Args:
            routes: List of routes to score
        
        Returns:
            Ranked routes with scores
        """
        return self.rank(routes)
