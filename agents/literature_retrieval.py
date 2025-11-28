import numpy as np
try:
    from sentence_transformers import SentenceTransformer
except ImportError:
    SentenceTransformer = None

from typing import List, Dict, Any

class LiteratureRetrievalAgent:
    """
    Agent 3: Literature Retrieval
    Finds papers, patents, protocols, and extracts reaction snippets.
    """
    def __init__(self):
        # Load a pre-trained model for embeddings
        # In a real app, this would be hosted or cached
        try:
            if SentenceTransformer:
                self.model = SentenceTransformer('all-MiniLM-L6-v2')
            else:
                self.model = None
        except:
            print("Warning: Could not load SentenceTransformer. Using mock embeddings.")
            self.model = None

        # Mock Vector DB
        self.documents = [
            {"id": "doc1", "text": "Synthesis of Aspirin using salicylic acid and acetic anhydride.", "yield": "85%"},
            {"id": "doc2", "text": "Preparation of Ibuprofen via Friedel-Crafts acylation.", "yield": "90%"},
            {"id": "doc3", "text": "Green synthesis of Aspirin with microwave irradiation.", "yield": "92%"},
        ]
        self.embeddings = self._compute_embeddings([d["text"] for d in self.documents])

    def _compute_embeddings(self, texts: List[str]) -> np.ndarray:
        if self.model:
            return self.model.encode(texts)
        else:
            return np.random.rand(len(texts), 384) # Mock 384-dim embeddings

    def search(self, query: str, top_k: int = 3) -> List[Dict[str, Any]]:
        """
        Retrieves relevant documents based on semantic similarity.
        """
        if not self.model:
            return self.documents[:top_k]

        query_embedding = self.model.encode([query])
        
        # Cosine similarity
        similarities = np.dot(self.embeddings, query_embedding.T).flatten()
        
        # Sort by similarity
        top_indices = np.argsort(similarities)[::-1][:top_k]
        
        results = []
        for idx in top_indices:
            doc = self.documents[idx].copy()
            doc["score"] = float(similarities[idx])
            results.append(doc)
            
        return results

if __name__ == "__main__":
    agent = LiteratureRetrievalAgent()
    results = agent.search("How to make Aspirin?")
    print(results)
