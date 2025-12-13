try:
    from pymongo import MongoClient
    mongodb_available = True
except ImportError:
    mongodb_available = False
    MongoClient = None

try:
    import mlflow
    mlflow_available = True
except ImportError:
    mlflow_available = False

try:
    import boto3
    s3_available = True
except ImportError:
    s3_available = False
    boto3 = None

try:
    from qdrant_client import QdrantClient
    from qdrant_client.models import Distance, VectorParams, PointStruct
    qdrant_available = True
except ImportError:
    qdrant_available = False
    QdrantClient = None

from typing import List, Dict, Any, Optional
from datetime import datetime
import json
import os
import uuid

class DataCurationAgent:
    """
    Agent 11: Data Curation & Provenance Logger
    Tracks data lineage and stores research artifacts using:
    - MongoDB for append-only provenance logs
    - MLflow for model versioning and experiment tracking
    - MongoDB/S3 for storing PDFs and log files
    - Vector DB (Qdrant/Pinecone) for embeddings storage
    """
    def __init__(self, 
                 mongo_uri: str = "mongodb://localhost:27017/",
                 mlflow_uri: str = "sqlite:///mlflow.db",
                 s3_bucket: str = "chemora-documents",
                 vector_db_url: str = "http://localhost:6333"):
        
        # Initialize MongoDB for provenance logs
        self.mongo_client = None
        self.provenance_collection = None
        self.documents_collection = None
        
        if mongodb_available:
            try:
                self.mongo_client = MongoClient(mongo_uri)
                db = self.mongo_client["chemora"]
                self.provenance_collection = db["provenance_logs"]
                self.documents_collection = db["documents"]
                
                # Create indexes for efficient querying
                self.provenance_collection.create_index("timestamp")
                self.provenance_collection.create_index("user_id")
                self.provenance_collection.create_index("experiment_id")
                
                print("MongoDB: Connected for provenance logging")
            except Exception as e:
                print(f"MongoDB Init Error: {e}")
        
        # Initialize MLflow for model versioning
        if mlflow_available:
            try:
                mlflow.set_tracking_uri(mlflow_uri)
                print(f"MLflow: Tracking URI set to {mlflow_uri}")
            except Exception as e:
                print(f"MLflow Init Error: {e}")
        
        # Initialize S3 for document storage
        self.s3_client = None
        self.s3_bucket = s3_bucket
        
        if s3_available:
            try:
                self.s3_client = boto3.client('s3')
                # Create bucket if it doesn't exist
                # self.s3_client.create_bucket(Bucket=s3_bucket)
                print(f"S3: Ready for document storage (bucket: {s3_bucket})")
            except Exception as e:
                print(f"S3 Init Error: {e}")
        
        # Initialize Vector DB for embeddings
        self.vector_client = None
        self.collection_name = "chemical_embeddings"
        
        if qdrant_available:
            try:
                self.vector_client = QdrantClient(url=vector_db_url)
                # Create collection if it doesn't exist
                # self.vector_client.create_collection(
                #     collection_name=self.collection_name,
                #     vectors_config=VectorParams(size=768, distance=Distance.COSINE)
                # )
                print(f"Qdrant Vector DB: Connected at {vector_db_url}")
            except Exception as e:
                print(f"Qdrant Init Error: {e}")

    def log_provenance(self, event_type: str, data: Dict[str, Any], user_id: str = "user_123") -> str:
        """
        Logs provenance event to MongoDB (append-only).
        
        Args:
            event_type: Type of event (e.g., 'synthesis_request', 'route_generated')
            data: Event data
            user_id: User identifier
        
        Returns:
            Event ID
        """
        event_id = str(uuid.uuid4())
        
        provenance_entry = {
            "event_id": event_id,
            "event_type": event_type,
            "timestamp": datetime.utcnow().isoformat(),
            "user_id": user_id,
            "data": data,
            "version": "1.0"
        }
        
        if self.provenance_collection is not None:
            try:
                result = self.provenance_collection.insert_one(provenance_entry)
                print(f"Provenance logged: {event_id} (MongoDB ID: {result.inserted_id})")
                return event_id
            except Exception as e:
                print(f"MongoDB logging error: {e}")
        else:
            print(f"Provenance logged (mock): {event_id}")
        
        return event_id

    def get_provenance_chain(self, experiment_id: str) -> List[Dict[str, Any]]:
        """
        Retrieves complete provenance chain for an experiment.
        """
        if self.provenance_collection is None:
            return []
        
        try:
            events = list(self.provenance_collection.find(
                {"data.experiment_id": experiment_id},
                {"_id": 0}  # Exclude MongoDB internal ID
            ).sort("timestamp", 1))
            
            return events
        except Exception as e:
            print(f"Provenance retrieval error: {e}")
            return []

    def log_model_version(self, model_name: str, model_params: Dict[str, Any], 
                         metrics: Dict[str, float], artifacts_path: Optional[str] = None) -> str:
        """
        Logs model version using MLflow.
        
        Args:
            model_name: Name of the model
            model_params: Hyperparameters and configuration
            metrics: Performance metrics
            artifacts_path: Path to model artifacts
        
        Returns:
            Run ID
        """
        if not mlflow_available:
            return "mlflow_unavailable"
        
        try:
            with mlflow.start_run(run_name=model_name):
                # Log parameters
                for key, value in model_params.items():
                    mlflow.log_param(key, value)
                
                # Log metrics
                for key, value in metrics.items():
                    mlflow.log_metric(key, value)
                
                # Log artifacts if provided
                if artifacts_path and os.path.exists(artifacts_path):
                    mlflow.log_artifacts(artifacts_path)
                
                # Get run ID
                run = mlflow.active_run()
                run_id = run.info.run_id
                
                print(f"MLflow: Model '{model_name}' logged (Run ID: {run_id})")
                return run_id
                
        except Exception as e:
            print(f"MLflow logging error: {e}")
            return "error"

    def store_document_mongodb(self, document_name: str, content: bytes, 
                               metadata: Dict[str, Any]) -> str:
        """
        Stores document (PDF, log file) in MongoDB GridFS.
        """
        if self.documents_collection is None:
            return "mongodb_unavailable"
        
        try:
            document = {
                "document_id": str(uuid.uuid4()),
                "name": document_name,
                "content": content,  # Store as binary
                "metadata": metadata,
                "uploaded_at": datetime.utcnow().isoformat()
            }
            
            result = self.documents_collection.insert_one(document)
            print(f"Document stored in MongoDB: {document_name}")
            return str(result.inserted_id)
            
        except Exception as e:
            print(f"MongoDB document storage error: {e}")
            return "error"

    def store_document_s3(self, document_name: str, file_path: str, 
                         metadata: Dict[str, Any]) -> str:
        """
        Stores document in S3.
        
        Args:
            document_name: Name of the document
            file_path: Local path to file
            metadata: Document metadata
        
        Returns:
            S3 object key
        """
        if self.s3_client is None:
            return "s3_unavailable"
        
        try:
            object_key = f"documents/{datetime.utcnow().strftime('%Y/%m/%d')}/{document_name}"
            
            # Upload file
            self.s3_client.upload_file(
                file_path,
                self.s3_bucket,
                object_key,
                ExtraArgs={"Metadata": metadata}
            )
            
            print(f"Document stored in S3: s3://{self.s3_bucket}/{object_key}")
            return object_key
            
        except Exception as e:
            print(f"S3 storage error: {e}")
            return "error"

    def retrieve_document_s3(self, object_key: str, download_path: str) -> bool:
        """
        Retrieves document from S3.
        """
        if self.s3_client is None:
            return False
        
        try:
            self.s3_client.download_file(
                self.s3_bucket,
                object_key,
                download_path
            )
            print(f"Document retrieved from S3: {object_key}")
            return True
            
        except Exception as e:
            print(f"S3 retrieval error: {e}")
            return False

    def store_embedding(self, embedding_id: str, vector: List[float], 
                       payload: Dict[str, Any]) -> bool:
        """
        Stores embedding in Vector DB (Qdrant).
        
        Args:
            embedding_id: Unique identifier
            vector: Embedding vector
            payload: Associated metadata
        
        Returns:
            Success status
        """
        if self.vector_client is None:
            return False
        
        try:
            point = PointStruct(
                id=embedding_id,
                vector=vector,
                payload=payload
            )
            
            self.vector_client.upsert(
                collection_name=self.collection_name,
                points=[point]
            )
            
            print(f"Embedding stored in Vector DB: {embedding_id}")
            return True
            
        except Exception as e:
            print(f"Vector DB storage error: {e}")
            return False

    def search_embeddings(self, query_vector: List[float], limit: int = 10) -> List[Dict[str, Any]]:
        """
        Searches for similar embeddings in Vector DB.
        """
        if self.vector_client is None:
            return []
        
        try:
            results = self.vector_client.search(
                collection_name=self.collection_name,
                query_vector=query_vector,
                limit=limit
            )
            
            return [
                {
                    "id": hit.id,
                    "score": hit.score,
                    "payload": hit.payload
                }
                for hit in results
            ]
            
        except Exception as e:
            print(f"Vector DB search error: {e}")
            return []

    def curate_experiment(self, experiment_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Comprehensive experiment curation workflow.
        
        Workflow:
        1. Log provenance to MongoDB
        2. Store documents to S3/MongoDB
        3. Log model version to MLflow
        4. Store embeddings to Vector DB
        5. Return curation report
        """
        experiment_id = experiment_data.get("experiment_id", str(uuid.uuid4()))
        
        curation_report = {
            "experiment_id": experiment_id,
            "curation_timestamp": datetime.utcnow().isoformat(),
            "provenance_id": None,
            "model_run_id": None,
            "documents_stored": [],
            "embeddings_stored": 0,
            "status": "success"
        }
        
        # 1. Log provenance
        provenance_id = self.log_provenance(
            event_type="experiment_curated",
            data={
                "experiment_id": experiment_id,
                "target_molecule": experiment_data.get("target_molecule"),
                "route_selected": experiment_data.get("route_id")
            }
        )
        curation_report["provenance_id"] = provenance_id
        
        # 2. Store documents (if any)
        documents = experiment_data.get("documents", [])
        for doc in documents:
            if doc.get("storage") == "s3" and doc.get("file_path"):
                s3_key = self.store_document_s3(
                    doc["name"],
                    doc["file_path"],
                    {"experiment_id": experiment_id}
                )
                curation_report["documents_stored"].append({"name": doc["name"], "s3_key": s3_key})
        
        # 3. Log model version (if model data present)
        if "model_data" in experiment_data:
            model_data = experiment_data["model_data"]
            run_id = self.log_model_version(
                model_name=model_data.get("name", "unnamed_model"),
                model_params=model_data.get("params", {}),
                metrics=model_data.get("metrics", {})
            )
            curation_report["model_run_id"] = run_id
        
        # 4. Store embeddings (if present)
        embeddings = experiment_data.get("embeddings", [])
        for emb in embeddings:
            success = self.store_embedding(
                embedding_id=emb["id"],
                vector=emb["vector"],
                payload=emb.get("metadata", {})
            )
            if success:
                curation_report["embeddings_stored"] += 1
        
        # Log the curation itself
        self.log_provenance(
            event_type="data_curated",
            data=curation_report
        )
        
        return curation_report

    def get_experiment_artifacts(self, experiment_id: str) -> Dict[str, Any]:
        """
        Retrieves all artifacts for an experiment.
        """
        return {
            "experiment_id": experiment_id,
            "provenance_chain": self.get_provenance_chain(experiment_id),
            "documents": "See MongoDB/S3",
            "embeddings": "See Vector DB",
            "model_runs": "See MLflow"
        }
