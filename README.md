# Chemora

Chemora is an advanced AI-powered chemistry assistant designed to help students, researchers, and enthusiasts with various chemistry-related tasks. It leverages Large Language Models (LLMs) like Gemini and OpenAI's GPT to provide intelligent responses and perform complex actions.

## Features

-   **Intelligent Chat**: Ask questions about chemical concepts, reactions, and mechanisms.
-   **Agentic Workflow**: Utilizes specialized agents for tasks such as:
    -   **Retrosynthesis**: Planning synthetic routes for molecules.
    -   **Safety Checks**: Assessing safety risks of chemical procedures.
    -   **Literature Retrieval**: Finding relevant scientific papers.
-   **Dual LLM Support**: Configurable to use either Google's Gemini (Free) or OpenAI's GPT models.
-   **Interactive Frontend**: A user-friendly web interface for interacting with the AI.

## Project Structure

-   **Chemora_Backend**: The backend server built with Python/FastAPI, handling the AI logic and agent coordination.
-   **Chemora_Frontend**: The frontend application (likely React/Next.js based on file structure hints) for the user interface.
-   **docs**: Documentation and architectural diagrams.

## Getting Started

### Backend Setup

1.  Navigate to `Chemora_Backend`.
2.  Install dependencies: `pip install -r requirements.txt`.
3.  Set up your API keys (see `GEMINI_SETUP.md` or `OPENAI_SETUP.md`).
4.  Run the server: `uvicorn main:app --reload`.

### Frontend Setup

1.  Navigate to `Chemora_Frontend`.
2.  Install dependencies: `npm install`.
3.  Run the development server: `npm run dev`.

## Contributing

Contributions are welcome! Please create a new branch for your features and submit a pull request.
