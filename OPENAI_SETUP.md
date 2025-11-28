# Setting Up ChatGPT-Style Chemistry AI

## OpenAI API Configuration

To enable the full conversational AI capabilities, you need to set up your OpenAI API key.

### Option 1: Environment Variable (Recommended)
```bash
export OPENAI_API_KEY="sk-your-api-key-here"
```

Add this to your `~/.zshrc` or `~/.bashrc` to make it permanent.

### Option 2: In Python Code
```python
import openai
openai.api_key = "sk-your-api-key-here"
```

### Option 3: .env File
Create a `.env` file in the backend directory:
```
OPENAI_API_KEY=sk-your-api-key-here
```

Then load it in `main.py`:
```python
from dotenv import load_dotenv
load_dotenv()
```

## Getting an API Key

1. Go to https://platform.openai.com/
2. Sign up or log in
3. Navigate to API Keys section
4. Click "Create new secret key"
5. Copy the key (starts with `sk-`)
6. Set it as environment variable

## Cost Considerations

**GPT-4o-mini** (currently used):
- ~$0.15 per 1M input tokens
- ~$0.60 per 1M output tokens
- Very affordable for chat usage

**GPT-4** (optional upgrade):
- Better quality responses
- More expensive (~10x cost)
- Change model in chat_agent.py line 561

## Testing

```bash
# Test if API key is set
python3 -c "import openai; print('API Key:', openai.api_key[:10] + '...' if openai.api_key else 'Not set')"

# Start backend
uvicorn main:app --reload --port 8000
```

## Usage Examples

Once configured, the chat will answer ANY chemistry question:

**Conceptual Questions:**
```
USER: "Why is benzene aromatic?"
AI: [Explains Hückel's rule, resonance, stability...]
```

**Mechanisms:**
```  
USER: "Explain the mechanism of E2 elimination"
AI: [Shows concerted mechanism, anti-periplanar geometry...]
```

**Calculations:**
```
USER: "Calculate pH of 0.1M acetic acid (pKa = 4.76)"
AI: [Shows Henderson-Hasselbalch calculation...]
```

**Troubleshooting:**
```
USER: "My recrystallization keeps oiling out, help!"
AI: [Suggests solutions: seed crystals, scratch, cooling rate...]
```

**Still Uses Agents for:**
- Synthesis planning
- Safety assessments  
- Literature searches
- Protocol generation

## Fallback Mode

Without API key, the chat still works but only for agent-powered tasks:
- Synthesis: ✅ (Uses RetrosynthesisAgent)
- Safety: ✅ (Uses SafetyCheckerAgent)
- Literature: ✅ (Uses LiteratureRetrievalAgent)
- General Q&A: ❌ (Needs OpenAI)

---

**With OpenAI configured, Chemora becomes "ChatGPT for Chemistry"!** 🧪💬
