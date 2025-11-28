# Using Gemini API for Chat (FREE!)

## Why Gemini?

✅ **FREE** - Generous free tier (60 requests/minute!)
✅ **Fast** - Gemini 1.5 Flash is optimized for speed
✅ **Smart** - Excellent chemistry knowledge
✅ **Easy** - No credit card required

## Setup (2 Minutes)

### Step 1: Get Free API Key

1. Go to **https://aistudio.google.com/app/apikey**
2. Sign in with Google account
3. Click **"Create API Key"**
4. Copy the key (starts with `AIza...`)

### Step 2: Install Gemini Library

```bash
cd /Users/rishwantharyansripathi/Desktop/Major\ project/Chemora_Backend
pip install google-generativeai
```

### Step 3: Set Environment Variable

```bash
export GEMINI_API_KEY="AIza_your_key_here"
```

**Make it permanent:**
Add to your `~/.zshrc`:
```bash
echo 'export GEMINI_API_KEY="AIza_your_key_here"' >> ~/.zshrc
source ~/.zshrc
```

### Step 4: Restart Backend

```bash
# Stop current backend (Ctrl+C)
# Then restart:
uvicorn main:app --reload --port 8000
```

You should see: `Chat Agent LLM: gemini` in the console!

## That's It! 🎉

Now your chat can answer **ANY chemistry question**:

### Try These:

**Conceptual:**
```
"Why is benzene aromatic?"
"Explain molecular orbital theory"
```

**Mechanisms:**
```
"Show me the mechanism of aldol condensation"
"How does SN2 reaction work?"
```

**Calculations:**
```
"Calculate pH of 0.1M acetic acid with pKa 4.76"
"What's the theoretical yield if I start with 5g?"
```

**Troubleshooting:**
```
"My recrystallization keeps oiling out, help!"
"Why is my NMR showing extra peaks?"
```

**Synthesis (Hybrid - Uses Both Gemini + Your Agents):**
```
"How do I synthesize aspirin?"
→ Uses RetrosynthesisAgent for routes
→ Uses Gemini to explain the chemistry
```

## Features

| Feature | Free Tier |
|---------|-----------|
| Requests/minute | 60 |
| Requests/day | 1,500 |
| Cost | **$0.00** |
| Model | Gemini 1.5 Flash |
| Context window | 1M tokens |

## Gemini vs OpenAI

| | Gemini 1.5 Flash | GPT-4o-mini |
|---|---|---|
| **Cost** | FREE | $0.15/1M tokens |
| **Speed** | Very Fast | Fast |
| **Chemistry** | Excellent | Excellent |
| **Free Tier** | YES (generous) | $5 credit (limited) |
| **Best For** | Development, students | Production at scale |

## Verification

Test if it's working:

```bash
# In terminal
python3 << EOF
import google.generativeai as genai
import os
genai.configure(api_key=os.getenv("GEMINI_API_KEY"))
model = genai.GenerativeModel('gemini-1.5-flash')
response = model.generate_content("What is benzene?")
print("Success! Gemini is working:")
print(response.text[:100] + "...")
EOF
```

## Troubleshooting

**"Module not found"**
```bash
pip install google-generativeai
```

**"API key not valid"**
- Double-check the key from https://aistudio.google.com/app/apikey
- Make sure it starts with `AIza`
- Verify it's set: `echo $GEMINI_API_KEY`

**"Quota exceeded"**
- Free tier: 60 requests/min, 1500/day
- Wait a minute or upgrade to paid (still very cheap)

## Upgrading to Gemini Pro

For even better responses, change in `chat_agent.py` line 619:
```python
model_name='gemini-1.5-pro',  # More capable, still free tier!
```

---

**Recommended: Use Gemini for development (free) and consider OpenAI for production if needed.**

Enjoy your free ChatGPT-for-Chemistry! 🧪✨
