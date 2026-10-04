#!/usr/bin/env python3
"""Interactive Gemini-assisted review of PDAC and obesity gene signatures."""
import json
import time
import os
import re
import sys

try:
    import google.generativeai as genai
    from google.generativeai.types import HarmCategory, HarmBlockThreshold
except ImportError:
    print(" Critical Error: 'google-generativeai' library not found.")
    print("   Please install it: pip install -q -U google-generativeai")
    exit(1)

GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
if not GEMINI_API_KEY:
    raise ValueError("Set GEMINI_API_KEY as environment variable")
genai.configure(api_key=GEMINI_API_KEY)

JSON_PATH = os.getenv("SIGNATURE_INPUT", "/content/drive/MyDrive/new_dysfunction/signature/Corrected_Signatures_02112025.json")
OUTPUT_PATH = os.getenv("SIGNATURE_OUTPUT", "/content/drive/MyDrive/new_dysfunction/signature/Final_Refined_Signatures_Gemini.json")
LOG_PATH = os.getenv("SIGNATURE_LOG", "/content/drive/MyDrive/new_dysfunction/signature/Interactive_Log.txt")

DEDUPLICATION_OVERLAP_THRESHOLD = 0.50
NEW_SIG_OVERLAP_THRESHOLD = 0.20
MIN_GENES = 8
MAX_GENES = 12

ACTIVE_MODEL_NAME = "models/gemini-2.5-flash"

def clean_json_response(text):
    """Removes markdown wrappers"""
    text = text.strip()
    text = re.sub(r'^```(json)?\s*', '', text)
    text = re.sub(r'\s*```$', '', text)
    return text

def overlap_coefficient(set1, set2):
    if not set1 or not set2: return 0.0
    intersection = len(set1 & set2)
    min_size = min(len(set1), len(set2))
    return intersection / min_size

def print_gene_breakdown(shared, unique1, unique2, name1, name2):
    """Visualizes the gene overlap"""
    print("\nGENE BREAKDOWN:")

    s_list = sorted(list(shared))
    s_str = ", ".join(s_list[:15]) + ("..." if len(s_list) > 15 else "")
    print(f"SHARED ({len(s_list)}): {s_str}")

    u1_list = sorted(list(unique1))
    if u1_list:
        u1_str = ", ".join(u1_list[:15]) + ("..." if len(u1_list) > 15 else "")
        print(f"UNIQUE to '{name1}' ({len(u1_list)}): {u1_str}")
    else:
        print(f"UNIQUE to '{name1}': (none) - SUBSET")

    u2_list = sorted(list(unique2))
    if u2_list:
        u2_str = ", ".join(u2_list[:15]) + ("..." if len(u2_list) > 15 else "")
        print(f"UNIQUE to '{name2}' ({len(u2_list)}): {u2_str}")
    else:
        print(f"UNIQUE to '{name2}': (none) - SUBSET")

def get_manual_dedup_choice(ai_recommendation):
    """Interactive menu for deduplication"""
    print("\n" + "-"*80)
    print(f"AI RECOMMENDATION: Choice {ai_recommendation['choice']}")
    print(f"   Reason: {ai_recommendation['reasoning']}")
    print("-" * 80)
    print("   1. Keep Sig 1 ONLY (Delete Sig 2)")
    print("   2. Keep Sig 2 ONLY (Delete Sig 1)")
    print("   3. Keep BOTH")
    print("   4. MERGE into Sig 1 Name")
    print("   5. MERGE into Sig 2 Name")
    print("   S. Skip")
    print("-" * 80)
    
    while True:
        choice = input("👉 Your Decision (1-5, S): ").strip().upper()
        if choice in ['1', '2', '3', '4', '5']: return int(choice)
        if choice == 'S': return None

generation_config={
    "response_mime_type": "application/json",
    "temperature": 0.2
}

def call_gemini_with_retry(prompt, retries=3):
    """Calls Gemini with retry logic"""
    if not ACTIVE_MODEL_NAME: return None
    model = genai.GenerativeModel(ACTIVE_MODEL_NAME)
    
    safety_settings = {
        HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
    }

    for attempt in range(retries):
        try:
            response = model.generate_content(
                prompt, 
                generation_config=generation_config,
                safety_settings=safety_settings
            )
            cleaned_text = clean_json_response(response.text)
            return json.loads(cleaned_text)

        except Exception as e:
            error_str = str(e)
            if "429" in error_str:
                time.sleep((attempt + 1) * 2)
            elif "404" in error_str:
                print(f"Model 404 Error: {ACTIVE_MODEL_NAME}")
                return None
            else:
                time.sleep(1)
    
    return None

def ask_gemini_deduplication(pair_info):
    prompt = f"""
    Resolve duplicate gene signatures.
    Context: Pancreatic Cancer (PDAC) + Obesity.
    
    Sig 1: "{pair_info['sig1']}"
    Sig 2: "{pair_info['sig2']}"
    Overlap: {pair_info['overlap_coef']:.1%}
    
    Respond JSON only: {{ "choice": 1-5, "reasoning": "sentence" }}
    Choices: 1=Keep 1, 2=Keep 2, 3=Keep Both, 4=Merge into 1, 5=Merge into 2.
    """
    return call_gemini_with_retry(prompt)

def suggest_missing_signatures(cell_type, current_signatures):
    existing_names = ", ".join(list(current_signatures.keys()))
    prompt = f"""
    You are an expert in Pancreatic Adenocarcinoma (PDAC) and Obesity.
    Analyze the cell type '{cell_type}'.
    
    EXISTING SIGNATURES: {existing_names}
    
    TASK: Suggest up to 3 MISSING gene signatures that specifically drive OBESITY-INDUCED DYSFUNCTION in PDAC for this cell type.
    
    CRITERIA:
    1. Must NOT overlap >20% with existing signatures.
    2. Must contain {MIN_GENES}-{MAX_GENES} genes.
    3. Must be biologically critical for obesity-driven cancer progression.
    
    Respond JSON only:
    {{ 
      "suggestions": [
        {{ "name": "Standard Name", "genes": ["G1", "G2"...], "rationale": "Why this specific signature matters for obesity in PDAC" }}
      ]
    }}
    Or empty list if none.
    """
    result = call_gemini_with_retry(prompt)
    if result: return result.get('suggestions', [])
    return []

def run_deduplication(signatures, log_file):
    print(f"\n{'='*80}\n⚙️  PHASE 1: DEDUPLICATION REVIEW\n{'='*80}")
    
    pairs = []
    for cell, sigs in signatures.items():
        names = list(sigs.keys())
        for i, s1 in enumerate(names):
            for s2 in names[i+1:]:
                g1 = set(str(g).upper() for g in sigs[s1])
                g2 = set(str(g).upper() for g in sigs[s2])
                ov = overlap_coefficient(g1, g2)
                if ov > DEDUPLICATION_OVERLAP_THRESHOLD:
                    pairs.append({
                        'cell_type': cell, 'sig1': s1, 'sig2': s2,
                        'genes1': g1, 'genes2': g2, 'overlap_coef': ov,
                        'shared': g1 & g2
                    })

    if not pairs:
        print("No duplicates found.")
        return

    for i, pair in enumerate(pairs, 1):
        cell = pair['cell_type']
        s1, s2 = pair['sig1'], pair['sig2']
        
        if s1 not in signatures[cell] or s2 not in signatures[cell]: continue

        print(f"\n[{i}/{len(pairs)}] 📍 {cell}: '{s1}' ⟷ '{s2}'")
        print(f"   Overlap: {pair['overlap_coef']:.1%}")

        unique1 = pair['genes1'] - pair['genes2']
        unique2 = pair['genes2'] - pair['genes1']
        print_gene_breakdown(pair['shared'], unique1, unique2, s1, s2)

        ai_rec = ask_gemini_deduplication(pair)
        if not ai_rec:
            print("AI Error. Skipping.")
            continue
            
        choice = get_manual_dedup_choice(ai_rec)

        if choice == 1:
            del signatures[cell][s2]
            log_file.write(f"DEDUP {cell}: Kept {s1}, Deleted {s2}\n")
            print("   Action: Deleted Signature 2")
        elif choice == 2:
            del signatures[cell][s1]
            log_file.write(f"DEDUP {cell}: Kept {s2}, Deleted {s1}\n")
            print(" Action: Deleted Signature 1")
        elif choice == 3:
            log_file.write(f"DEDUP {cell}: Kept Both\n")
            print("   Action: Kept Both")
        elif choice == 4:
            combined = sorted(list(pair['genes1'] | pair['genes2']))
            signatures[cell][s1] = combined
            del signatures[cell][s2]
            log_file.write(f"DEDUP {cell}: Merged into {s1}\n")
            print("   Action: Merged into Sig 1")
        elif choice == 5:
            combined = sorted(list(pair['genes1'] | pair['genes2']))
            signatures[cell][s2] = combined
            del signatures[cell][s1]
            log_file.write(f"DEDUP {cell}: Merged into {s2}\n")
            print("   Action: Merged into Sig 2")

def run_discovery(signatures, log_file):
    print(f"\n{'='*80}\n PHASE 2: MISSING SIGNATURE DISCOVERY\n{'='*80}")
    print("   Searching for: Obesity-Induced Dysfunction in Pancreatic Cancer")
    
    for idx, cell_type in enumerate(signatures.keys(), 1):
        print(f"\n Analyzing {cell_type} ({idx}/{len(signatures)})...")
        
        suggestions = suggest_missing_signatures(cell_type, signatures[cell_type])
        
        if not suggestions:
            print("   (No missing signatures detected)")
            continue
            
        for sig in suggestions:
            name = sig['name']
            genes = sorted(sig['genes'])
            rationale = sig.get('rationale', 'No rationale provided')

            if name in signatures[cell_type]: continue
            
            is_valid = True
            current_genes = [set(str(g).upper() for g in gs) for gs in signatures[cell_type].values()]
            new_gene_set = set(str(g).upper() for g in genes)
            
            for existing_set in current_genes:
                if overlap_coefficient(new_gene_set, existing_set) > NEW_SIG_OVERLAP_THRESHOLD:
                    is_valid = False
                    break
            
            if not is_valid: continue 

            print(f"\n💡 SUGGESTION FOR {cell_type}:")
            print(f"   Name: {name}")
            print(f"   Rationale: {rationale}")
            print(f"   🧬 Genes ({len(genes)}): {', '.join(genes)}")
            
            while True:
                choice = input(f"   👉 Add this signature? [y/n]: ").strip().lower()
                if choice in ['y', 'yes', 'n', 'no']: break
            
            if choice in ['y', 'yes']:
                signatures[cell_type][name] = genes
                print(f" Added '{name}'")
                log_file.write(f"ADDED {cell_type}: {name}\n")
            else:
                print(f" Rejected")

if __name__ == "__main__":
    
    print(f" Loading {JSON_PATH}...")
    with open(JSON_PATH, 'r') as f:
        data = json.load(f)
        
    with open(LOG_PATH, 'w') as log:
        log.write(f"Session: {time.ctime()}\nModel: {ACTIVE_MODEL_NAME}\n\n")
        
        try:
            run_deduplication(data, log)
            run_discovery(data, log)
        except KeyboardInterrupt:
            print("\n\n Interrupted by user. Saving current progress...")

    print(f"\n Saving to {OUTPUT_PATH}...")
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(data, f, indent=2)
        
    print(" Done!")

