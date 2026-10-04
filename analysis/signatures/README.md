# Interactive signature curation

[curate_signatures.py](curate_signatures.py) provides human-reviewed Gemini assistance for a JSON catalogue of cell-type-specific gene signatures in the PDAC and obesity context. It starts from an existing catalogue. It does not create the single-cell expression reference used by BayesPrism.

The input schema is `{cell_type: {signature_name: [gene_symbols]}}`. The script reviews redundancy and then proposes additions. Retention, deletion, merging and additions require a console decision; suggestions alone do not alter the catalogue.

## Run

Install [requirements.txt](requirements.txt), provide `GEMINI_API_KEY` in the environment, and configure paths:

```bash
python -m pip install -r analysis/signatures/requirements.txt
export SIGNATURE_INPUT=/path/to/input_signatures.json
export SIGNATURE_OUTPUT=/path/to/new_refined_signatures.json
export SIGNATURE_LOG=/path/to/new_curation_log.txt
python analysis/signatures/curate_signatures.py
```

The output and log directories must already exist. Use new output and log filenames: the script writes both in overwrite mode. The default model remains `models/gemini-2.5-flash`, with JSON responses and temperature 0.2. The retained legacy SDK and model require service availability. No credentials are stored in the repository.

## Implemented rules

- Redundancy uses the overlap coefficient, `|A ∩ B| / min(|A|, |B|)`, on uppercased gene symbols within each cell type. Pairs exceeding 0.50 enter review.
- The deduplication prompt passes signature names and overlap percentage. Gene lists are displayed to the human reviewer, but are not included in that prompt.
- The addition prompt requests up to three signatures containing 8–12 genes. The script rejects existing names and proposals with an overlap coefficient greater than 0.20 against any current signature in that cell type; exactly 0.20 passes. The 8–12 gene limit is a prompt instruction, not an implemented length check.
- Candidate duplicate pairs are computed once before review. The script preserves that original candidate-list behaviour, including during successive merges.
- The log records applied actions and the model name. It does not retain complete prompts, model responses or independent literature validation. Biological suitability and supporting evidence require human review.

## Provenance

This workflow comes from [obese-model at fe14ec8](https://github.com/arunviswanathan91/obese-model/tree/fe14ec848bfeac469528dc87c5aea8aacb14a87c/03_LLM_signature_curation_with). Prompts, overlap rules and human-decision logic are retained. Two SDK integration issues were repaired: configuring the supplied API key and passing temperature inside `generation_config`, as required by the [SDK interface](https://github.com/google-gemini/deprecated-generative-ai-python/blob/main/docs/api/google/generativeai/GenerativeModel.md). Paths accept environment overrides. See [provenance.json](provenance.json).

The script documents the repository-hosted curation process. The frozen catalogue used for the manuscript must be obtained separately; rerunning an interactive LLM workflow is not a deterministic reconstruction of it. No live LLM requests were made during code verification.
