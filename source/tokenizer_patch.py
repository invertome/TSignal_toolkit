# Patch to fix tokenizer compatibility issue
# The saved model has a tokenizer from transformers 3.3.1
# which is missing attributes like tokens_trie needed by newer versions

import torch
from transformers import BertTokenizer

def patch_model_tokenizer(model):
    """
    Replace the deserialized tokenizer with a fresh one to ensure compatibility.
    """
    # Get the original tokenizer's vocab/config
    tokenizer_name = getattr(model.hparams, 'tokenizer_name', 'Rostlab/prot_bert_bfd')

    # Create a fresh tokenizer with the same config
    try:
        new_tokenizer = BertTokenizer.from_pretrained(tokenizer_name, do_lower_case=False)
        model.tokenizer = new_tokenizer
        print(f"Patched tokenizer: reinitialized from {tokenizer_name}")
    except Exception as e:
        print(f"Warning: Could not reinitialize tokenizer: {e}")

    return model

def load_model_with_patched_tokenizer(model_path):
    """
    Load model and patch the tokenizer for compatibility.
    """
    model = torch.load(model_path, map_location='cpu')
    model = patch_model_tokenizer(model)
    return model
