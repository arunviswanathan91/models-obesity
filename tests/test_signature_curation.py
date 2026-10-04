"""Offline checks for SDK configuration and human-reviewed signature changes."""
import contextlib
import io
import os
from pathlib import Path
import runpy
import sys
import types
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / 'analysis/signatures/curate_signatures.py'


class SignatureCurationTests(unittest.TestCase):
    def setUp(self):
        google = types.ModuleType('google')
        sdk = types.ModuleType('google.generativeai')
        api_types = types.ModuleType('google.generativeai.types')
        api_types.HarmCategory = types.SimpleNamespace(**{name: name for name in (
            'HARM_CATEGORY_HATE_SPEECH', 'HARM_CATEGORY_HARASSMENT',
            'HARM_CATEGORY_SEXUALLY_EXPLICIT', 'HARM_CATEGORY_DANGEROUS_CONTENT')})
        api_types.HarmBlockThreshold = types.SimpleNamespace(BLOCK_NONE='BLOCK_NONE')
        self.calls = []
        calls = self.calls

        class Model:
            def __init__(self, name):
                calls.append(('model', name))

            def generate_content(self, prompt, *, generation_config, safety_settings):
                calls.append(('request', generation_config.copy()))
                return types.SimpleNamespace(text='```json\n{"choice": 3, "reasoning": "distinct"}\n```')

        sdk.configure = lambda **kwargs: calls.append(('configure', kwargs))
        sdk.GenerativeModel = Model
        sdk.types = api_types
        google.generativeai = sdk
        with patch.dict(sys.modules, {'google': google, 'google.generativeai': sdk,
                                      'google.generativeai.types': api_types}), \
             patch.dict(os.environ, {'GEMINI_API_KEY': 'offline-test-value'}):
            self.module = runpy.run_path(str(SCRIPT), run_name='offline_test')

    def test_sdk_key_and_generation_configuration(self):
        result = self.module['call_gemini_with_retry']('test prompt')
        self.assertEqual(result['choice'], 3)
        self.assertIn(('configure', {'api_key': 'offline-test-value'}), self.calls)
        self.assertIn(('request', {'response_mime_type': 'application/json', 'temperature': 0.2}), self.calls)

    def test_addition_overlap_boundary_and_human_acceptance(self):
        run = self.module['run_discovery']
        namespace = run.__globals__
        namespace['suggest_missing_signatures'] = lambda *_: [
            {'name': 'rejected_overlap', 'genes': ['A', 'B', 'C', 'X', 'Y']},
            {'name': 'accepted_boundary', 'genes': ['A', 'V', 'W', 'X', 'Y']},
            {'name': 'human_rejected', 'genes': ['I', 'J', 'K', 'L', 'M']}]
        data = {'Cell': {'existing': ['A', 'B', 'C', 'D', 'E']}}
        with patch('builtins.input', side_effect=['y', 'n']), contextlib.redirect_stdout(io.StringIO()):
            run(data, io.StringIO())
        self.assertEqual(set(data['Cell']), {'existing', 'accepted_boundary'})

    def test_merge_uses_reviewers_choice(self):
        run = self.module['run_deduplication']
        namespace = run.__globals__
        namespace['ask_gemini_deduplication'] = lambda _: {'choice': 1, 'reasoning': 'test'}
        namespace['get_manual_dedup_choice'] = lambda _: 4
        data = {'Cell': {'first': ['A', 'B', 'C'], 'second': ['A', 'B', 'D']}}
        with contextlib.redirect_stdout(io.StringIO()):
            run(data, io.StringIO())
        self.assertEqual(data, {'Cell': {'first': ['A', 'B', 'C', 'D']}})


if __name__ == '__main__':
    unittest.main()
