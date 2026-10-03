"""Bundle the pinned Python plotting runtime into the Pages deployment."""
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import urlopen

BASE = 'https://cdn.jsdelivr.net/pyodide/v314.0.7/full/'
CORE = {
    'pyodide.js': '3141b814715a72e59b51b1b18b9ceae5bf19f7c852417e431bb0a34feadf825c',
    'pyodide.asm.mjs': 'f7cdc8ece80678ceb712f8e65ebe6d3a83203a180c399865f49612a051693635',
    'pyodide.asm.wasm': 'cc36e3cab04fdfc9a63ff13eb52eae2b911bf46c025cc7b281f394bd3de1d5e6',
    'python_stdlib.zip': 'fa1957e5777068fc4f7437f96d860ae2fbe9c19732ba06c84e004ec16dd7dd7a',
    'pyodide-lock.json': '5dc2fc119108bc148c7457dc86e7675b5c87e1cafd420b9c34c1eaef7b36c010',
}
SEABORN = 'https://files.pythonhosted.org/packages/83/11/00d3c3dfc25ad54e731d91449895a79e4bf2384dc3ac01809010ba88f6d5/seaborn-0.13.2-py3-none-any.whl'

def main(destination):
    root = Path(destination) / 'runtime'
    root.mkdir(parents=True, exist_ok=True)
    def fetch(name, sha, url=None):
        target = root / name
        if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest() == sha:
            return
        with urlopen(url or BASE + name, timeout=120) as response:
            data = response.read()
        if hashlib.sha256(data).hexdigest() != sha:
            raise RuntimeError('Checksum mismatch: ' + name)
        target.write_bytes(data)
        print('Bundled', name, len(data), flush=True)
    for name, sha in CORE.items():
        fetch(name, sha)
    packages = json.loads((root / 'pyodide-lock.json').read_text())['packages']
    visited = set()
    def package(name):
        if name in visited:
            return
        visited.add(name)
        item = packages[name]
        for dependency in item['depends']:
            package(dependency)
        fetch(item['file_name'], item['sha256'])
    for name in ['numpy', 'pandas', 'matplotlib', 'micropip']:
        package(name)
    fetch(SEABORN.rsplit('/', 1)[-1], '636f8336facf092165e27924f223d3c62ca560b1f2bb5dff7ab7fad265361987', SEABORN)
    print('Runtime ready:', len(visited), 'packages plus seaborn')

if __name__ == '__main__':
    main(sys.argv[1])
