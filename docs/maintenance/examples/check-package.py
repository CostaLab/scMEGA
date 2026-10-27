"""Run actual install/build/check and retain every outcome, including failure."""
import json
import os
from pathlib import Path
import subprocess
import sys

root = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
out.mkdir(parents=True, exist_ok=True)
lib = out / 'library'
lib.mkdir(exist_ok=True)
env = dict(os.environ, R_LIBS_USER='/tmp/scmega-r-lib')
env.pop('SCMEGA_SOURCE_TESTS', None)
env.pop('SCMEGA_SOURCE_ROOT', None)
commands = [('install', ['R', 'CMD', 'INSTALL', '--library='+str(lib), str(root)]),
            ('build', ['R', 'CMD', 'build', '--no-build-vignettes', '--no-manual', str(root)])]
results = {}
for name, command in commands:
    with (out / (name+'.log')).open('w') as log:
        p = subprocess.run(command, cwd=out, env=env, stdout=log, stderr=subprocess.STDOUT)
    results[name] = {'command': command, 'exit_code': p.returncode}
if results['build']['exit_code'] == 0:
    tarballs = list(out.glob('scMEGA_*.tar.gz'))
    assert len(tarballs) == 1
    command = ['R', 'CMD', 'check', '--no-manual', '--no-build-vignettes', str(tarballs[0])]
    with (out / 'check.log').open('w') as log:
        p = subprocess.run(command, cwd=out, env=env, stdout=log, stderr=subprocess.STDOUT)
    results['check'] = {'command': command, 'exit_code': p.returncode}
(out / 'results.json').write_text(json.dumps(results, indent=2)+'\n')
print(json.dumps(results))
