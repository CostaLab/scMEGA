"""Validate one issue with the separate metadata proposal, without merging branches."""
import io
from pathlib import Path
import subprocess
import sys
import tarfile

ref, target = sys.argv[1:3]
target = Path(target)
target.mkdir(parents=True, exist_ok=False)
archive = subprocess.check_output(['git','archive',ref])
with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
    tar.extractall(target, filter='data')
patch = subprocess.check_output(['git','diff','e6ce4f9','ec24fdf','--','DESCRIPTION','NAMESPACE','R/imports.R'])
subprocess.run(['git','apply','--check','-'],input=patch,cwd=target,check=True)
subprocess.run(['git','apply','-'],input=patch,cwd=target,check=True)
print('Prepared scratch validation copy of',ref,'plus metadata ec24fdf at',target)
