import sys
import subprocess

with open('fasta/prochlorococcus.txt') as f:
    r = f.readlines()
    r = [i.strip() for i in r]

for i in r:
    subprocess.call(
        [
            "prokka",
            "--outdir",
            f"prokka_runs/{i[:-6]}",
            "--genus",
            "Prochloroccus",
            f"fasta/prochlorococcus/{i}"
        ]
    )
