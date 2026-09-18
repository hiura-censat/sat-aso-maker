#!/usr/bin/env python3
"""Run fast workflow contracts before any long-running stage."""
import argparse,json,subprocess
from pathlib import Path

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',required=True);parser.add_argument('--python',required=True);args=parser.parse_args()
    root=Path(__file__).resolve().parents[1];run=(root/args.run).resolve()
    subprocess.run([args.python,'-m','unittest','discover','-s','scripts','-p','test_*.py'],cwd=run,check=True)
    marker=run/'.workflow/preflight.done.json';marker.parent.mkdir(exist_ok=True)
    marker.write_text(json.dumps({'stage':'preflight','status':'PASS','checks':['unit_tests','schema_contracts']},indent=2)+'\n')
if __name__=='__main__':main()
