"""Stage-level workflow for the existing Satellite ASO maker analysis."""
import glob
import os
import re

configfile: 'workflow/config.yaml'

RUN_ID = config['run_id']
if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', RUN_ID):
    raise ValueError('run_id may contain only letters, digits, _, . and -')
RUN = 'runs/' + RUN_ID
PY = config['python']
T = config['threads']
for name in ('step0', 'step1_discovery', 'step1_counts', 'step2', 'step3_mismatch', 'step4_mismatch'):
    if not isinstance(T[name], int) or T[name] < 1:
        raise ValueError('threads.' + name + ' must be a positive integer')

def done(stage):
    return RUN + '/.workflow/' + stage + '.done.json'

rule all:
    input:
        done('landscape_finish'),
        RUN + '/step4/all72_landscape/COMPLETE.json'

rule bootstrap:
    input:
        'workflow/config.yaml',
        'workflow/bootstrap.py',
        'workflow/stage.py',
        *sorted(p for p in glob.glob('scripts/*') if os.path.isfile(p) and not p.endswith('.pyc'))
    output:
        RUN + '/.bootstrap.json'
    shell:
        '{PY} workflow/bootstrap.py --config workflow/config.yaml --run {RUN}'

rule step0:
    input: rules.bootstrap.output
    output: done('step0')
    threads: T['step0']
    shell: '{PY} workflow/stage.py step0 --run {RUN} --threads {threads}'

rule step1_discovery:
    input: rules.step0.output
    output: done('step1_discovery')
    threads: T['step1_discovery']
    shell: '{PY} workflow/stage.py step1_discovery --run {RUN} --threads {threads}'

rule step1_counts:
    input: rules.step1_discovery.output
    output: done('step1_counts')
    threads: T['step1_counts']
    shell: '{PY} workflow/stage.py step1_counts --run {RUN} --threads {threads}'

rule step1_finish:
    input: rules.step1_counts.output
    output: done('step1_finish')
    shell: '{PY} workflow/stage.py step1_finish --run {RUN}'

rule sweep:
    input: rules.step1_finish.output
    output: done('sweep')
    shell: '{PY} workflow/stage.py sweep --run {RUN}'

rule step2_counts:
    input: rules.sweep.output
    output: done('step2_counts')
    threads: T['step2']
    shell: '{PY} workflow/stage.py step2_counts --run {RUN} --threads {threads}'

rule step2_finish:
    input: rules.step2_counts.output
    output: done('step2_finish')
    shell: '{PY} workflow/stage.py step2_finish --run {RUN}'

rule step3_prepare:
    input: rules.step2_finish.output
    output: done('step3_prepare')
    shell: '{PY} workflow/stage.py step3_prepare --run {RUN}'

rule step3_run:
    input: rules.step3_prepare.output
    output: done('step3_run')
    threads: T['step3_mismatch']
    shell: '{PY} workflow/stage.py step3_run --run {RUN} --threads {threads}'

rule step4_rank:
    input: rules.step3_run.output
    output: done('step4_rank')
    shell: '{PY} workflow/stage.py step4_rank --run {RUN}'

rule step4_neighbors:
    input: rules.step4_rank.output
    output: done('step4_neighbors')
    shell: '{PY} workflow/stage.py step4_neighbors --run {RUN}'

rule step4_mismatch:
    input: rules.step4_neighbors.output
    output: done('step4_mismatch')
    threads: T['step4_mismatch']
    shell: '{PY} workflow/stage.py step4_mismatch --run {RUN} --threads {threads}'

rule step4_finish:
    input: rules.step4_mismatch.output
    output: done('step4_finish')
    shell: '{PY} workflow/stage.py step4_finish --run {RUN}'

rule top3:
    input: rules.step4_finish.output
    output: done('top3')
    shell: '{PY} workflow/stage.py top3 --run {RUN}'

rule scan72:
    input: rules.step4_finish.output
    output: done('scan72')
    threads: 1
    shell: '{PY} workflow/stage.py scan72 --run {RUN}'

rule landscape:
    input: rules.scan72.output
    output: done('landscape')
    shell: '{PY} workflow/stage.py landscape --run {RUN}'

rule landscape_plots:
    input: rules.landscape.output
    output: done('landscape_plots')
    shell: '{PY} workflow/stage.py landscape_plots --run {RUN}'

rule landscape_finish:
    input: rules.landscape_plots.output, rules.top3.output
    output: done('landscape_finish'), RUN + '/step4/all72_landscape/COMPLETE.json'
    shell: '{PY} workflow/stage.py landscape_finish --run {RUN}'
