#! /bin/env python

import os
import argparse
import stat

parser = argparse.ArgumentParser(description='Parse the jobid for use in naming the outfiles', add_help=True)
parser.add_argument('outdir', help='output directory')
parser.add_argument('infiles', nargs='+', help='input files')
args = parser.parse_args()

conda_path = os.environ['CONDA_PREFIX'].replace('envs/mlenv', '')

with open('batch.sh', 'w') as batch:
    batch.write('''#! /bin/bash

. {CONDAPATH}/etc/profile.d/conda.sh
export PATH={CONDAPATH}/bin:$PATH

set -o errexit
set -o nounset

SANDBOX=$PWD
conda activate mlenv

{SCRIPTPATH}/preprocess-ntuples.py $@

mv *.tfrecord {OUTDIR}/.
    '''.format(
        SCRIPTPATH=os.path.abspath(os.path.dirname(__file__)),
        OUTDIR=os.path.abspath(args.outdir),
        CONDAPATH=conda_path,
    )
            )

os.chmod(
    'batch.sh', 
    os.stat('batch.sh').st_mode | stat.S_IEXEC
)

with open('condor.jdl', 'w') as jdl:
    jdl.write('''
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
getenv = True
executable = ./batch.sh
+MaxRuntime = 160600
requirements = (OpSysAndVer =?= "CentOS7")
    ''')

    for idx, loc_infile in enumerate(args.infiles):
        infile = os.path.abspath(loc_infile)
        jdl.write('''
Output = {PROCID}.out
Error = {PROCID}.err
Log = {PROCID}.log
Arguments = {INFILE} {OUTFILE}
Queue
        '''.format(PROCID=idx, INFILE=infile, OUTFILE=os.path.basename(infile).replace('.root', '.tfrecord')))
