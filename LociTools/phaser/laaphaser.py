
from __future__ import print_function

import csv
import logging
import os
import os.path

from collections import namedtuple
from subprocess import Popen, PIPE
from shutil import rmtree
from tempfile import mkdtemp

from pbcore.io import FastqReader

LaaRecord = namedtuple("LaaRecord", ["id", "sequence", "qualities", "quality", "barcode", "cluster", "phase", "coverage", "noise", "chimera", "converged", "subreads"])

def which(exe, env=os.environ):
    def isExe(p):
        return os.path.isfile(p) and os.access(p, os.X_OK)
    if os.path.split(exe)[0] and isExe(exe):
        return exe
    for d in env.get("PATH", "").split(os.pathsep):
        exePath = os.path.join(d.strip("\""), exe)
        if isExe(exePath):
            return exePath
    return None

class LaaPhaser(object):

    def __init__(self, barcode, dataset, options=""):
        opts = options.split(' ') if options else []
        illegalOpts = set(["--doBc", "--resultFile", "--reportsFile", "--subreadsReportPrefix"])
        invalidOpts = set(opts) & illegalOpts
        if invalidOpts:
            raise RuntimeError("invalid options to laa: '{0}'".format(" ".join(invalidOpts)))
        self.__laaOpts = opts
        self.__barcode = barcode
        self.__dataset = dataset
        self.__records = None

    def __enter__(self):
        tmpdir = mkdtemp()
        try:
            laa = which("laa")
            if not laa:
                raise RuntimeError("laa not on PATH")
            cmd = [laa, "--doBc", self.__barcode]
            cmd.extend(self.__laaOpts)
            cmd.append(self.__dataset)
            logging.info("running `{0}` in '{1}'".format(" ".join(cmd), tmpdir))
            proc = Popen(cmd, cwd=tmpdir, stderr=PIPE, close_fds=True)
            proc.wait()
            if proc.returncode != 0:
                raise RuntimeError("`{0}` failed with exit code {1}:\n{2}".format(' '.join(cmd), proc.returncode, proc.stderr.read()))
            with open(os.path.join(tmpdir, "amplicon_analysis_summary.csv")) as summ:
                rdr = csv.reader(summ)
                hdr = next(rdr)
                cluster = hdr.index("CoarseCluster") - 2
                phase = hdr.index("Phase") - 2
                coverage = hdr.index("TotalCoverage") - 2
                readQuality = hdr.index("PredictedAccuracy") - 2
                didConverge = hdr.index("ConsensusConverged") - 2
                isChimera = hdr.index("IsChimera") - 2
                recData = dict((row[1], row[2:]) for row in rdr)
            with open(os.path.join(tmpdir, "amplicon_analysis_subreads.{0}.csv".format(self.__barcode))) as rds:
                rdr = csv.reader(rds)
                hdr = next(rdr)
                srData = dict((name, dict()) for name in hdr[1:])
                for row in rdr:
                    for i in range(1, len(row)):
                        weight = float(row[i])
                        if weight > 0.0:
                            srData[hdr[i]][row[0]] = weight
            records = []
            for fname, isNoise in (("amplicon_analysis.fastq", False), ("amplicon_analysis_chimeras_noise.fastq", True)):
                rdr = FastqReader(os.path.join(tmpdir, fname))
                for rec in rdr:
                    attrs = recData[rec.id]
                    subreads = srData[rec.id]
                    records.append(LaaRecord(rec.id, rec.sequence, rec.quality, float(attrs[readQuality]),
                                             self.__barcode, int(attrs[cluster]), int(attrs[phase]), int(attrs[coverage]),
                                             isNoise, bool(int(attrs[isChimera])), bool(int(attrs[didConverge])), subreads))
            self.__records = records
        finally:
            rmtree(tmpdir)
        return self

    def __exit__(self, typ, val, traceback):
        self.__records = None

    def __iter__(self):
        if self.__records is None:
            raise RuntimeError("LaaPhaser is a context object! Use it as such to generate records")
        return iter(self.__records)

# for testing purposes
if __name__ == "__main__":
    import sys
    bc, ds = sys.argv[1:]
    logging.basicConfig(level=logging.INFO)
    with LaaPhaser(bc, os.path.abspath(ds)) as lp:
        for r in lp:
            print(repr(r))
    sys.exit(0)
