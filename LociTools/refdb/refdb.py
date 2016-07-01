
import logging
import os.path

from glob import glob

class RefDb(object):

    def ___init__(self, dbPath):
        fastas = []
        for suffix in (".fa", ".fna", ".fasta"):
            fastas.extend(glob(os.path.join(dbPath, "*.{0}".format(suffix))))
        refs = dict()
        for fasta in fastas:
            bn = os.path.basename(fasta)
            loci, _ = os.path.splitext(bn)
            suffixArray = fasta + ".sa"
            if not os.path.exists(suffixArray):
                logging.warn("missing suffix array for : '{0}'".format(fasta))
                suffixArray = None
            refs[loci] = (fasta, suffixArray)
        self.__refs = refs

    def __iter__(self):
        return iter(self.__refs)
