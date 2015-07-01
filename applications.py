from params import *
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
import os, sys
import collections
import numpy as np
import subprocess, threading

class TerminationPipe(object):
    #Background process class
    def __init__(self, cmd, timeout, silent=True):
        self.cmd = cmd
        self.timeout = timeout
        self.process = None
        self.output = None
        self.failure = False
        self.stderr = 'EMPTY'
        self.stdout = 'EMPTY'
        self.silent = silent

    def run(self, silent=None, changeDir=False):
        def silentTarget():
            if sys.platform == 'win32':
                if changeDir:
                    self.process = subprocess.Popen("requires\\" + self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
                else:
                    self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE, shell=True, stderr=subprocess.PIPE)
            else:
                if changeDir:
                    self.process = subprocess.Popen("./requires/" + self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
                else:
                    self.process = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,shell=True, stderr=subprocess.PIPE)
            self.output = self.process.communicate()

        def loudTarget():
            if sys.platform == 'win32':
                if changeDir:
                    self.process = subprocess.Popen("requires\\" + self.cmd, shell=False)
                else:
                    self.process = subprocess.Popen(self.cmd, shell=False)
            else:
                if changeDir:
                    self.process = subprocess.Popen("./requires/" + self.cmd, shell=False)
                else:
                    self.process = subprocess.Popen(self.cmd, shell=False)
            self.output=self.process.communicate()

        if silent: self.silent = silent
        if self.silent:
            thread = threading.Thread(target=silentTarget)
        else:
            thread = threading.Thread(target=loudTarget)
        thread.start()
        thread.join(self.timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
            self.failure = True


def executeCMD(cmd, prog, silent=True, useos=False):
    """Execute a command line in the shell"""

    if useos:
        os.system(cmd)
        return None
    else :
        p = subprocess.Popen(
            cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        out, err = p.communicate()
        if(not silent):
            print "%s : \n----------------- %s" % (prog, out)
        return err


def call_seqprog(commandLine, outputFile, timeout, outputFormat='fasta', silent=True):
    """clean call to alignment or alignment edition"""
    geneOutput = None
    retry = True
    aligned = False
    while retry :
        pipe = TerminationPipe(commandLine, timeout)
        pipe.run(silent=silent)
        if not pipe.failure:
            try:
                geneOutput = AlignIO.read(outputFile, outputFormat)
            except Exception as e:
                print e
                raise RuntimeError(commandLine+"\nunable to run:\n\tHave you write access?\tCheck you gene selection")
            aligned = True
        else:
            print("Action not completed in time allowed (%d)"%timeout)
            t = raw_input("Retry with new timeout (enter to skip gene) ? ")
            t = t.strip()
            if(t and isInt(t)):
                timeout = int(t)
            else:
                retry = False

        return aligned, geneOutput


def alignstat(aligndict, allist, gapType='-'):
    output = collections.defaultdict(dict)
    for key in allist:
        if key in aligndict.keys():
            align = aligndict[key]
            gapLength = [align.get_alignment_length()]*len(align)
            unGapLength = [len(x.seq.ungap(gapType)) for x in align]
            gapNumber = []
            for g, u in zip(gapLength, unGapLength):
                gapNumber.append(g-u)
            output[key]['gapmean'] = np.mean(gapNumber)
            output[key]['gapmedian'] = np.median(gapNumber)
            output[key]['gapstd'] = np.std(gapNumber)
            output[key]['gapmax'] = max(gapNumber)
            output[key]['gapmin'] = min(gapNumber)
            output[key]['length'] = gapLength[0]
    return output



def run_orthomcl(spec2genes, indir, outdir, confile=None):
    """Run orthomcl-pipeline"""

    def group_file():
        return len(os.listdir(os.path.join(outdir,"groups"))) > 0

    if not confile:
        confile =  def_conf_file
    if not os.path.exists(indir):
        os.mkdir(indir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for spec in spec2genes.keys():
        SeqIO.write(spec2genes[spec], os.path.join(indir, spec.replace(' ', '_')+".fasta"), format="fasta")
    cmd = orthomclpipe + " -i "+indir + " -o " +outdir + " -m "+ confile +" --yes"
    pipe = TerminationPipe(cmd, 9999999)
    pipe.run()
    if pipe.failure or not group_file():
        cmd +=" --nocompliant"
        pipe = TerminationPipe(cmd, 9999999)
        pipe.run()
        if pipe.failure or not group_file():
            raise RuntimeError("Either orthomcl-pipeline failed, or it ran out of time")

    print "Orthomcl done"
    clusters = []
    with open(os.path.join(outdir, "groups/groups.txt")) as IN:
        for line in IN :
            data = (line.split(': ')[1]).split()
            d = ()
            for couple in data:
                genespec = (couple.split('|')[1])
                d += (genespec,)
            clusters.append(d)
    return clusters
