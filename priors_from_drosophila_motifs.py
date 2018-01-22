# Written by Alex Rives

import os, sys
import cStringIO as StringIO
import numpy as np
import scipy.sparse as sp
import pandas

ndeveaux_reference="/mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster"
genome = os.path.join(ndeveaux_reference, "dm6.fa")
gtf = os.path.join(ndeveaux_reference, "genes.gtf")

# Copied over survey meme from /mnt/home/victle/Drosophila_inf/combined.meme.txt
annotationfile = os.path.join(ndeveaux_reference, "survey.meme.txt")

class Annotation(object):

    def __init__(self):
        self.data = {}
        self.keys = {}
        self.scan = {
            'MOTIF_ID': 'motif',
            'TF_ID': 'eid',
            'TF_NAME': 'name',
        }
        self.primarykey = 'motif'
        self.errorlines = []

    def get_id(self, motif):
        if motif in self.data:
            return self.data[motif]['eid']
        else:
            return None

    def parse_header(self, line):
        line = line.strip()
        line = line.split()
        for i, colname in enumerate(line):
            if colname in self.scan:
                if self.scan[colname] not in self.keys:
                    self.keys[self.scan[colname]] = i
        error = False
        for colname, key in self.scan.iteritems():
            if key not in self.keys:
                print "Missing: %s" % (colname)
                error = True
        if error:
            raise Exception("Error parsing header")


    def parse_line(self, line):
        line = line.rstrip('\n')
        line = line.split('\t')
        if len(line) == 0:
            return
        d = {}
        for key, i in self.keys.iteritems():
            try:
                d[key] = line[i]
            except:
                print line
        key = d[self.primarykey]
        if key in self.data:
            self.errorlines.append(line)
#             raise Exception("Primary key %s is already in db, %d" % (key, len(self.data)))
        self.data[key] = d


class TranscriptionFactor(object):

    def __init__(self, uid):
        self.uid = uid
        self.motifs = {}
        self.targets = {}

    def add_target(self, motif, target_id, target_dict):
        if motif not in self.motifs:
            self.motifs[motif] = {}
        self.motifs[motif][target_id] = True
        if target_id in self.targets.keys():
            self.targets[target_id]['count'] += 1
        else:
            self.targets[target_id] = target_dict
            self.targets[target_id]['count'] = 1

    def __str__(self):
        return "<TF:%s>M:%dT:%d" % (self.uid, len(self.motifs), len(self.targets))


class TFNetwork(object):

    def __init__(self, annotation_path):
        self.txn_factors = {}
        self.targets = {}
        self.motif_annotation_path = annotation_path
        self.motif_annotation = Annotation()

    def load_motif_annotation(self):
        with open(self.motif_annotation_path, 'r') as handle:
            header = handle.readline()
            print header
            self.motif_annotation.parse_header(header)
            for line in handle:
                self.motif_annotation.parse_line(line)

    def get_txn_factor(self, uid):
        if uid in self.txn_factors:
            return self.txn_factors[uid]
        else:
            tf = TranscriptionFactor(uid)
            self.txn_factors[uid] = tf
            return tf


    def add_target(self, motif, gtf_dict):
#         print motif
        uid = self.motif_annotation.get_id(motif)
        tf = self.get_txn_factor(uid)
        target_id = gtf_dict['gene_id']
        tf.add_target(motif, target_id, gtf_dict)
        self.targets[target_id] = True


    def parse_line(self, line):
        line = line.strip(';\n')
        line = line.split('\t')
        motif = line[3]
        target_line = line[9]
        target_dict = {'gene_id': target_line}
#         print motif
#         print target_dict
        self.add_target(motif, target_dict)

def main():

    nw = TFNetwork(annotationfile)
    nw.load_motif_annotation()

    files = os.listdir('targets')
    total = len(files)
    i = 1
    for f in files:
        if i % 100 == 0:
            print i
        with open(os.path.join('targets', f), 'r') as handle:
            for line in handle:
                nw.parse_line(line)
        i+=1


    colnames = nw.txn_factors.keys()
    rownames = nw.targets.keys()

    colkeys = {}
    rowkeys = {}
    for i, gid in enumerate(colnames):
        colkeys[gid] = i
    for i, gid in enumerate(rownames):
        rowkeys[gid] = i


    edges = 0
    row_indices = []
    col_indices = []
    data = []
    w = 0
    for tf_id, tf in nw.txn_factors.iteritems():
        for target_id in tf.targets:
            tf_index = colkeys[tf_id]
            target_index = rowkeys[target_id]
            row_indices.append(target_index)
            col_indices.append(tf_index)
            data.append(tf.targets[target_id]['count'])
            edges+=1


    print "Edges: %d" % (edges)

    row_indices = np.array(row_indices, dtype=int)
    col_indices = np.array(col_indices, dtype=int)
    data = np.array(data, dtype=int)

    # from scipy.sparse import coo_matrix
    n_genes = len(nw.targets)
    n_txn_factors = len(nw.txn_factors)
    M = sp.coo_matrix((data, (row_indices, col_indices)), shape=(n_genes, n_txn_factors), dtype=np.int8)

    df = pandas.DataFrame(M.toarray())
    df.index = map(lambda x: x.split(".")[0], rownames)
    df.columns = colnames

    output = StringIO.StringIO()
    df.to_csv(output, sep='\t')

    out_dir = 'inferelator_input'
    os.mkdir(out_dir)

    with open(os.path.join(out_dir, "priors.tsv"), 'w') as handle:
        handle.write(output.getvalue()[1:])

    # save transcription factors
    # to txt file
    with open(os.path.join(out_dir, 'txn_factor_names.tsv'), 'w') as handle:
        for geneid in colnames:
            print >> handle, geneid


if __name__ == '__main__':
    main()