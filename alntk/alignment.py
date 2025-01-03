import numpy as np
from numpy.lib.stride_tricks import as_strided
from Bio import SeqIO

class Alignment:
    def __init__(self):
        np.random.seed(42)

    def import_from_fasta(self, fasta_file):
        """
        Import an array of descriptions and an array of sequences from a FASTA file.
        Save them as attributes of the Alignment instance.
    
        Parameters
        ----------
        fasta_file: str
            filename of FASTA file, typically `.fasta`, `.faa` or `.fa`.
        """
        seqs = []
        descs = []
        
        for record in SeqIO.parse(fasta_file, 'fasta'):
            seqs.append(str(record.seq))
            descs.append(record.description.split('/')[0])
    
        self.descs_arr = np.array(descs)
        self.seqs_arr = np.array([[residue for residue in seq] for seq in seqs])
        self.seq_idxs0 = np.arange(0, self.seqs_arr.shape[0])
        self.pos_idxs0 = np.arange(0, self.seqs_arr.shape[1])
        self.seq_idxs = self.seq_idxs0
        self.pos_idxs = self.pos_idxs0
        self.seq_untrim_idxs = self.seq_idxs0

    def subsample(self, n_subsample):
        """
        Randomly choose a subsample for faster computation of statistics.
        Sequences are not deleted, just the list of indices `seq_idxs` is updated.
        Generator is seeded at 42 for reproducibility.

        Parameters
        ----------
        n_subsample: int
            number of sequences that we want to subsample from the original alignment
        """
        self.seq_idxs = np.random.randint(0, self.seqs_arr.shape[0], size=n_subsample)

    def drop(self, del_idxs):
        aln_seqs = self.get_seqs()
        drop_seq_trim_idxs = np.where(np.sum(aln_seqs[:, del_idxs[1]] != '-', axis=1) > 0)[0]
        self.seq_untrim_idxs = np.delete(self.seq_untrim_idxs, drop_seq_trim_idxs)
        self.seq_idxs = np.delete(self.seq_idxs, del_idxs[0])
        self.pos_idxs = np.delete(self.pos_idxs, del_idxs[1])

    def reset_seq(self):
        self.seq_idxs = self.seq_idxs0

    def reset_pos(self):
        self.pos_idxs = self.pos_idxs0
        self.seq_untrim_idxs = self.seq_idxs0

    def get_seqs(self):
        return self.seqs_arr[self.seq_idxs, :][:, self.pos_idxs]

    def get_descs(self):
        return self.descs_arr[self.seq_idxs]

    def print_report(self, trimmed=False):
        print(f"Number of sequences: {len(self.seq_idxs)}")
        print(f"Number of positions: {len(self.pos_idxs)}")
        if trimmed:
            print(f"Number of untrimmed sequences: {len(self.seq_untrim_idxs)}")

    def get_seq_gap(self):
        """
        """
        aln_seqs = self.get_seqs()
        gaps_per_seq = np.sum(aln_seqs == '-', axis=1) / aln_seqs.shape[1]
    
        return gaps_per_seq
    
    def get_pos_gap(self):
        """
        """
        aln_seqs = self.get_seqs()
        gaps_per_pos = np.sum(aln_seqs == '-', axis=0) / aln_seqs.shape[0]
    
        return gaps_per_pos

def filter_gappy(aln, threshold=0.8):
    """
    """
    aln_seqs = aln.get_seqs()

    col_tbd = np.where(np.sum(aln_seqs == '-', axis=0) / aln_seqs.shape[0] > threshold)
    return [], col_tbd

def filter_length(aln, min_len, max_len):
    """
    """
    aln_seqs = aln.get_seqs()
    
    cond_min = np.sum(aln_seqs != '-', axis=1) > min_len
    cond_max = np.sum(aln_seqs != '-', axis=1) < max_len
    seq_tbd = np.where(cond_min & cond_max == False)[0]
    return seq_tbd, []

def find_pattern(text, pattern):
    """
    """
    pattern = np.asarray(list(pattern))
    
    if len(pattern) > len(text):
        return np.array([], dtype=int)
    
    # Create sliding window view of text
    windows = np.lib.stride_tricks.as_strided(
        text, 
        shape=(len(text) - len(pattern) + 1, len(pattern)),
        strides=(text.strides[0], text.strides[0])
    )
    
    # Find indices where windows match the pattern
    match_indices = np.where(np.all(windows == pattern, axis=1))[0]
    
    return match_indices

def find_sequence(aln, accession):
    aln_seqs = aln.get_seqs()
    aln_descs = aln.get_descs()

    ref_idx = np.where(np.array([e.split('|')[2] for e in aln_descs]) == accession)[0][0]

    return aln_seqs[ref_idx]

def filter_residue(aln, ref_seq_acc, ref_pattern, residue):
    """
    """
    aln_seqs = aln.get_seqs()
    aln_descs = aln.get_descs()

    ref_idx = np.where(np.array([e.split('|')[2] for e in aln_descs]) == ref_seq_acc)[0][0]
    ref_seq = aln_seqs[ref_idx]
    ref_pattern_idxs = find_pattern(ref_seq, ref_pattern)

    # delete sequences that do not have a specific residue in the first position of the pattern
    seq_tbd = np.where(aln_seqs[:, ref_pattern_idxs[0]] != residue)[0]

    return seq_tbd, []

def filter_residues_any(aln, ref_seq_acc, ref_pattern):
    """
    """
    aln_seqs = aln.get_seqs()
    aln_descs = aln.get_descs()

    ref_idx = np.where(np.array([e.split('|')[2] for e in aln_descs]) == ref_seq_acc)[0][0]
    ref_seq = aln_seqs[ref_idx]
    ref_pattern_idxs = find_pattern(ref_seq, ref_pattern)

    # delete sequences that have a gap in the columns corresponding
    # to the first occurence of the reference pattern in the reference sequence
    seq_tbd = np.where(np.sum(aln_seqs[:, ref_pattern_idxs[0]:ref_pattern_idxs[0]+len(ref_pattern)] == '-', axis=1) > 0)[0]
    pos_tbd = np.arange(0, ref_pattern_idxs[0])

    return seq_tbd, pos_tbd

def filter_ambiguous(aln):
    """
    """
    aln_seqs = aln.get_seqs()
    seq_tbd = np.where(np.sum((aln_seqs == 'B') + (aln_seqs == 'J') + (aln_seqs == 'X') + (aln_seqs == 'Z'), axis=1) > 0)[0]
    return seq_tbd, []

def filter_compact(aln, gap_threshold_ratio=0.95):
    """
    Compactify the alignment by deleting sequences and consequently
    deleting entirely gapped columns.
    """
    aln_seqs = aln.get_seqs()

    # Get the gap threshold in terms of integer number of gaps from the ratio
    gap_threshold = int(gap_threshold_ratio*aln_seqs.shape[0])

    # Get the gap / no gap boolean matrix
    aln_seqs_isgap = (aln_seqs == '-')
    
    # Count the number of gaps in each column
    gaps_per_column = np.sum(aln_seqs_isgap, axis=0)

    # Get the index of all the columns that have more gaps than the threshold
    # We call these "positions to be deleted"
    pos_tbd = np.where(gaps_per_column > gap_threshold)[0]

    # Count how many gaps do sequences have in "gappy positions"
    # We call these counts of "unusual amino acids"
    # The sequences that have these amino acids are "sequences to be deleted"
    unusual_aa_per_sequence = np.sum(aln_seqs[:, pos_tbd] != '-', axis=1)
    seq_tbd = np.where(unusual_aa_per_sequence > 0)[0]

    return seq_tbd, pos_tbd 

def write_to_fasta(aln, fasta_file):
    """
    Write alignment to a FASTA file.
    """
    aln_descs = aln.get_descs()
    aln_seqs = aln.get_seqs()
    with open(fasta_file, 'w') as f:
        for d, s in zip(aln_descs, aln_seqs):
            s = ''.join(s)
            f.write(f'>{d}\n')
            f.write(s)
            f.write('\n\n')