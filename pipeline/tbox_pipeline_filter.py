#tbox_pipeline_filter.py
#By Merrick Pierson Smela
#Filters output columns of T-box predictions
#As currently configured, doesn't remove anything

import pandas as pd
import sys

col_list = ["Name", "FASTA_sequence", "Rank", "E_value", "Score", "Bias", "Tbox_start", "Tbox_end", "CM_accuracy", "GC", "Sequence", "Structure", "s1_start", "s1_loop_start", "s1_loop_end", "s1_end", "antiterm_start", "antiterm_end", "term_start", "term_end", "codon_start", "codon_end", "codon", "codon_region", "discrim_start", "discrim_end", "discriminator", "warnings", "type", "source", "whole_antiterm_structure", "other_stems", "whole_antiterm_warnings", "term_sequence", "term_structure", "terminator_energy", "term_errors", "antiterm_term_sequence", "infernal_antiterminator_structure", "vienna_antiterminator_structure", "vienna_antiterminator_energy", "vienna_antiterminator_errors", "terminator_structure", "terminator_errors", "new_term_structure", "new_term_energy", "new_term_errors", "whole_term_structure", "folded_antiterm_structure", "Trimmed_sequence", "Trimmed_antiterm_struct", "Trimmed_term_struct", "hash_string", "unique_name", "accession_url", "accession_name", "locus_start", "tbox_length", "locus_end", "locus_view_start", "locus_view_end", "deltadelta_g", "TaxId", "GBSeq_organism", "phylum", "class", "order", "family", "genus", "downstream_protein", "downstream_protein_id", "downstream_protein_EC", "protein_desc", "refine_codon", "refine_codon_io", "refine_codon_code", "refine_codon_top", "refine_codon_alt_1", "refine_codon_alt_2", "refine_codon_num", "amino_acid_top", "trna_family_top", "trna_seq_top", "trna_struc_top", "amino_acid_alt_1", "trna_family_alt_1", "trna_seq_alt_1", "trna_struc_alt_1", "amino_acid_alt_2", "trna_family_alt_2", "trna_seq_alt_2", "trna_struc_alt_2"]

tboxes = pd.read_csv(sys.argv[1], usecols = col_list)
tboxes.to_csv(sys.argv[2], index = False)