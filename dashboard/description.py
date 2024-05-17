amino_acid_features_hover_text = """Various protein annotation tools were used to predict sequence features.  
The amino acids are shown, and colored according to the tool prediction.  
    S = Signal Peptide Prediction  
    C = Cytoplasm Localized  
    I = Transmembrane segment inside cell  
    O = Segment localized outside cell  
    M = Transmembrane segment within membrane"""
hmmer_meta_features_text = """By default, we trim insertions with respect to the query (human) out of alignments.
To see full alignment with gaps and surrounding amino acid context, toggle off trimmed alignment.

    Heatmap features: \n
    match_to_query=Search found a perfect match to human transcript.\n
    tree_exists=Parent gene of query sequence has annotated orthologue relationships available.\n
    match_to_homologue=If tree_exists, HMMER match has an annotated orthologue relationship with human gene.\n
    start_in_window=Methionine found within 15 aa of human start codon in alignment.\n  
    stop_in_window=Stop found within 20 aa of human stop in alignment.\n
    fraction_identity=0.0-1.0, identity between sORF and putative orthologue."""