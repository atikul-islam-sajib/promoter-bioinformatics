single_nucleosides = ["A", "C", "G", "T"]

di_nucleosides = [
    "AA",
    "AC",
    "AG",
    "AT",
    "CA",
    "CC",
    "CG",
    "CT",
    "GA",
    "GC",
    "GG",
    "GT",
    "TA",
    "TC",
    "TG",
    "TT",
]

tri_nucleosides = [
    "AAA",
    "AAC",
    "AAG",
    "AAT",
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "AGA",
    "AGC",
    "AGG",
    "AGT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CAA",
    "CAC",
    "CAG",
    "CAT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CGA",
    "CGC",
    "CGG",
    "CGT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GAA",
    "GAC",
    "GAG",
    "GAT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GGA",
    "GGC",
    "GGG",
    "GGT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TAA",
    "TAC",
    "TAG",
    "TAT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TGA",
    "TGC",
    "TGG",
    "TGT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
]


tetra_nucleosides = [
    "AAAA",
    "AAAC",
    "AAAG",
    "AAAT",
    "AACA",
    "AACC",
    "AACG",
    "AACT",
    "AAGA",
    "AAGC",
    "AAGG",
    "AAGT",
    "AATA",
    "AATC",
    "AATG",
    "AATT",
    "ACAA",
    "ACAC",
    "ACAG",
    "ACAT",
    "ACCA",
    "ACCC",
    "ACCG",
    "ACCT",
    "ACGA",
    "ACGC",
    "ACGG",
    "ACGT",
    "ACTA",
    "ACTC",
    "ACTG",
    "ACTT",
    "AGAA",
    "AGAC",
    "AGAG",
    "AGAT",
    "AGCA",
    "AGCC",
    "AGCG",
    "AGCT",
    "AGGA",
    "AGGC",
    "AGGG",
    "AGGT",
    "AGTA",
    "AGTC",
    "AGTG",
    "AGTT",
    "ATAA",
    "ATAC",
    "ATAG",
    "ATAT",
    "ATCA",
    "ATCC",
    "ATCG",
    "ATCT",
    "ATGA",
    "ATGC",
    "ATGG",
    "ATGT",
    "ATTA",
    "ATTC",
    "ATTG",
    "ATTT",
    "CAAA",
    "CAAC",
    "CAAG",
    "CAAT",
    "CACA",
    "CACC",
    "CACG",
    "CACT",
    "CAGA",
    "CAGC",
    "CAGG",
    "CAGT",
    "CATA",
    "CATC",
    "CATG",
    "CATT",
    "CCAA",
    "CCAC",
    "CCAG",
    "CCAT",
    "CCCA",
    "CCCC",
    "CCCG",
    "CCCT",
    "CCGA",
    "CCGC",
    "CCGG",
    "CCGT",
    "CCTA",
    "CCTC",
    "CCTG",
    "CCTT",
    "CGAA",
    "CGAC",
    "CGAG",
    "CGAT",
    "CGCA",
    "CGCC",
    "CGCG",
    "CGCT",
    "CGGA",
    "CGGC",
    "CGGG",
    "CGGT",
    "CGTA",
    "CGTC",
    "CGTG",
    "CGTT",
    "CTAA",
    "CTAC",
    "CTAG",
    "CTAT",
    "CTCA",
    "CTCC",
    "CTCG",
    "CTCT",
    "CTGA",
    "CTGC",
    "CTGG",
    "CTGT",
    "CTTA",
    "CTTC",
    "CTTG",
    "CTTT",
    "GAAA",
    "GAAC",
    "GAAG",
    "GAAT",
    "GACA",
    "GACC",
    "GACG",
    "GACT",
    "GAGA",
    "GAGC",
    "GAGG",
    "GAGT",
    "GATA",
    "GATC",
    "GATG",
    "GATT",
    "GCAA",
    "GCAC",
    "GCAG",
    "GCAT",
    "GCCA",
    "GCCC",
    "GCCG",
    "GCCT",
    "GCGA",
    "GCGC",
    "GCGG",
    "GCGT",
    "GCTA",
    "GCTC",
    "GCTG",
    "GCTT",
    "GGAA",
    "GGAC",
    "GGAG",
    "GGAT",
    "GGCA",
    "GGCC",
    "GGCG",
    "GGCT",
    "GGGA",
    "GGGC",
    "GGGG",
    "GGGT",
    "GGTA",
    "GGTC",
    "GGTG",
    "GGTT",
    "GTAA",
    "GTAC",
    "GTAG",
    "GTAT",
    "GTCA",
    "GTCC",
    "GTCG",
    "GTCT",
    "GTGA",
    "GTGC",
    "GTGG",
    "GTGT",
    "GTTA",
    "GTTC",
    "GTTG",
    "GTTT",
    "TAAA",
    "TAAC",
    "TAAG",
    "TAAT",
    "TACA",
    "TACC",
    "TACG",
    "TACT",
    "TAGA",
    "TAGC",
    "TAGG",
    "TAGT",
    "TATA",
    "TATC",
    "TATG",
    "TATT",
    "TCAA",
    "TCAC",
    "TCAG",
    "TCAT",
    "TCCA",
    "TCCC",
    "TCCG",
    "TCCT",
    "TCGA",
    "TCGC",
    "TCGG",
    "TCGT",
    "TCTA",
    "TCTC",
    "TCTG",
    "TCTT",
    "TGAA",
    "TGAC",
    "TGAG",
    "TGAT",
    "TGCA",
    "TGCC",
    "TGCG",
    "TGCT",
    "TGGA",
    "TGGC",
    "TGGG",
    "TGGT",
    "TGTA",
    "TGTC",
    "TGTG",
    "TGTT",
    "TTAA",
    "TTAC",
    "TTAG",
    "TTAT",
    "TTCA",
    "TTCC",
    "TTCG",
    "TTCT",
    "TTGA",
    "TTGC",
    "TTGG",
    "TTGT",
    "TTTA",
    "TTTC",
    "TTTG",
    "TTTT",
]
