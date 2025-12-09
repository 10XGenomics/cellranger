//! Mappings from symbol sample index set names to sample index sequences
//! for official 10x sample index plate products.
#![expect(missing_docs)]

use lazy_static::lazy_static;
use std::collections::HashMap;

lazy_static! {
    pub static ref SAMPLE_INDEX_MAP: HashMap<&'static str, [&'static str; 4]> = {
        let mut m = HashMap::default();
        m.insert("220027", ["TCGCCATA", "GTATACAC", "AATGGTGG", "CGCATGCT"]);
        m.insert("220028", ["TATCCTCG", "GCGAGGTC", "CGCTTCAA", "ATAGAAGT"]);
        m.insert("220029", ["TGACGTCG", "CTTGTGTA", "ACGACCGT", "GACTAAAC"]);
        m.insert("220030", ["ATCTAGCT", "GAGCGTAC", "TCAGCCTG", "CGTATAGA"]);
        m.insert("220031", ["CCGTTCCC", "ATACAGTT", "TGTAGTAA", "GACGCAGG"]);
        m.insert("220032", ["TCAATTGG", "AGTTAGAA", "GAGCGCTT", "CTCGCACC"]);
        m.insert("220033", ["CTGCCTTG", "ACTAGCCC", "GGCTAGAT", "TAAGTAGA"]);
        m.insert("220034", ["GGCAGAAA", "ACGGTTCT", "CATTCGTC", "TTACACGG"]);
        m.insert("SI-001", ["TCGCCATA", "GTATACAC", "AATGGTGG", "CGCATGCT"]);
        m.insert("SI-002", ["TATCCTCG", "GCGAGGTC", "CGCTTCAA", "ATAGAAGT"]);
        m.insert("SI-003", ["TGACGTCG", "CTTGTGTA", "ACGACCGT", "GACTAAAC"]);
        m.insert("SI-004", ["ATCTAGCT", "GAGCGTAC", "TCAGCCTG", "CGTATAGA"]);
        m.insert("SI-005", ["CCGTTCCC", "ATACAGTT", "TGTAGTAA", "GACGCAGG"]);
        m.insert("SI-006", ["TCAATTGG", "AGTTAGAA", "GAGCGCTT", "CTCGCACC"]);
        m.insert("SI-007", ["CTGCCTTG", "ACTAGCCC", "GGCTAGAT", "TAAGTAGA"]);
        m.insert("SI-008", ["GGCAGAAA", "ACGGTTCT", "CATTCGTC", "TTACACGG"]);
        m.insert("SI-3A-A1", ["AAACGGCG", "CCTACCAT", "GGCGTTTC", "TTGTAAGA"]);
        m.insert(
            "SI-3A-A10",
            ["ACAGCAAC", "CGCAATTT", "GAGTTGCG", "TTTCGCGA"],
        );
        m.insert(
            "SI-3A-A11",
            ["ACCAGTCC", "CTTTCCTT", "GGACAGGG", "TAGGTAAA"],
        );
        m.insert(
            "SI-3A-A12",
            ["ACTACTGT", "CGGGAACG", "GACCTCTC", "TTATGGAA"],
        );
        m.insert("SI-3A-A2", ["AGCCCTTT", "CAAGTCCA", "GTGAGAAG", "TCTTAGGC"]);
        m.insert("SI-3A-A3", ["AAAGCATA", "CTGCAGCC", "GCCTTTAT", "TGTAGCGG"]);
        m.insert("SI-3A-A4", ["AGAACGCC", "CATGGCAG", "GTCTTTGA", "TCGCAATT"]);
        m.insert("SI-3A-A5", ["ATTGGGAA", "CAGTCTGG", "GGCATACT", "TCACACTC"]);
        m.insert("SI-3A-A6", ["ACGGGACT", "CTTTCGAC", "GAACATGA", "TGCATCTG"]);
        m.insert("SI-3A-A7", ["AGGTCATA", "CTCATCAT", "GCTGAGGG", "TAACGTCC"]);
        m.insert("SI-3A-A8", ["ATGATACG", "CCACAGAA", "GACTGTTC", "TGTGCCGT"]);
        m.insert("SI-3A-A9", ["ACAACTTG", "CTCCAACA", "GAGTGCGT", "TGTGTGAC"]);
        m.insert("SI-3A-B1", ["AGGCTACC", "CTAGCTGT", "GCCAACAA", "TATTGGTG"]);
        m.insert(
            "SI-3A-B10",
            ["ACCATTAA", "CTGGACGT", "GAACGGTC", "TGTTCACG"],
        );
        m.insert(
            "SI-3A-B11",
            ["ATGGTCGC", "CGACATAG", "GATTCGCT", "TCCAGATA"],
        );
        m.insert(
            "SI-3A-B12",
            ["ACGCTTGG", "CGCTACAT", "GAAAGACA", "TTTGCGTC"],
        );
        m.insert("SI-3A-B2", ["AAGTTGAT", "CCCACCCA", "GGTCGAGC", "TTAGATTG"]);
        m.insert("SI-3A-B3", ["ATTGGACG", "CAGCTTAC", "GGCAAGGA", "TCATCCTT"]);
        m.insert("SI-3A-B4", ["AGGGACTG", "CCTCTAAC", "GACAGGCT", "TTATCTGA"]);
        m.insert("SI-3A-B5", ["ATCGTACT", "CATCAGTG", "GGGACTAC", "TCATGCGA"]);
        m.insert("SI-3A-B6", ["AACGCGAA", "CTATTTGG", "GCGCACCT", "TGTAGATC"]);
        m.insert("SI-3A-B7", ["AGGGATGA", "CTTCTGTT", "GAATGCAC", "TCCACACG"]);
        m.insert("SI-3A-B8", ["ACGTTCAC", "CAAGGTCT", "GTTAAGTG", "TGCCCAGA"]);
        m.insert("SI-3A-B9", ["AAGCGTGT", "CTTGACCG", "GCCACGTA", "TGATTAAC"]);
        m.insert("SI-3A-C1", ["AGACTTTC", "CCGAGGCA", "GATGCAGT", "TTCTACAG"]);
        m.insert(
            "SI-3A-C10",
            ["ATCTGATC", "CGTGCTAA", "GAGAAGGG", "TCACTCCT"],
        );
        m.insert(
            "SI-3A-C11",
            ["ACCGAACA", "CGACTCTT", "GTTTGTGG", "TAGACGAC"],
        );
        m.insert(
            "SI-3A-C12",
            ["ATCCGGCA", "CCGTTATG", "GGTAATGT", "TAAGCCAC"],
        );
        m.insert("SI-3A-C2", ["AATCACTA", "CCGAGAAC", "GTAGTGCG", "TGCTCTGT"]);
        m.insert("SI-3A-C3", ["ACGTTACA", "CGTAGGTT", "GACGACGG", "TTACCTAC"]);
        m.insert("SI-3A-C4", ["ACATTGGC", "CTTAGTCA", "GAGCCCAT", "TGCGAATG"]);
        m.insert("SI-3A-C5", ["ATGCATTC", "CACTGACT", "GGTACGGG", "TCAGTCAA"]);
        m.insert("SI-3A-C6", ["ACTCAGAC", "CGCTCAGG", "GAGGTTTA", "TTAAGCCT"]);
        m.insert("SI-3A-C7", ["ACACCGGG", "CATAATCC", "GGCGGAAT", "TTGTTCTA"]);
        m.insert("SI-3A-C8", ["AGCTCGAG", "CAGGAAGA", "GCACGTTT", "TTTATCCC"]);
        m.insert("SI-3A-C9", ["AGATCGGT", "CATCGTCG", "GTCATATA", "TCGGACAC"]);
        m.insert("SI-3A-D1", ["AGCTGCGT", "CAACCATC", "GTGGAGCA", "TCTATTAG"]);
        m.insert(
            "SI-3A-D10",
            ["AGATAACA", "CTTATTTG", "GCGGGCAT", "TACCCGGC"],
        );
        m.insert(
            "SI-3A-D11",
            ["ATATGAGA", "CACCTCAG", "GCTACTTC", "TGGGAGCT"],
        );
        m.insert(
            "SI-3A-D12",
            ["AGAAACGT", "CACTCAAC", "GCTGTGTA", "TTGCGTCG"],
        );
        m.insert("SI-3A-D2", ["ACATTCCG", "CTGCGGTA", "GACACAAT", "TGTGATGC"]);
        m.insert("SI-3A-D3", ["ACTTCACT", "CGAAGTTG", "GAGCACGC", "TTCGTGAA"]);
        m.insert("SI-3A-D4", ["AAATCGTC", "CTTCGAAT", "GCGATCGG", "TGCGATCA"]);
        m.insert("SI-3A-D5", ["AGACGGAT", "CCTTTAGA", "GTCGACTC", "TAGACTCG"]);
        m.insert("SI-3A-D6", ["ATGCCAAA", "CCTTATCG", "GAAGTCTT", "TGCAGGGC"]);
        m.insert("SI-3A-D7", ["AACTTAGA", "CCGGATCC", "GGTCGCAT", "TTAACGTG"]);
        m.insert("SI-3A-D8", ["AATCTTTG", "CTCAAGAC", "GGATGAGT", "TCGGCCCA"]);
        m.insert("SI-3A-D9", ["ACCTACTG", "CAAGGGAC", "GGGACACA", "TTTCTTGT"]);
        m.insert("SI-3A-E1", ["ACGAAAGC", "CGCCCGTA", "GTTTGCCT", "TAAGTTAG"]);
        m.insert(
            "SI-3A-E10",
            ["ATTGTTTC", "CGCAGGAG", "GCACCAGT", "TAGTACCA"],
        );
        m.insert(
            "SI-3A-E11",
            ["ATCGCCAT", "CATAAAGG", "GGGTTTCC", "TCACGGTA"],
        );
        m.insert(
            "SI-3A-E12",
            ["ACGCGGAA", "CGCTATCC", "GTTGCATG", "TAAATCGT"],
        );
        m.insert("SI-3A-E2", ["AGGCTGGT", "CACAACTA", "GTTGGTCC", "TCATCAAG"]);
        m.insert("SI-3A-E3", ["AACAAGTC", "CGGCTCCA", "GTATGTAT", "TCTGCAGG"]);
        m.insert("SI-3A-E4", ["AGCTGACG", "CCGGTGTC", "GTAAACAT", "TATCCTGA"]);
        m.insert("SI-3A-E5", ["ATCCAAGG", "CCGTTGAA", "GGAAGCTC", "TATGCTCT"]);
        m.insert("SI-3A-E6", ["ATTGAAAC", "CAGCCCGA", "GCCATTTG", "TGATGGCT"]);
        m.insert("SI-3A-E7", ["AAGACGTG", "CCATGTGT", "GTTCACAA", "TGCGTACC"]);
        m.insert("SI-3A-E8", ["AGCCTATG", "CTAACGCA", "GCTTACAT", "TAGGGTGC"]);
        m.insert("SI-3A-E9", ["AGTAAGCA", "CCGGTAAT", "GTATCTTC", "TACCGCGG"]);
        m.insert("SI-3A-F1", ["ATCGCTCC", "CCGTACAG", "GATAGGTA", "TGACTAGT"]);
        m.insert(
            "SI-3A-F10",
            ["ATTCGTGC", "CGCGTGCA", "GAATACTG", "TCGACAAT"],
        );
        m.insert(
            "SI-3A-F11",
            ["AGCAGTTA", "CTTGTACC", "GAACCCGG", "TCGTAGAT"],
        );
        m.insert(
            "SI-3A-F12",
            ["AATTGAAC", "CCAGTGGA", "GTCCATTG", "TGGACCCT"],
        );
        m.insert("SI-3A-F2", ["ATGGTTAG", "CATTGATA", "GCAAACGC", "TGCCCGCT"]);
        m.insert("SI-3A-F3", ["AGTCTGTA", "CAGAATAG", "GCCTCCGT", "TTAGGACC"]);
        m.insert("SI-3A-F4", ["AACGACAC", "CGTCCTCT", "GCATGATA", "TTGATGGG"]);
        m.insert("SI-3A-F5", ["AACACAGC", "CGGTTTAG", "GTACGGCT", "TCTGACTA"]);
        m.insert("SI-3A-F6", ["ATGCCGGC", "CCTAATTA", "GACTTCCT", "TGAGGAAG"]);
        m.insert("SI-3A-F7", ["ACCCGAGA", "CAAACTTT", "GGTTAGAC", "TTGGTCCG"]);
        m.insert("SI-3A-F8", ["AGTTGGGA", "CCAGAAAG", "GTGCCCTC", "TACATTCT"]);
        m.insert("SI-3A-F9", ["AGTTAGTT", "CACGCACG", "GTACTTAA", "TCGAGCGC"]);
        m.insert("SI-3A-G1", ["ATGCGATT", "CATATGCG", "GGATACGA", "TCCGCTAC"]);
        m.insert(
            "SI-3A-G10",
            ["ATACTGAG", "CGGAGACT", "GATGCCTC", "TCCTATGA"],
        );
        m.insert(
            "SI-3A-G11",
            ["AGGGCGTT", "CTATACGC", "GCTCGTCA", "TACATAAG"],
        );
        m.insert(
            "SI-3A-G12",
            ["ACCCGCAC", "CATGCGTA", "GTGATAGT", "TGATATCG"],
        );
        m.insert("SI-3A-G2", ["ATAACCTA", "CGGTGAGC", "GATCTTAT", "TCCGAGCG"]);
        m.insert("SI-3A-G3", ["ATGTCCAG", "CGACGTCA", "GCTATAGC", "TACGAGTT"]);
        m.insert("SI-3A-G4", ["AGCTTCTC", "CCTGCGGT", "GTACAACG", "TAGAGTAA"]);
        m.insert("SI-3A-G5", ["ATGAAGTA", "CGCCGAAC", "GAAGCTCG", "TCTTTCGT"]);
        m.insert("SI-3A-G6", ["AGCACTGG", "CATTACAC", "GTGCGACA", "TCAGTGTT"]);
        m.insert("SI-3A-G7", ["ATTACCGG", "CAGTAATT", "GCCGGTAA", "TGACTGCC"]);
        m.insert("SI-3A-G8", ["AAGTACTC", "CTTGGAGA", "GGAACTCT", "TCCCTGAG"]);
        m.insert("SI-3A-G9", ["AGTCTCAG", "CAATGGCA", "GCCGAAGT", "TTGACTTC"]);
        m.insert("SI-3A-H1", ["AAACTCAT", "CGGGAGTA", "GTCACAGG", "TCTTGTCC"]);
        m.insert(
            "SI-3A-H10",
            ["ATTTCAGC", "CGAGTGAT", "GACCGCCA", "TCGAATTG"],
        );
        m.insert(
            "SI-3A-H11",
            ["AGGATCGA", "CACGATTC", "GTATCGAG", "TCTCGACT"],
        );
        m.insert(
            "SI-3A-H12",
            ["ACGAGTAG", "CAATCCCT", "GTCCAGGC", "TGTGTATA"],
        );
        m.insert("SI-3A-H2", ["AACGGTCA", "CCGAACTC", "GGTCCAAG", "TTATTGGT"]);
        m.insert("SI-3A-H3", ["ACACCTAA", "CGTTTGGG", "GACAAACC", "TTGGGCTT"]);
        m.insert("SI-3A-H4", ["ACTGGAGC", "CGGTCGTG", "GAAATCAA", "TTCCATCT"]);
        m.insert("SI-3A-H5", ["ATAGTATG", "CCGCGTCT", "GGCTCCAC", "TATAAGGA"]);
        m.insert("SI-3A-H6", ["AAGCATAA", "CCCATCGC", "GGTTGATG", "TTAGCGCT"]);
        m.insert("SI-3A-H7", ["AACGGGTG", "CTAATTCT", "GCTTCAAC", "TGGCACGA"]);
        m.insert("SI-3A-H8", ["AAGAGCGG", "CTTGTTAT", "GGCCCATC", "TCATAGCA"]);
        m.insert("SI-3A-H9", ["ACCTGCCA", "CTTCATAC", "GGAATATG", "TAGGCGGT"]);
        m.insert("SI-GA-A1", ["GGTTTACT", "CTAAACGG", "TCGGCGTC", "AACCGTAA"]);
        m.insert(
            "SI-GA-A10",
            ["GAAACCCT", "TTTCTGTC", "CCGTGTGA", "AGCGAAAG"],
        );
        m.insert(
            "SI-GA-A11",
            ["GTCCGGTC", "AAGATCAT", "CCTGAAGG", "TGATCTCA"],
        );
        m.insert(
            "SI-GA-A12",
            ["AGTGGAAC", "GTCTCCTT", "TCACATCA", "CAGATGGG"],
        );
        m.insert("SI-GA-A2", ["TTTCATGA", "ACGTCCCT", "CGCATGTG", "GAAGGAAC"]);
        m.insert("SI-GA-A3", ["CAGTACTG", "AGTAGTCT", "GCAGTAGA", "TTCCCGAC"]);
        m.insert("SI-GA-A4", ["TATGATTC", "CCCACAGT", "ATGCTGAA", "GGATGCCG"]);
        m.insert("SI-GA-A5", ["CTAGGTGA", "TCGTTCAG", "AGCCAATT", "GATACGCC"]);
        m.insert("SI-GA-A6", ["CGCTATGT", "GCTGTCCA", "TTGAGATC", "AAACCGAG"]);
        m.insert("SI-GA-A7", ["ACAGAGGT", "TATAGTTG", "CGGTCCCA", "GTCCTAAC"]);
        m.insert("SI-GA-A8", ["GCATCTCC", "TGTAAGGT", "CTGCGATG", "AACGTCAA"]);
        m.insert("SI-GA-A9", ["TCTTAAAG", "CGAGGCTC", "GTCCTTCT", "AAGACGGA"]);
        m.insert("SI-GA-B1", ["GTAATCTT", "TCCGGAAG", "AGTTCGGC", "CAGCATCA"]);
        m.insert(
            "SI-GA-B10",
            ["ACCGTATG", "GATTAGAT", "CTGACTGA", "TGACGCCC"],
        );
        m.insert(
            "SI-GA-B11",
            ["GTTCCTCA", "AGGTACGC", "TAAGTATG", "CCCAGGAT"],
        );
        m.insert(
            "SI-GA-B12",
            ["TACCACCA", "CTAAGTTT", "GGGTCAAG", "ACTGTGGC"],
        );
        m.insert("SI-GA-B2", ["TACTCTTC", "CCTGTGCG", "GGACACGT", "ATGAGAAA"]);
        m.insert("SI-GA-B3", ["GTGTATTA", "TGTGCGGG", "ACCATAAC", "CAACGCCT"]);
        m.insert("SI-GA-B4", ["ACTTCATA", "GAGATGAC", "TGCCGTGG", "CTAGACCT"]);
        m.insert("SI-GA-B5", ["AATAATGG", "CCAGGGCA", "TGCCTCAT", "GTGTCATC"]);
        m.insert("SI-GA-B6", ["CGTTAATC", "GCCACGCT", "TTACTCAG", "AAGGGTGA"]);
        m.insert("SI-GA-B7", ["AAACCTCA", "GCCTTGGT", "CTGGACTC", "TGTAGAAG"]);
        m.insert("SI-GA-B8", ["AAAGTGCT", "GCTACCTG", "TGCTGTAA", "CTGCAAGC"]);
        m.insert("SI-GA-B9", ["CTGTAACT", "TCTAGCGA", "AGAGTGTG", "GACCCTAC"]);
        m.insert("SI-GA-C1", ["CCACTTAT", "AACTGGCG", "TTGGCATA", "GGTAACGC"]);
        m.insert(
            "SI-GA-C10",
            ["TCTCAGTG", "GAGACTAT", "CGCTTAGC", "ATAGGCCA"],
        );
        m.insert(
            "SI-GA-C11",
            ["GAGGATCT", "AGACCATA", "TCCTGCGC", "CTTATGAG"],
        );
        m.insert(
            "SI-GA-C12",
            ["TCTCGTTT", "GGCTAGCG", "ATGACCGC", "CAAGTAAA"],
        );
        m.insert("SI-GA-C2", ["CCTAGACC", "ATCTCTGT", "TAGCTCTA", "GGAGAGAG"]);
        m.insert("SI-GA-C3", ["TCAGCCGT", "CAGAGGCC", "GGTCAATA", "ATCTTTAG"]);
        m.insert("SI-GA-C4", ["ACAATTCA", "TGCGCAGC", "CATCACTT", "GTGTGGAG"]);
        m.insert("SI-GA-C5", ["CGACTTGA", "TACAGACT", "ATTGCGTG", "GCGTACAC"]);
        m.insert("SI-GA-C6", ["ATTACTTC", "TGCGAACT", "GCATTCGG", "CAGCGGAA"]);
        m.insert("SI-GA-C7", ["GTCTCTCG", "AATCTCTC", "CGGAGGGA", "TCAGAAAT"]);
        m.insert("SI-GA-C8", ["GTTGAGAA", "AGATCTGG", "TCGATACT", "CACCGCTC"]);
        m.insert("SI-GA-C9", ["GCGCAGAA", "ATCTTACC", "TATGGTGT", "CGAACCTG"]);
        m.insert("SI-GA-D1", ["CACTCGGA", "GCTGAATT", "TGAAGTAC", "ATGCTCCG"]);
        m.insert(
            "SI-GA-D10",
            ["CAATACCC", "TGTCTATG", "ACCACGAA", "GTGGGTGT"],
        );
        m.insert(
            "SI-GA-D11",
            ["CTTTGCGG", "TGCACAAA", "AAGCAGTC", "GCAGTTCT"],
        );
        m.insert(
            "SI-GA-D12",
            ["GCACAATG", "CTTGGTAC", "TGCACCGT", "AAGTTGCA"],
        );
        m.insert("SI-GA-D2", ["TAACAAGG", "GGTTCCTC", "ATCATGCA", "CCGGGTAT"]);
        m.insert("SI-GA-D3", ["ACATTACT", "TTTGGGTA", "CAGCCCAC", "GGCAATGG"]);
        m.insert("SI-GA-D4", ["CCCTAACA", "ATTCCGAT", "TGGATTGC", "GAAGGCTG"]);
        m.insert("SI-GA-D5", ["CTCGTCAC", "GATCAGCA", "ACAACAGG", "TGGTGTTT"]);
        m.insert("SI-GA-D6", ["CATGCGAT", "TGATATTC", "GTGATCGA", "ACCCGACG"]);
        m.insert("SI-GA-D7", ["ATTTGCTA", "TAGACACC", "CCACAGGG", "GGCGTTAT"]);
        m.insert("SI-GA-D8", ["GCAACAAA", "TAGTTGTC", "CGCCATCG", "ATTGGCGT"]);
        m.insert("SI-GA-D9", ["AGGAGATG", "GATGTGGT", "CTACATCC", "TCCTCCAA"]);
        m.insert("SI-GA-E1", ["TGGTAAAC", "GAAAGGGT", "ACTGCTCG", "CTCCTCTA"]);
        m.insert(
            "SI-GA-E10",
            ["AAATGTGC", "GGGCAAAT", "TCTATCCG", "CTCGCGTA"],
        );
        m.insert(
            "SI-GA-E11",
            ["AAGCGCTG", "CGTTTGAT", "GTAGCACA", "TCCAATGC"],
        );
        m.insert(
            "SI-GA-E12",
            ["ACCGGCTC", "GAGTTAGT", "CGTCCTAG", "TTAAAGCA"],
        );
        m.insert("SI-GA-E2", ["GTGGTACC", "TACTATAG", "ACAAGGTA", "CGTCCCGT"]);
        m.insert("SI-GA-E3", ["AGGTATTG", "CTCCTAGT", "TCAAGGCC", "GATGCCAA"]);
        m.insert("SI-GA-E4", ["TTCGCCCT", "GGATGGGC", "AATCAATG", "CCGATTAA"]);
        m.insert("SI-GA-E5", ["CATTAGCG", "TTCGCTGA", "ACAAGAAT", "GGGCTCTC"]);
        m.insert("SI-GA-E6", ["CTGCGGCT", "GACTCAAA", "AGAAACTC", "TCTGTTGG"]);
        m.insert("SI-GA-E7", ["CACGCCTT", "GTATATAG", "TCTCGGGC", "AGGATACA"]);
        m.insert("SI-GA-E8", ["ATAGTTAC", "TGCTGAGT", "CCTACGTA", "GAGCACCG"]);
        m.insert("SI-GA-E9", ["TTGTTTCC", "GGAGGAGG", "CCTAACAA", "AACCCGTT"]);
        m.insert("SI-GA-F1", ["GTTGCAGC", "TGGAATTA", "CAATGGAG", "ACCCTCCT"]);
        m.insert(
            "SI-GA-F10",
            ["GCTTGGCT", "AAACAAAC", "CGGGCTTA", "TTCATCGG"],
        );
        m.insert(
            "SI-GA-F11",
            ["GCGAGAGT", "TACGTTCA", "AGTCCCAC", "CTATAGTG"],
        );
        m.insert(
            "SI-GA-F12",
            ["TGATGCAT", "GCTACTGA", "CACCTGCC", "ATGGAATG"],
        );
        m.insert("SI-GA-F2", ["TTTACATG", "CGCGATAC", "ACGCGGGT", "GAATTCCA"]);
        m.insert("SI-GA-F3", ["TTCAGGTG", "ACGGACAT", "GATCTTGA", "CGATCACC"]);
        m.insert("SI-GA-F4", ["CCCAATAG", "GTGTCGCT", "AGAGTCGC", "TATCGATA"]);
        m.insert("SI-GA-F5", ["GACTACGT", "CTAGCGAG", "TCTATATC", "AGGCGTCA"]);
        m.insert("SI-GA-F6", ["CGGAGCAC", "GACCTATT", "ACTTAGGA", "TTAGCTCG"]);
        m.insert("SI-GA-F7", ["CGTGCAGA", "AACAAGAT", "TCGCTTCG", "GTATGCTC"]);
        m.insert("SI-GA-F8", ["CATGAACA", "TCACTCGC", "AGCTGGAT", "GTGACTTG"]);
        m.insert("SI-GA-F9", ["CAAGCTCC", "GTTCACTG", "TCGTGAAA", "AGCATGGT"]);
        m.insert("SI-GA-G1", ["ATGAATCT", "GATCTCAG", "CCAGGAGC", "TGCTCGTA"]);
        m.insert(
            "SI-GA-G10",
            ["TCGCCAGC", "AATGTTAG", "CGATAGCT", "GTCAGCTA"],
        );
        m.insert(
            "SI-GA-G11",
            ["TTATCGTT", "AGCAGAGC", "CATCTCCA", "GCGGATAG"],
        );
        m.insert(
            "SI-GA-G12",
            ["ATTCTAAG", "CCCGATTA", "TGGAGGCT", "GAATCCGC"],
        );
        m.insert("SI-GA-G2", ["TGATTCTA", "ACTAGGAG", "CAGCCACT", "GTCGATGC"]);
        m.insert("SI-GA-G3", ["CCTCATTC", "AGCATCCG", "GTGGCAAT", "TAATGGGA"]);
        m.insert("SI-GA-G4", ["GCGATGTG", "AGATACAA", "TTTCCACT", "CACGGTGC"]);
        m.insert("SI-GA-G5", ["GAGCAAGA", "TCTGTGAT", "CGCAGTTC", "ATATCCCG"]);
        m.insert("SI-GA-G6", ["CTGACGCG", "GGTCGTAC", "TCCTTCTT", "AAAGAAGA"]);
        m.insert("SI-GA-G7", ["GGTATGCA", "CTCGAAAT", "ACACCTTC", "TAGTGCGG"]);
        m.insert("SI-GA-G8", ["TATGAGCT", "CCGATAGC", "ATACCCAA", "GGCTGTTG"]);
        m.insert("SI-GA-G9", ["TAGGACGT", "ATCCCACA", "GGAATGTC", "CCTTGTAG"]);
        m.insert("SI-GA-H1", ["GTATGTCA", "TGTCAGAC", "CACGTCGG", "ACGACATT"]);
        m.insert(
            "SI-GA-H10",
            ["GTAATTGC", "AGTCGCTT", "CACGAGAA", "TCGTCACG"],
        );
        m.insert(
            "SI-GA-H11",
            ["GGCGAGTA", "ACTTCTAT", "CAAATACG", "TTGCGCGC"],
        );
        m.insert(
            "SI-GA-H12",
            ["GACAGCAT", "TTTGTACA", "AGGCCGTG", "CCATATGC"],
        );
        m.insert("SI-GA-H2", ["TAATGACC", "ATGCCTTA", "GCCGAGAT", "CGTATCGG"]);
        m.insert("SI-GA-H3", ["CCAAGATG", "AGGCCCGA", "TACGTGAC", "GTTTATCT"]);
        m.insert("SI-GA-H4", ["GCCATTCC", "CAAGAATT", "TTGCCGGA", "AGTTGCAG"]);
        m.insert("SI-GA-H5", ["CCACTACA", "GATTCTGG", "TGCGGCTT", "ATGAAGAC"]);
        m.insert("SI-GA-H6", ["TAGGATAA", "CCTTTGTC", "GTACGCGG", "AGCACACT"]);
        m.insert("SI-GA-H7", ["AGCTATCA", "CATATAAC", "TCAGGGTG", "GTGCCCGT"]);
        m.insert("SI-GA-H8", ["TTGTTGAT", "GCTCAACC", "CAAAGTGG", "AGCGCCTA"]);
        m.insert("SI-GA-H9", ["ACACTGTT", "CAGGATGG", "GGCTGAAC", "TTTACCCA"]);
        m.insert(
            "SI-P01-A1",
            ["TTGTAAGA", "GGCGTTTC", "CCTACCAT", "AAACGGCG"],
        );
        m.insert(
            "SI-P01-A10",
            ["ACAGCAAC", "TTTCGCGA", "CGCAATTT", "GAGTTGCG"],
        );
        m.insert(
            "SI-P01-A11",
            ["CTTTCCTT", "TAGGTAAA", "ACCAGTCC", "GGACAGGG"],
        );
        m.insert(
            "SI-P01-A12",
            ["TTATGGAA", "ACTACTGT", "CGGGAACG", "GACCTCTC"],
        );
        m.insert(
            "SI-P01-A2",
            ["AGCCCTTT", "TCTTAGGC", "GTGAGAAG", "CAAGTCCA"],
        );
        m.insert(
            "SI-P01-A3",
            ["AAAGCATA", "GCCTTTAT", "CTGCAGCC", "TGTAGCGG"],
        );
        m.insert(
            "SI-P01-A4",
            ["CATGGCAG", "AGAACGCC", "GTCTTTGA", "TCGCAATT"],
        );
        m.insert(
            "SI-P01-A5",
            ["CAGTCTGG", "TCACACTC", "ATTGGGAA", "GGCATACT"],
        );
        m.insert(
            "SI-P01-A6",
            ["CTTTCGAC", "ACGGGACT", "TGCATCTG", "GAACATGA"],
        );
        m.insert(
            "SI-P01-A7",
            ["CTCATCAT", "TAACGTCC", "AGGTCATA", "GCTGAGGG"],
        );
        m.insert(
            "SI-P01-A8",
            ["GACTGTTC", "ATGATACG", "CCACAGAA", "TGTGCCGT"],
        );
        m.insert(
            "SI-P01-A9",
            ["GAGTGCGT", "CTCCAACA", "ACAACTTG", "TGTGTGAC"],
        );
        m.insert(
            "SI-P01-B1",
            ["CTAGCTGT", "GCCAACAA", "AGGCTACC", "TATTGGTG"],
        );
        m.insert(
            "SI-P01-B10",
            ["ACCATTAA", "CTGGACGT", "GAACGGTC", "TGTTCACG"],
        );
        m.insert(
            "SI-P01-B11",
            ["TCCAGATA", "GATTCGCT", "CGACATAG", "ATGGTCGC"],
        );
        m.insert(
            "SI-P01-B12",
            ["GAAAGACA", "CGCTACAT", "ACGCTTGG", "TTTGCGTC"],
        );
        m.insert(
            "SI-P01-B2",
            ["GGTCGAGC", "TTAGATTG", "CCCACCCA", "AAGTTGAT"],
        );
        m.insert(
            "SI-P01-B3",
            ["TCATCCTT", "ATTGGACG", "CAGCTTAC", "GGCAAGGA"],
        );
        m.insert(
            "SI-P01-B4",
            ["GACAGGCT", "CCTCTAAC", "AGGGACTG", "TTATCTGA"],
        );
        m.insert(
            "SI-P01-B5",
            ["TCATGCGA", "ATCGTACT", "CATCAGTG", "GGGACTAC"],
        );
        m.insert(
            "SI-P01-B6",
            ["GCGCACCT", "AACGCGAA", "CTATTTGG", "TGTAGATC"],
        );
        m.insert(
            "SI-P01-B7",
            ["TCCACACG", "CTTCTGTT", "GAATGCAC", "AGGGATGA"],
        );
        m.insert(
            "SI-P01-B8",
            ["ACGTTCAC", "TGCCCAGA", "CAAGGTCT", "GTTAAGTG"],
        );
        m.insert(
            "SI-P01-B9",
            ["AAGCGTGT", "CTTGACCG", "TGATTAAC", "GCCACGTA"],
        );
        m.insert(
            "SI-P01-C1",
            ["GATGCAGT", "AGACTTTC", "TTCTACAG", "CCGAGGCA"],
        );
        m.insert(
            "SI-P01-C10",
            ["CGTGCTAA", "TCACTCCT", "ATCTGATC", "GAGAAGGG"],
        );
        m.insert(
            "SI-P01-C11",
            ["GTTTGTGG", "ACCGAACA", "TAGACGAC", "CGACTCTT"],
        );
        m.insert(
            "SI-P01-C12",
            ["TAAGCCAC", "CCGTTATG", "GGTAATGT", "ATCCGGCA"],
        );
        m.insert(
            "SI-P01-C2",
            ["CCGAGAAC", "TGCTCTGT", "GTAGTGCG", "AATCACTA"],
        );
        m.insert(
            "SI-P01-C3",
            ["ACGTTACA", "TTACCTAC", "GACGACGG", "CGTAGGTT"],
        );
        m.insert(
            "SI-P01-C4",
            ["ACATTGGC", "GAGCCCAT", "CTTAGTCA", "TGCGAATG"],
        );
        m.insert(
            "SI-P01-C5",
            ["TCAGTCAA", "CACTGACT", "ATGCATTC", "GGTACGGG"],
        );
        m.insert(
            "SI-P01-C6",
            ["CGCTCAGG", "GAGGTTTA", "ACTCAGAC", "TTAAGCCT"],
        );
        m.insert(
            "SI-P01-C7",
            ["GGCGGAAT", "ACACCGGG", "CATAATCC", "TTGTTCTA"],
        );
        m.insert(
            "SI-P01-C8",
            ["TTTATCCC", "GCACGTTT", "CAGGAAGA", "AGCTCGAG"],
        );
        m.insert(
            "SI-P01-C9",
            ["AGATCGGT", "CATCGTCG", "GTCATATA", "TCGGACAC"],
        );
        m.insert(
            "SI-P01-D1",
            ["AGCTGCGT", "GTGGAGCA", "TCTATTAG", "CAACCATC"],
        );
        m.insert(
            "SI-P01-D10",
            ["CTTATTTG", "GCGGGCAT", "AGATAACA", "TACCCGGC"],
        );
        m.insert(
            "SI-P01-D11",
            ["GCTACTTC", "CACCTCAG", "ATATGAGA", "TGGGAGCT"],
        );
        m.insert(
            "SI-P01-D12",
            ["GCTGTGTA", "AGAAACGT", "CACTCAAC", "TTGCGTCG"],
        );
        m.insert(
            "SI-P01-D2",
            ["ACATTCCG", "GACACAAT", "CTGCGGTA", "TGTGATGC"],
        );
        m.insert(
            "SI-P01-D3",
            ["GAGCACGC", "CGAAGTTG", "TTCGTGAA", "ACTTCACT"],
        );
        m.insert(
            "SI-P01-D4",
            ["AAATCGTC", "GCGATCGG", "CTTCGAAT", "TGCGATCA"],
        );
        m.insert(
            "SI-P01-D5",
            ["GTCGACTC", "AGACGGAT", "CCTTTAGA", "TAGACTCG"],
        );
        m.insert(
            "SI-P01-D6",
            ["GAAGTCTT", "TGCAGGGC", "ATGCCAAA", "CCTTATCG"],
        );
        m.insert(
            "SI-P01-D7",
            ["CCGGATCC", "GGTCGCAT", "TTAACGTG", "AACTTAGA"],
        );
        m.insert(
            "SI-P01-D8",
            ["AATCTTTG", "GGATGAGT", "CTCAAGAC", "TCGGCCCA"],
        );
        m.insert(
            "SI-P01-D9",
            ["CAAGGGAC", "ACCTACTG", "GGGACACA", "TTTCTTGT"],
        );
        m.insert(
            "SI-P01-E1",
            ["CGCCCGTA", "GTTTGCCT", "TAAGTTAG", "ACGAAAGC"],
        );
        m.insert(
            "SI-P01-E10",
            ["GCACCAGT", "CGCAGGAG", "TAGTACCA", "ATTGTTTC"],
        );
        m.insert(
            "SI-P01-E11",
            ["ATCGCCAT", "TCACGGTA", "GGGTTTCC", "CATAAAGG"],
        );
        m.insert(
            "SI-P01-E12",
            ["CGCTATCC", "ACGCGGAA", "TAAATCGT", "GTTGCATG"],
        );
        m.insert(
            "SI-P01-E2",
            ["TCATCAAG", "GTTGGTCC", "AGGCTGGT", "CACAACTA"],
        );
        m.insert(
            "SI-P01-E3",
            ["TCTGCAGG", "CGGCTCCA", "AACAAGTC", "GTATGTAT"],
        );
        m.insert(
            "SI-P01-E4",
            ["GTAAACAT", "TATCCTGA", "AGCTGACG", "CCGGTGTC"],
        );
        m.insert(
            "SI-P01-E5",
            ["CCGTTGAA", "TATGCTCT", "ATCCAAGG", "GGAAGCTC"],
        );
        m.insert(
            "SI-P01-E6",
            ["TGATGGCT", "GCCATTTG", "ATTGAAAC", "CAGCCCGA"],
        );
        m.insert(
            "SI-P01-E7",
            ["AAGACGTG", "CCATGTGT", "GTTCACAA", "TGCGTACC"],
        );
        m.insert(
            "SI-P01-E8",
            ["GCTTACAT", "TAGGGTGC", "AGCCTATG", "CTAACGCA"],
        );
        m.insert(
            "SI-P01-E9",
            ["AGTAAGCA", "TACCGCGG", "CCGGTAAT", "GTATCTTC"],
        );
        m.insert(
            "SI-P01-F1",
            ["TGACTAGT", "GATAGGTA", "CCGTACAG", "ATCGCTCC"],
        );
        m.insert(
            "SI-P01-F10",
            ["TCGACAAT", "GAATACTG", "ATTCGTGC", "CGCGTGCA"],
        );
        m.insert(
            "SI-P01-F11",
            ["GAACCCGG", "AGCAGTTA", "TCGTAGAT", "CTTGTACC"],
        );
        m.insert(
            "SI-P01-F12",
            ["AATTGAAC", "TGGACCCT", "CCAGTGGA", "GTCCATTG"],
        );
        m.insert(
            "SI-P01-F2",
            ["TGCCCGCT", "GCAAACGC", "CATTGATA", "ATGGTTAG"],
        );
        m.insert(
            "SI-P01-F3",
            ["TTAGGACC", "AGTCTGTA", "GCCTCCGT", "CAGAATAG"],
        );
        m.insert(
            "SI-P01-F4",
            ["GCATGATA", "CGTCCTCT", "AACGACAC", "TTGATGGG"],
        );
        m.insert(
            "SI-P01-F5",
            ["TCTGACTA", "GTACGGCT", "CGGTTTAG", "AACACAGC"],
        );
        m.insert(
            "SI-P01-F6",
            ["GACTTCCT", "TGAGGAAG", "ATGCCGGC", "CCTAATTA"],
        );
        m.insert(
            "SI-P01-F7",
            ["GGTTAGAC", "CAAACTTT", "ACCCGAGA", "TTGGTCCG"],
        );
        m.insert(
            "SI-P01-F8",
            ["AGTTGGGA", "TACATTCT", "CCAGAAAG", "GTGCCCTC"],
        );
        m.insert(
            "SI-P01-F9",
            ["AGTTAGTT", "GTACTTAA", "CACGCACG", "TCGAGCGC"],
        );
        m.insert(
            "SI-P01-G1",
            ["CATATGCG", "ATGCGATT", "TCCGCTAC", "GGATACGA"],
        );
        m.insert(
            "SI-P01-G10",
            ["CGGAGACT", "TCCTATGA", "ATACTGAG", "GATGCCTC"],
        );
        m.insert(
            "SI-P01-G11",
            ["AGGGCGTT", "CTATACGC", "TACATAAG", "GCTCGTCA"],
        );
        m.insert(
            "SI-P01-G12",
            ["CATGCGTA", "ACCCGCAC", "TGATATCG", "GTGATAGT"],
        );
        m.insert(
            "SI-P01-G2",
            ["CGGTGAGC", "ATAACCTA", "TCCGAGCG", "GATCTTAT"],
        );
        m.insert(
            "SI-P01-G3",
            ["TACGAGTT", "ATGTCCAG", "GCTATAGC", "CGACGTCA"],
        );
        m.insert(
            "SI-P01-G4",
            ["CCTGCGGT", "GTACAACG", "AGCTTCTC", "TAGAGTAA"],
        );
        m.insert(
            "SI-P01-G5",
            ["ATGAAGTA", "GAAGCTCG", "TCTTTCGT", "CGCCGAAC"],
        );
        m.insert(
            "SI-P01-G6",
            ["GTGCGACA", "TCAGTGTT", "AGCACTGG", "CATTACAC"],
        );
        m.insert(
            "SI-P01-G7",
            ["GCCGGTAA", "TGACTGCC", "ATTACCGG", "CAGTAATT"],
        );
        m.insert(
            "SI-P01-G8",
            ["AAGTACTC", "GGAACTCT", "TCCCTGAG", "CTTGGAGA"],
        );
        m.insert(
            "SI-P01-G9",
            ["TTGACTTC", "GCCGAAGT", "CAATGGCA", "AGTCTCAG"],
        );
        m.insert(
            "SI-P01-H1",
            ["TCTTGTCC", "CGGGAGTA", "GTCACAGG", "AAACTCAT"],
        );
        m.insert(
            "SI-P01-H10",
            ["GACCGCCA", "TCGAATTG", "ATTTCAGC", "CGAGTGAT"],
        );
        m.insert(
            "SI-P01-H11",
            ["TCTCGACT", "AGGATCGA", "CACGATTC", "GTATCGAG"],
        );
        m.insert(
            "SI-P01-H12",
            ["TGTGTATA", "GTCCAGGC", "CAATCCCT", "ACGAGTAG"],
        );
        m.insert(
            "SI-P01-H2",
            ["CCGAACTC", "AACGGTCA", "TTATTGGT", "GGTCCAAG"],
        );
        m.insert(
            "SI-P01-H3",
            ["TTGGGCTT", "GACAAACC", "ACACCTAA", "CGTTTGGG"],
        );
        m.insert(
            "SI-P01-H4",
            ["TTCCATCT", "ACTGGAGC", "CGGTCGTG", "GAAATCAA"],
        );
        m.insert(
            "SI-P01-H5",
            ["ATAGTATG", "TATAAGGA", "GGCTCCAC", "CCGCGTCT"],
        );
        m.insert(
            "SI-P01-H6",
            ["AAGCATAA", "CCCATCGC", "TTAGCGCT", "GGTTGATG"],
        );
        m.insert(
            "SI-P01-H7",
            ["TGGCACGA", "AACGGGTG", "CTAATTCT", "GCTTCAAC"],
        );
        m.insert(
            "SI-P01-H8",
            ["AAGAGCGG", "TCATAGCA", "GGCCCATC", "CTTGTTAT"],
        );
        m.insert(
            "SI-P01-H9",
            ["GGAATATG", "ACCTGCCA", "CTTCATAC", "TAGGCGGT"],
        );
        m.insert(
            "SI-P02-A1",
            ["GGTTTACT", "CTAAACGG", "TCGGCGTC", "AACCGTAA"],
        );
        m.insert(
            "SI-P02-A10",
            ["GAAACCCT", "TTTCTGTC", "CCGTGTGA", "AGCGAAAG"],
        );
        m.insert(
            "SI-P02-A11",
            ["GTCCGGTC", "AAGATCAT", "CCTGAAGG", "TGATCTCA"],
        );
        m.insert(
            "SI-P02-A12",
            ["AGTGGAAC", "GTCTCCTT", "TCACATCA", "CAGATGGG"],
        );
        m.insert(
            "SI-P02-A2",
            ["TTTCATGA", "ACGTCCCT", "CGCATGTG", "GAAGGAAC"],
        );
        m.insert(
            "SI-P02-A3",
            ["CAGTACTG", "AGTAGTCT", "GCAGTAGA", "TTCCCGAC"],
        );
        m.insert(
            "SI-P02-A4",
            ["TATGATTC", "CCCACAGT", "ATGCTGAA", "GGATGCCG"],
        );
        m.insert(
            "SI-P02-A5",
            ["CTAGGTGA", "TCGTTCAG", "AGCCAATT", "GATACGCC"],
        );
        m.insert(
            "SI-P02-A6",
            ["CGCTATGT", "GCTGTCCA", "TTGAGATC", "AAACCGAG"],
        );
        m.insert(
            "SI-P02-A7",
            ["ACAGAGGT", "TATAGTTG", "CGGTCCCA", "GTCCTAAC"],
        );
        m.insert(
            "SI-P02-A8",
            ["GCATCTCC", "TGTAAGGT", "CTGCGATG", "AACGTCAA"],
        );
        m.insert(
            "SI-P02-A9",
            ["TCTTAAAG", "CGAGGCTC", "GTCCTTCT", "AAGACGGA"],
        );
        m.insert(
            "SI-P02-B1",
            ["GTAATCTT", "TCCGGAAG", "AGTTCGGC", "CAGCATCA"],
        );
        m.insert(
            "SI-P02-B10",
            ["ACCGTATG", "GATTAGAT", "CTGACTGA", "TGACGCCC"],
        );
        m.insert(
            "SI-P02-B11",
            ["GTTCCTCA", "AGGTACGC", "TAAGTATG", "CCCAGGAT"],
        );
        m.insert(
            "SI-P02-B12",
            ["TACCACCA", "CTAAGTTT", "GGGTCAAG", "ACTGTGGC"],
        );
        m.insert(
            "SI-P02-B2",
            ["TACTCTTC", "CCTGTGCG", "GGACACGT", "ATGAGAAA"],
        );
        m.insert(
            "SI-P02-B3",
            ["GTGTATTA", "TGTGCGGG", "ACCATAAC", "CAACGCCT"],
        );
        m.insert(
            "SI-P02-B4",
            ["ACTTCATA", "GAGATGAC", "TGCCGTGG", "CTAGACCT"],
        );
        m.insert(
            "SI-P02-B5",
            ["AATAATGG", "CCAGGGCA", "TGCCTCAT", "GTGTCATC"],
        );
        m.insert(
            "SI-P02-B6",
            ["CGTTAATC", "GCCACGCT", "TTACTCAG", "AAGGGTGA"],
        );
        m.insert(
            "SI-P02-B7",
            ["AAACCTCA", "GCCTTGGT", "CTGGACTC", "TGTAGAAG"],
        );
        m.insert(
            "SI-P02-B8",
            ["AAAGTGCT", "GCTACCTG", "TGCTGTAA", "CTGCAAGC"],
        );
        m.insert(
            "SI-P02-B9",
            ["CTGTAACT", "TCTAGCGA", "AGAGTGTG", "GACCCTAC"],
        );
        m.insert(
            "SI-P02-C1",
            ["CCACTTAT", "AACTGGCG", "TTGGCATA", "GGTAACGC"],
        );
        m.insert(
            "SI-P02-C10",
            ["TCTCAGTG", "GAGACTAT", "CGCTTAGC", "ATAGGCCA"],
        );
        m.insert(
            "SI-P02-C11",
            ["GAGGATCT", "AGACCATA", "TCCTGCGC", "CTTATGAG"],
        );
        m.insert(
            "SI-P02-C12",
            ["TCTCGTTT", "GGCTAGCG", "ATGACCGC", "CAAGTAAA"],
        );
        m.insert(
            "SI-P02-C2",
            ["CCTAGACC", "ATCTCTGT", "TAGCTCTA", "GGAGAGAG"],
        );
        m.insert(
            "SI-P02-C3",
            ["TCAGCCGT", "CAGAGGCC", "GGTCAATA", "ATCTTTAG"],
        );
        m.insert(
            "SI-P02-C4",
            ["ACAATTCA", "TGCGCAGC", "CATCACTT", "GTGTGGAG"],
        );
        m.insert(
            "SI-P02-C5",
            ["CGACTTGA", "TACAGACT", "ATTGCGTG", "GCGTACAC"],
        );
        m.insert(
            "SI-P02-C6",
            ["ATTACTTC", "TGCGAACT", "GCATTCGG", "CAGCGGAA"],
        );
        m.insert(
            "SI-P02-C7",
            ["GTCTCTCG", "AATCTCTC", "CGGAGGGA", "TCAGAAAT"],
        );
        m.insert(
            "SI-P02-C8",
            ["GTTGAGAA", "AGATCTGG", "TCGATACT", "CACCGCTC"],
        );
        m.insert(
            "SI-P02-C9",
            ["GCGCAGAA", "ATCTTACC", "TATGGTGT", "CGAACCTG"],
        );
        m.insert(
            "SI-P02-D1",
            ["CACTCGGA", "GCTGAATT", "TGAAGTAC", "ATGCTCCG"],
        );
        m.insert(
            "SI-P02-D10",
            ["CAATACCC", "TGTCTATG", "ACCACGAA", "GTGGGTGT"],
        );
        m.insert(
            "SI-P02-D11",
            ["CTTTGCGG", "TGCACAAA", "AAGCAGTC", "GCAGTTCT"],
        );
        m.insert(
            "SI-P02-D12",
            ["GCACAATG", "CTTGGTAC", "TGCACCGT", "AAGTTGCA"],
        );
        m.insert(
            "SI-P02-D2",
            ["TAACAAGG", "GGTTCCTC", "ATCATGCA", "CCGGGTAT"],
        );
        m.insert(
            "SI-P02-D3",
            ["ACATTACT", "TTTGGGTA", "CAGCCCAC", "GGCAATGG"],
        );
        m.insert(
            "SI-P02-D4",
            ["CCCTAACA", "ATTCCGAT", "TGGATTGC", "GAAGGCTG"],
        );
        m.insert(
            "SI-P02-D5",
            ["CTCGTCAC", "GATCAGCA", "ACAACAGG", "TGGTGTTT"],
        );
        m.insert(
            "SI-P02-D6",
            ["CATGCGAT", "TGATATTC", "GTGATCGA", "ACCCGACG"],
        );
        m.insert(
            "SI-P02-D7",
            ["ATTTGCTA", "TAGACACC", "CCACAGGG", "GGCGTTAT"],
        );
        m.insert(
            "SI-P02-D8",
            ["GCAACAAA", "TAGTTGTC", "CGCCATCG", "ATTGGCGT"],
        );
        m.insert(
            "SI-P02-D9",
            ["AGGAGATG", "GATGTGGT", "CTACATCC", "TCCTCCAA"],
        );
        m.insert(
            "SI-P02-E1",
            ["TGGTAAAC", "GAAAGGGT", "ACTGCTCG", "CTCCTCTA"],
        );
        m.insert(
            "SI-P02-E10",
            ["AAATGTGC", "GGGCAAAT", "TCTATCCG", "CTCGCGTA"],
        );
        m.insert(
            "SI-P02-E11",
            ["AAGCGCTG", "CGTTTGAT", "GTAGCACA", "TCCAATGC"],
        );
        m.insert(
            "SI-P02-E12",
            ["ACCGGCTC", "GAGTTAGT", "CGTCCTAG", "TTAAAGCA"],
        );
        m.insert(
            "SI-P02-E2",
            ["GTGGTACC", "TACTATAG", "ACAAGGTA", "CGTCCCGT"],
        );
        m.insert(
            "SI-P02-E3",
            ["AGGTATTG", "CTCCTAGT", "TCAAGGCC", "GATGCCAA"],
        );
        m.insert(
            "SI-P02-E4",
            ["TTCGCCCT", "GGATGGGC", "AATCAATG", "CCGATTAA"],
        );
        m.insert(
            "SI-P02-E5",
            ["CATTAGCG", "TTCGCTGA", "ACAAGAAT", "GGGCTCTC"],
        );
        m.insert(
            "SI-P02-E6",
            ["CTGCGGCT", "GACTCAAA", "AGAAACTC", "TCTGTTGG"],
        );
        m.insert(
            "SI-P02-E7",
            ["CACGCCTT", "GTATATAG", "TCTCGGGC", "AGGATACA"],
        );
        m.insert(
            "SI-P02-E8",
            ["ATAGTTAC", "TGCTGAGT", "CCTACGTA", "GAGCACCG"],
        );
        m.insert(
            "SI-P02-E9",
            ["TTGTTTCC", "GGAGGAGG", "CCTAACAA", "AACCCGTT"],
        );
        m.insert(
            "SI-P02-F1",
            ["GTTGCAGC", "TGGAATTA", "CAATGGAG", "ACCCTCCT"],
        );
        m.insert(
            "SI-P02-F10",
            ["GCTTGGCT", "AAACAAAC", "CGGGCTTA", "TTCATCGG"],
        );
        m.insert(
            "SI-P02-F11",
            ["GCGAGAGT", "TACGTTCA", "AGTCCCAC", "CTATAGTG"],
        );
        m.insert(
            "SI-P02-F12",
            ["TGATGCAT", "GCTACTGA", "CACCTGCC", "ATGGAATG"],
        );
        m.insert(
            "SI-P02-F2",
            ["TTTACATG", "CGCGATAC", "ACGCGGGT", "GAATTCCA"],
        );
        m.insert(
            "SI-P02-F3",
            ["TTCAGGTG", "ACGGACAT", "GATCTTGA", "CGATCACC"],
        );
        m.insert(
            "SI-P02-F4",
            ["CCCAATAG", "GTGTCGCT", "AGAGTCGC", "TATCGATA"],
        );
        m.insert(
            "SI-P02-F5",
            ["GACTACGT", "CTAGCGAG", "TCTATATC", "AGGCGTCA"],
        );
        m.insert(
            "SI-P02-F6",
            ["CGGAGCAC", "GACCTATT", "ACTTAGGA", "TTAGCTCG"],
        );
        m.insert(
            "SI-P02-F7",
            ["CGTGCAGA", "AACAAGAT", "TCGCTTCG", "GTATGCTC"],
        );
        m.insert(
            "SI-P02-F8",
            ["CATGAACA", "TCACTCGC", "AGCTGGAT", "GTGACTTG"],
        );
        m.insert(
            "SI-P02-F9",
            ["CAAGCTCC", "GTTCACTG", "TCGTGAAA", "AGCATGGT"],
        );
        m.insert(
            "SI-P02-G1",
            ["ATGAATCT", "GATCTCAG", "CCAGGAGC", "TGCTCGTA"],
        );
        m.insert(
            "SI-P02-G10",
            ["TCGCCAGC", "AATGTTAG", "CGATAGCT", "GTCAGCTA"],
        );
        m.insert(
            "SI-P02-G11",
            ["TTATCGTT", "AGCAGAGC", "CATCTCCA", "GCGGATAG"],
        );
        m.insert(
            "SI-P02-G12",
            ["ATTCTAAG", "CCCGATTA", "TGGAGGCT", "GAATCCGC"],
        );
        m.insert(
            "SI-P02-G2",
            ["TGATTCTA", "ACTAGGAG", "CAGCCACT", "GTCGATGC"],
        );
        m.insert(
            "SI-P02-G3",
            ["CCTCATTC", "AGCATCCG", "GTGGCAAT", "TAATGGGA"],
        );
        m.insert(
            "SI-P02-G4",
            ["GCGATGTG", "AGATACAA", "TTTCCACT", "CACGGTGC"],
        );
        m.insert(
            "SI-P02-G5",
            ["GAGCAAGA", "TCTGTGAT", "CGCAGTTC", "ATATCCCG"],
        );
        m.insert(
            "SI-P02-G6",
            ["CTGACGCG", "GGTCGTAC", "TCCTTCTT", "AAAGAAGA"],
        );
        m.insert(
            "SI-P02-G7",
            ["GGTATGCA", "CTCGAAAT", "ACACCTTC", "TAGTGCGG"],
        );
        m.insert(
            "SI-P02-G8",
            ["TATGAGCT", "CCGATAGC", "ATACCCAA", "GGCTGTTG"],
        );
        m.insert(
            "SI-P02-G9",
            ["TAGGACGT", "ATCCCACA", "GGAATGTC", "CCTTGTAG"],
        );
        m.insert(
            "SI-P02-H1",
            ["GTATGTCA", "TGTCAGAC", "CACGTCGG", "ACGACATT"],
        );
        m.insert(
            "SI-P02-H10",
            ["GTAATTGC", "AGTCGCTT", "CACGAGAA", "TCGTCACG"],
        );
        m.insert(
            "SI-P02-H11",
            ["GGCGAGTA", "ACTTCTAT", "CAAATACG", "TTGCGCGC"],
        );
        m.insert(
            "SI-P02-H12",
            ["GACAGCAT", "TTTGTACA", "AGGCCGTG", "CCATATGC"],
        );
        m.insert(
            "SI-P02-H2",
            ["TAATGACC", "ATGCCTTA", "GCCGAGAT", "CGTATCGG"],
        );
        m.insert(
            "SI-P02-H3",
            ["CCAAGATG", "AGGCCCGA", "TACGTGAC", "GTTTATCT"],
        );
        m.insert(
            "SI-P02-H4",
            ["GCCATTCC", "CAAGAATT", "TTGCCGGA", "AGTTGCAG"],
        );
        m.insert(
            "SI-P02-H5",
            ["CCACTACA", "GATTCTGG", "TGCGGCTT", "ATGAAGAC"],
        );
        m.insert(
            "SI-P02-H6",
            ["TAGGATAA", "CCTTTGTC", "GTACGCGG", "AGCACACT"],
        );
        m.insert(
            "SI-P02-H7",
            ["AGCTATCA", "CATATAAC", "TCAGGGTG", "GTGCCCGT"],
        );
        m.insert(
            "SI-P02-H8",
            ["TTGTTGAT", "GCTCAACC", "CAAAGTGG", "AGCGCCTA"],
        );
        m.insert(
            "SI-P02-H9",
            ["ACACTGTT", "CAGGATGG", "GGCTGAAC", "TTTACCCA"],
        );
        m.insert(
            "SI-P03-A1",
            ["AAACGGCG", "CCTACCAT", "GGCGTTTC", "TTGTAAGA"],
        );
        m.insert(
            "SI-P03-A10",
            ["ACAGCAAC", "CGCAATTT", "GAGTTGCG", "TTTCGCGA"],
        );
        m.insert(
            "SI-P03-A11",
            ["ACCAGTCC", "CTTTCCTT", "GGACAGGG", "TAGGTAAA"],
        );
        m.insert(
            "SI-P03-A12",
            ["ACTACTGT", "CGGGAACG", "GACCTCTC", "TTATGGAA"],
        );
        m.insert(
            "SI-P03-A2",
            ["AGCCCTTT", "CAAGTCCA", "GTGAGAAG", "TCTTAGGC"],
        );
        m.insert(
            "SI-P03-A3",
            ["AAAGCATA", "CTGCAGCC", "GCCTTTAT", "TGTAGCGG"],
        );
        m.insert(
            "SI-P03-A4",
            ["AGAACGCC", "CATGGCAG", "GTCTTTGA", "TCGCAATT"],
        );
        m.insert(
            "SI-P03-A5",
            ["ATTGGGAA", "CAGTCTGG", "GGCATACT", "TCACACTC"],
        );
        m.insert(
            "SI-P03-A6",
            ["ACGGGACT", "CTTTCGAC", "GAACATGA", "TGCATCTG"],
        );
        m.insert(
            "SI-P03-A7",
            ["AGGTCATA", "CTCATCAT", "GCTGAGGG", "TAACGTCC"],
        );
        m.insert(
            "SI-P03-A8",
            ["ATGATACG", "CCACAGAA", "GACTGTTC", "TGTGCCGT"],
        );
        m.insert(
            "SI-P03-A9",
            ["ACAACTTG", "CTCCAACA", "GAGTGCGT", "TGTGTGAC"],
        );
        m.insert(
            "SI-P03-B1",
            ["AGGCTACC", "CTAGCTGT", "GCCAACAA", "TATTGGTG"],
        );
        m.insert(
            "SI-P03-B10",
            ["ACCATTAA", "CTGGACGT", "GAACGGTC", "TGTTCACG"],
        );
        m.insert(
            "SI-P03-B11",
            ["ATGGTCGC", "CGACATAG", "GATTCGCT", "TCCAGATA"],
        );
        m.insert(
            "SI-P03-B12",
            ["ACGCTTGG", "CGCTACAT", "GAAAGACA", "TTTGCGTC"],
        );
        m.insert(
            "SI-P03-B2",
            ["AAGTTGAT", "CCCACCCA", "GGTCGAGC", "TTAGATTG"],
        );
        m.insert(
            "SI-P03-B3",
            ["ATTGGACG", "CAGCTTAC", "GGCAAGGA", "TCATCCTT"],
        );
        m.insert(
            "SI-P03-B4",
            ["AGGGACTG", "CCTCTAAC", "GACAGGCT", "TTATCTGA"],
        );
        m.insert(
            "SI-P03-B5",
            ["ATCGTACT", "CATCAGTG", "GGGACTAC", "TCATGCGA"],
        );
        m.insert(
            "SI-P03-B6",
            ["AACGCGAA", "CTATTTGG", "GCGCACCT", "TGTAGATC"],
        );
        m.insert(
            "SI-P03-B7",
            ["AGGGATGA", "CTTCTGTT", "GAATGCAC", "TCCACACG"],
        );
        m.insert(
            "SI-P03-B8",
            ["ACGTTCAC", "CAAGGTCT", "GTTAAGTG", "TGCCCAGA"],
        );
        m.insert(
            "SI-P03-B9",
            ["AAGCGTGT", "CTTGACCG", "GCCACGTA", "TGATTAAC"],
        );
        m.insert(
            "SI-P03-C1",
            ["AGACTTTC", "CCGAGGCA", "GATGCAGT", "TTCTACAG"],
        );
        m.insert(
            "SI-P03-C10",
            ["ATCTGATC", "CGTGCTAA", "GAGAAGGG", "TCACTCCT"],
        );
        m.insert(
            "SI-P03-C11",
            ["ACCGAACA", "CGACTCTT", "GTTTGTGG", "TAGACGAC"],
        );
        m.insert(
            "SI-P03-C12",
            ["ATCCGGCA", "CCGTTATG", "GGTAATGT", "TAAGCCAC"],
        );
        m.insert(
            "SI-P03-C2",
            ["AATCACTA", "CCGAGAAC", "GTAGTGCG", "TGCTCTGT"],
        );
        m.insert(
            "SI-P03-C3",
            ["ACGTTACA", "CGTAGGTT", "GACGACGG", "TTACCTAC"],
        );
        m.insert(
            "SI-P03-C4",
            ["ACATTGGC", "CTTAGTCA", "GAGCCCAT", "TGCGAATG"],
        );
        m.insert(
            "SI-P03-C5",
            ["ATGCATTC", "CACTGACT", "GGTACGGG", "TCAGTCAA"],
        );
        m.insert(
            "SI-P03-C6",
            ["ACTCAGAC", "CGCTCAGG", "GAGGTTTA", "TTAAGCCT"],
        );
        m.insert(
            "SI-P03-C7",
            ["ACACCGGG", "CATAATCC", "GGCGGAAT", "TTGTTCTA"],
        );
        m.insert(
            "SI-P03-C8",
            ["AGCTCGAG", "CAGGAAGA", "GCACGTTT", "TTTATCCC"],
        );
        m.insert(
            "SI-P03-C9",
            ["AGATCGGT", "CATCGTCG", "GTCATATA", "TCGGACAC"],
        );
        m.insert(
            "SI-P03-D1",
            ["AGCTGCGT", "CAACCATC", "GTGGAGCA", "TCTATTAG"],
        );
        m.insert(
            "SI-P03-D10",
            ["AGATAACA", "CTTATTTG", "GCGGGCAT", "TACCCGGC"],
        );
        m.insert(
            "SI-P03-D11",
            ["ATATGAGA", "CACCTCAG", "GCTACTTC", "TGGGAGCT"],
        );
        m.insert(
            "SI-P03-D12",
            ["AGAAACGT", "CACTCAAC", "GCTGTGTA", "TTGCGTCG"],
        );
        m.insert(
            "SI-P03-D2",
            ["ACATTCCG", "CTGCGGTA", "GACACAAT", "TGTGATGC"],
        );
        m.insert(
            "SI-P03-D3",
            ["ACTTCACT", "CGAAGTTG", "GAGCACGC", "TTCGTGAA"],
        );
        m.insert(
            "SI-P03-D4",
            ["AAATCGTC", "CTTCGAAT", "GCGATCGG", "TGCGATCA"],
        );
        m.insert(
            "SI-P03-D5",
            ["AGACGGAT", "CCTTTAGA", "GTCGACTC", "TAGACTCG"],
        );
        m.insert(
            "SI-P03-D6",
            ["ATGCCAAA", "CCTTATCG", "GAAGTCTT", "TGCAGGGC"],
        );
        m.insert(
            "SI-P03-D7",
            ["AACTTAGA", "CCGGATCC", "GGTCGCAT", "TTAACGTG"],
        );
        m.insert(
            "SI-P03-D8",
            ["AATCTTTG", "CTCAAGAC", "GGATGAGT", "TCGGCCCA"],
        );
        m.insert(
            "SI-P03-D9",
            ["ACCTACTG", "CAAGGGAC", "GGGACACA", "TTTCTTGT"],
        );
        m.insert(
            "SI-P03-E1",
            ["ACGAAAGC", "CGCCCGTA", "GTTTGCCT", "TAAGTTAG"],
        );
        m.insert(
            "SI-P03-E10",
            ["ATTGTTTC", "CGCAGGAG", "GCACCAGT", "TAGTACCA"],
        );
        m.insert(
            "SI-P03-E11",
            ["ATCGCCAT", "CATAAAGG", "GGGTTTCC", "TCACGGTA"],
        );
        m.insert(
            "SI-P03-E12",
            ["ACGCGGAA", "CGCTATCC", "GTTGCATG", "TAAATCGT"],
        );
        m.insert(
            "SI-P03-E2",
            ["AGGCTGGT", "CACAACTA", "GTTGGTCC", "TCATCAAG"],
        );
        m.insert(
            "SI-P03-E3",
            ["AACAAGTC", "CGGCTCCA", "GTATGTAT", "TCTGCAGG"],
        );
        m.insert(
            "SI-P03-E4",
            ["AGCTGACG", "CCGGTGTC", "GTAAACAT", "TATCCTGA"],
        );
        m.insert(
            "SI-P03-E5",
            ["ATCCAAGG", "CCGTTGAA", "GGAAGCTC", "TATGCTCT"],
        );
        m.insert(
            "SI-P03-E6",
            ["ATTGAAAC", "CAGCCCGA", "GCCATTTG", "TGATGGCT"],
        );
        m.insert(
            "SI-P03-E7",
            ["AAGACGTG", "CCATGTGT", "GTTCACAA", "TGCGTACC"],
        );
        m.insert(
            "SI-P03-E8",
            ["AGCCTATG", "CTAACGCA", "GCTTACAT", "TAGGGTGC"],
        );
        m.insert(
            "SI-P03-E9",
            ["AGTAAGCA", "CCGGTAAT", "GTATCTTC", "TACCGCGG"],
        );
        m.insert(
            "SI-P03-F1",
            ["ATCGCTCC", "CCGTACAG", "GATAGGTA", "TGACTAGT"],
        );
        m.insert(
            "SI-P03-F10",
            ["ATTCGTGC", "CGCGTGCA", "GAATACTG", "TCGACAAT"],
        );
        m.insert(
            "SI-P03-F11",
            ["AGCAGTTA", "CTTGTACC", "GAACCCGG", "TCGTAGAT"],
        );
        m.insert(
            "SI-P03-F12",
            ["AATTGAAC", "CCAGTGGA", "GTCCATTG", "TGGACCCT"],
        );
        m.insert(
            "SI-P03-F2",
            ["ATGGTTAG", "CATTGATA", "GCAAACGC", "TGCCCGCT"],
        );
        m.insert(
            "SI-P03-F3",
            ["AGTCTGTA", "CAGAATAG", "GCCTCCGT", "TTAGGACC"],
        );
        m.insert(
            "SI-P03-F4",
            ["AACGACAC", "CGTCCTCT", "GCATGATA", "TTGATGGG"],
        );
        m.insert(
            "SI-P03-F5",
            ["AACACAGC", "CGGTTTAG", "GTACGGCT", "TCTGACTA"],
        );
        m.insert(
            "SI-P03-F6",
            ["ATGCCGGC", "CCTAATTA", "GACTTCCT", "TGAGGAAG"],
        );
        m.insert(
            "SI-P03-F7",
            ["ACCCGAGA", "CAAACTTT", "GGTTAGAC", "TTGGTCCG"],
        );
        m.insert(
            "SI-P03-F8",
            ["AGTTGGGA", "CCAGAAAG", "GTGCCCTC", "TACATTCT"],
        );
        m.insert(
            "SI-P03-F9",
            ["AGTTAGTT", "CACGCACG", "GTACTTAA", "TCGAGCGC"],
        );
        m.insert(
            "SI-P03-G1",
            ["ATGCGATT", "CATATGCG", "GGATACGA", "TCCGCTAC"],
        );
        m.insert(
            "SI-P03-G10",
            ["ATACTGAG", "CGGAGACT", "GATGCCTC", "TCCTATGA"],
        );
        m.insert(
            "SI-P03-G11",
            ["AGGGCGTT", "CTATACGC", "GCTCGTCA", "TACATAAG"],
        );
        m.insert(
            "SI-P03-G12",
            ["ACCCGCAC", "CATGCGTA", "GTGATAGT", "TGATATCG"],
        );
        m.insert(
            "SI-P03-G2",
            ["ATAACCTA", "CGGTGAGC", "GATCTTAT", "TCCGAGCG"],
        );
        m.insert(
            "SI-P03-G3",
            ["ATGTCCAG", "CGACGTCA", "GCTATAGC", "TACGAGTT"],
        );
        m.insert(
            "SI-P03-G4",
            ["AGCTTCTC", "CCTGCGGT", "GTACAACG", "TAGAGTAA"],
        );
        m.insert(
            "SI-P03-G5",
            ["ATGAAGTA", "CGCCGAAC", "GAAGCTCG", "TCTTTCGT"],
        );
        m.insert(
            "SI-P03-G6",
            ["AGCACTGG", "CATTACAC", "GTGCGACA", "TCAGTGTT"],
        );
        m.insert(
            "SI-P03-G7",
            ["ATTACCGG", "CAGTAATT", "GCCGGTAA", "TGACTGCC"],
        );
        m.insert(
            "SI-P03-G8",
            ["AAGTACTC", "CTTGGAGA", "GGAACTCT", "TCCCTGAG"],
        );
        m.insert(
            "SI-P03-G9",
            ["AGTCTCAG", "CAATGGCA", "GCCGAAGT", "TTGACTTC"],
        );
        m.insert(
            "SI-P03-H1",
            ["AAACTCAT", "CGGGAGTA", "GTCACAGG", "TCTTGTCC"],
        );
        m.insert(
            "SI-P03-H10",
            ["ATTTCAGC", "CGAGTGAT", "GACCGCCA", "TCGAATTG"],
        );
        m.insert(
            "SI-P03-H11",
            ["AGGATCGA", "CACGATTC", "GTATCGAG", "TCTCGACT"],
        );
        m.insert(
            "SI-P03-H12",
            ["ACGAGTAG", "CAATCCCT", "GTCCAGGC", "TGTGTATA"],
        );
        m.insert(
            "SI-P03-H2",
            ["AACGGTCA", "CCGAACTC", "GGTCCAAG", "TTATTGGT"],
        );
        m.insert(
            "SI-P03-H3",
            ["ACACCTAA", "CGTTTGGG", "GACAAACC", "TTGGGCTT"],
        );
        m.insert(
            "SI-P03-H4",
            ["ACTGGAGC", "CGGTCGTG", "GAAATCAA", "TTCCATCT"],
        );
        m.insert(
            "SI-P03-H5",
            ["ATAGTATG", "CCGCGTCT", "GGCTCCAC", "TATAAGGA"],
        );
        m.insert(
            "SI-P03-H6",
            ["AAGCATAA", "CCCATCGC", "GGTTGATG", "TTAGCGCT"],
        );
        m.insert(
            "SI-P03-H7",
            ["AACGGGTG", "CTAATTCT", "GCTTCAAC", "TGGCACGA"],
        );
        m.insert(
            "SI-P03-H8",
            ["AAGAGCGG", "CTTGTTAT", "GGCCCATC", "TCATAGCA"],
        );
        m.insert(
            "SI-P03-H9",
            ["ACCTGCCA", "CTTCATAC", "GGAATATG", "TAGGCGGT"],
        );
        m.insert("SI-P2-A1", ["GGTTTACT", "CTAAACGG", "TCGGCGTC", "AACCGTAA"]);
        m.insert(
            "SI-P2-A10",
            ["GAAACCCT", "TTTCTGTC", "CCGTGTGA", "AGCGAAAG"],
        );
        m.insert(
            "SI-P2-A11",
            ["GTCCGGTC", "AAGATCAT", "CCTGAAGG", "TGATCTCA"],
        );
        m.insert(
            "SI-P2-A12",
            ["AGTGGAAC", "GTCTCCTT", "TCACATCA", "CAGATGGG"],
        );
        m.insert("SI-P2-A2", ["TTTCATGA", "ACGTCCCT", "CGCATGTG", "GAAGGAAC"]);
        m.insert("SI-P2-A3", ["CAGTACTG", "AGTAGTCT", "GCAGTAGA", "TTCCCGAC"]);
        m.insert("SI-P2-A4", ["TATGATTC", "CCCACAGT", "ATGCTGAA", "GGATGCCG"]);
        m.insert("SI-P2-A5", ["CTAGGTGA", "TCGTTCAG", "AGCCAATT", "GATACGCC"]);
        m.insert("SI-P2-A6", ["CGCTATGT", "GCTGTCCA", "TTGAGATC", "AAACCGAG"]);
        m.insert("SI-P2-A7", ["ACAGAGGT", "TATAGTTG", "CGGTCCCA", "GTCCTAAC"]);
        m.insert("SI-P2-A8", ["GCATCTCC", "TGTAAGGT", "CTGCGATG", "AACGTCAA"]);
        m.insert("SI-P2-A9", ["TCTTAAAG", "CGAGGCTC", "GTCCTTCT", "AAGACGGA"]);
        m.insert("SI-P2-B1", ["GTAATCTT", "TCCGGAAG", "AGTTCGGC", "CAGCATCA"]);
        m.insert(
            "SI-P2-B10",
            ["ACCGTATG", "GATTAGAT", "CTGACTGA", "TGACGCCC"],
        );
        m.insert(
            "SI-P2-B11",
            ["GTTCCTCA", "AGGTACGC", "TAAGTATG", "CCCAGGAT"],
        );
        m.insert(
            "SI-P2-B12",
            ["TACCACCA", "CTAAGTTT", "GGGTCAAG", "ACTGTGGC"],
        );
        m.insert("SI-P2-B2", ["TACTCTTC", "CCTGTGCG", "GGACACGT", "ATGAGAAA"]);
        m.insert("SI-P2-B3", ["GTGTATTA", "TGTGCGGG", "ACCATAAC", "CAACGCCT"]);
        m.insert("SI-P2-B4", ["ACTTCATA", "GAGATGAC", "TGCCGTGG", "CTAGACCT"]);
        m.insert("SI-P2-B5", ["AATAATGG", "CCAGGGCA", "TGCCTCAT", "GTGTCATC"]);
        m.insert("SI-P2-B6", ["CGTTAATC", "GCCACGCT", "TTACTCAG", "AAGGGTGA"]);
        m.insert("SI-P2-B7", ["AAACCTCA", "GCCTTGGT", "CTGGACTC", "TGTAGAAG"]);
        m.insert("SI-P2-B8", ["AAAGTGCT", "GCTACCTG", "TGCTGTAA", "CTGCAAGC"]);
        m.insert("SI-P2-B9", ["CTGTAACT", "TCTAGCGA", "AGAGTGTG", "GACCCTAC"]);
        m.insert("SI-P2-C1", ["CCACTTAT", "AACTGGCG", "TTGGCATA", "GGTAACGC"]);
        m.insert(
            "SI-P2-C10",
            ["TCTCAGTG", "GAGACTAT", "CGCTTAGC", "ATAGGCCA"],
        );
        m.insert(
            "SI-P2-C11",
            ["GAGGATCT", "AGACCATA", "TCCTGCGC", "CTTATGAG"],
        );
        m.insert(
            "SI-P2-C12",
            ["TCTCGTTT", "GGCTAGCG", "ATGACCGC", "CAAGTAAA"],
        );
        m.insert("SI-P2-C2", ["CCTAGACC", "ATCTCTGT", "TAGCTCTA", "GGAGAGAG"]);
        m.insert("SI-P2-C3", ["TCAGCCGT", "CAGAGGCC", "GGTCAATA", "ATCTTTAG"]);
        m.insert("SI-P2-C4", ["ACAATTCA", "TGCGCAGC", "CATCACTT", "GTGTGGAG"]);
        m.insert("SI-P2-C5", ["CGACTTGA", "TACAGACT", "ATTGCGTG", "GCGTACAC"]);
        m.insert("SI-P2-C6", ["ATTACTTC", "TGCGAACT", "GCATTCGG", "CAGCGGAA"]);
        m.insert("SI-P2-C7", ["GTCTCTCG", "AATCTCTC", "CGGAGGGA", "TCAGAAAT"]);
        m.insert("SI-P2-C8", ["GTTGAGAA", "AGATCTGG", "TCGATACT", "CACCGCTC"]);
        m.insert("SI-P2-C9", ["GCGCAGAA", "ATCTTACC", "TATGGTGT", "CGAACCTG"]);
        m.insert("SI-P2-D1", ["CACTCGGA", "GCTGAATT", "TGAAGTAC", "ATGCTCCG"]);
        m.insert(
            "SI-P2-D10",
            ["CAATACCC", "TGTCTATG", "ACCACGAA", "GTGGGTGT"],
        );
        m.insert(
            "SI-P2-D11",
            ["CTTTGCGG", "TGCACAAA", "AAGCAGTC", "GCAGTTCT"],
        );
        m.insert(
            "SI-P2-D12",
            ["GCACAATG", "CTTGGTAC", "TGCACCGT", "AAGTTGCA"],
        );
        m.insert("SI-P2-D2", ["TAACAAGG", "GGTTCCTC", "ATCATGCA", "CCGGGTAT"]);
        m.insert("SI-P2-D3", ["ACATTACT", "TTTGGGTA", "CAGCCCAC", "GGCAATGG"]);
        m.insert("SI-P2-D4", ["CCCTAACA", "ATTCCGAT", "TGGATTGC", "GAAGGCTG"]);
        m.insert("SI-P2-D5", ["CTCGTCAC", "GATCAGCA", "ACAACAGG", "TGGTGTTT"]);
        m.insert("SI-P2-D6", ["CATGCGAT", "TGATATTC", "GTGATCGA", "ACCCGACG"]);
        m.insert("SI-P2-D7", ["ATTTGCTA", "TAGACACC", "CCACAGGG", "GGCGTTAT"]);
        m.insert("SI-P2-D8", ["GCAACAAA", "TAGTTGTC", "CGCCATCG", "ATTGGCGT"]);
        m.insert("SI-P2-D9", ["AGGAGATG", "GATGTGGT", "CTACATCC", "TCCTCCAA"]);
        m.insert("SI-P2-E1", ["TGGTAAAC", "GAAAGGGT", "ACTGCTCG", "CTCCTCTA"]);
        m.insert(
            "SI-P2-E10",
            ["AAATGTGC", "GGGCAAAT", "TCTATCCG", "CTCGCGTA"],
        );
        m.insert(
            "SI-P2-E11",
            ["AAGCGCTG", "CGTTTGAT", "GTAGCACA", "TCCAATGC"],
        );
        m.insert(
            "SI-P2-E12",
            ["ACCGGCTC", "GAGTTAGT", "CGTCCTAG", "TTAAAGCA"],
        );
        m.insert("SI-P2-E2", ["GTGGTACC", "TACTATAG", "ACAAGGTA", "CGTCCCGT"]);
        m.insert("SI-P2-E3", ["AGGTATTG", "CTCCTAGT", "TCAAGGCC", "GATGCCAA"]);
        m.insert("SI-P2-E4", ["TTCGCCCT", "GGATGGGC", "AATCAATG", "CCGATTAA"]);
        m.insert("SI-P2-E5", ["CATTAGCG", "TTCGCTGA", "ACAAGAAT", "GGGCTCTC"]);
        m.insert("SI-P2-E6", ["CTGCGGCT", "GACTCAAA", "AGAAACTC", "TCTGTTGG"]);
        m.insert("SI-P2-E7", ["CACGCCTT", "GTATATAG", "TCTCGGGC", "AGGATACA"]);
        m.insert("SI-P2-E8", ["ATAGTTAC", "TGCTGAGT", "CCTACGTA", "GAGCACCG"]);
        m.insert("SI-P2-E9", ["TTGTTTCC", "GGAGGAGG", "CCTAACAA", "AACCCGTT"]);
        m.insert("SI-P2-F1", ["GTTGCAGC", "TGGAATTA", "CAATGGAG", "ACCCTCCT"]);
        m.insert(
            "SI-P2-F10",
            ["GCTTGGCT", "AAACAAAC", "CGGGCTTA", "TTCATCGG"],
        );
        m.insert(
            "SI-P2-F11",
            ["GCGAGAGT", "TACGTTCA", "AGTCCCAC", "CTATAGTG"],
        );
        m.insert(
            "SI-P2-F12",
            ["TGATGCAT", "GCTACTGA", "CACCTGCC", "ATGGAATG"],
        );
        m.insert("SI-P2-F2", ["TTTACATG", "CGCGATAC", "ACGCGGGT", "GAATTCCA"]);
        m.insert("SI-P2-F3", ["TTCAGGTG", "ACGGACAT", "GATCTTGA", "CGATCACC"]);
        m.insert("SI-P2-F4", ["CCCAATAG", "GTGTCGCT", "AGAGTCGC", "TATCGATA"]);
        m.insert("SI-P2-F5", ["GACTACGT", "CTAGCGAG", "TCTATATC", "AGGCGTCA"]);
        m.insert("SI-P2-F6", ["CGGAGCAC", "GACCTATT", "ACTTAGGA", "TTAGCTCG"]);
        m.insert("SI-P2-F7", ["CGTGCAGA", "AACAAGAT", "TCGCTTCG", "GTATGCTC"]);
        m.insert("SI-P2-F8", ["CATGAACA", "TCACTCGC", "AGCTGGAT", "GTGACTTG"]);
        m.insert("SI-P2-F9", ["CAAGCTCC", "GTTCACTG", "TCGTGAAA", "AGCATGGT"]);
        m.insert("SI-P2-G1", ["ATGAATCT", "GATCTCAG", "CCAGGAGC", "TGCTCGTA"]);
        m.insert(
            "SI-P2-G10",
            ["TCGCCAGC", "AATGTTAG", "CGATAGCT", "GTCAGCTA"],
        );
        m.insert(
            "SI-P2-G11",
            ["TTATCGTT", "AGCAGAGC", "CATCTCCA", "GCGGATAG"],
        );
        m.insert(
            "SI-P2-G12",
            ["ATTCTAAG", "CCCGATTA", "TGGAGGCT", "GAATCCGC"],
        );
        m.insert("SI-P2-G2", ["TGATTCTA", "ACTAGGAG", "CAGCCACT", "GTCGATGC"]);
        m.insert("SI-P2-G3", ["CCTCATTC", "AGCATCCG", "GTGGCAAT", "TAATGGGA"]);
        m.insert("SI-P2-G4", ["GCGATGTG", "AGATACAA", "TTTCCACT", "CACGGTGC"]);
        m.insert("SI-P2-G5", ["GAGCAAGA", "TCTGTGAT", "CGCAGTTC", "ATATCCCG"]);
        m.insert("SI-P2-G6", ["CTGACGCG", "GGTCGTAC", "TCCTTCTT", "AAAGAAGA"]);
        m.insert("SI-P2-G7", ["GGTATGCA", "CTCGAAAT", "ACACCTTC", "TAGTGCGG"]);
        m.insert("SI-P2-G8", ["TATGAGCT", "CCGATAGC", "ATACCCAA", "GGCTGTTG"]);
        m.insert("SI-P2-G9", ["TAGGACGT", "ATCCCACA", "GGAATGTC", "CCTTGTAG"]);
        m.insert("SI-P2-H1", ["GTATGTCA", "TGTCAGAC", "CACGTCGG", "ACGACATT"]);
        m.insert(
            "SI-P2-H10",
            ["GTAATTGC", "AGTCGCTT", "CACGAGAA", "TCGTCACG"],
        );
        m.insert(
            "SI-P2-H11",
            ["GGCGAGTA", "ACTTCTAT", "CAAATACG", "TTGCGCGC"],
        );
        m.insert(
            "SI-P2-H12",
            ["GACAGCAT", "TTTGTACA", "AGGCCGTG", "CCATATGC"],
        );
        m.insert("SI-P2-H2", ["TAATGACC", "ATGCCTTA", "GCCGAGAT", "CGTATCGG"]);
        m.insert("SI-P2-H3", ["CCAAGATG", "AGGCCCGA", "TACGTGAC", "GTTTATCT"]);
        m.insert("SI-P2-H4", ["GCCATTCC", "CAAGAATT", "TTGCCGGA", "AGTTGCAG"]);
        m.insert("SI-P2-H5", ["CCACTACA", "GATTCTGG", "TGCGGCTT", "ATGAAGAC"]);
        m.insert("SI-P2-H6", ["TAGGATAA", "CCTTTGTC", "GTACGCGG", "AGCACACT"]);
        m.insert("SI-P2-H7", ["AGCTATCA", "CATATAAC", "TCAGGGTG", "GTGCCCGT"]);
        m.insert("SI-P2-H8", ["TTGTTGAT", "GCTCAACC", "CAAAGTGG", "AGCGCCTA"]);
        m.insert("SI-P2-H9", ["ACACTGTT", "CAGGATGG", "GGCTGAAC", "TTTACCCA"]);
        m.insert("SI-T2-1", ["GGGTGATC", "TTACCGAT", "AATGACGA", "CCCATTCG"]);
        m.insert("SI-T2-2", ["GGGTCGAA", "ATCCGCCC", "TCTATAGT", "CAAGATTG"]);
        m.insert("SI-T2-3", ["GCTGATAT", "TGCCGAGC", "AAATTGCG", "CTGACCTA"]);
        m.insert("SI-T2-4", ["ACTTCTGA", "TTCATCTT", "CGACGACG", "GAGGAGAC"]);
        m.insert("SI-T2-5", ["GAATACAA", "AGCATACC", "TCGGGTTT", "CTTCCGGG"]);
        m.insert("SI-T2-6", ["TATTGAGA", "GTAGTCAG", "CGCCATTC", "ACGACGCT"]);
        m.insert("SI-T2-7", ["AAATCTGT", "GTCCAACC", "TCTGGCTG", "CGGATGAA"]);
        m.insert("SI-T2-8", ["CCTTGAAC", "GAAATCGG", "TGGCCTCT", "ATCGAGTA"]);
        m
    };
}
