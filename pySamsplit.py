import click
import pysam
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

PATTERN_DICT = {
    "auto": None,
    "scopeV2.0.0": "C8L16C8L16C8U8T18",
    "scopeV2.0.1": "C8L16C8L16C8L1U8T18",
    "scopeV2.1.0": "C8L16C8L16C8U12T18",
    "scopeV2.1.1": "C8L16C8L16C8L1U12T18",
    "scopeV2.2.1": "C8L16C8L16C8L1U12T18",
    "scopeV3.0.1": "C9L16C9L16C9L1U12T18",
    "flv_rna": "C8L16C8L16C8U9L6",
    "flv": "U9C8L16C8L16C8",
    "bulk_vdj": "L18C6U16",
    "bulk_rna": "C9U12",
    "customized": None,
    # "rna_5p": "U8C9L4C9L4C9",
    # "rna_3p": "C9L4C9L4C9U8",
    "scope_5p3p": "C9C9C9U8"
}


bclist_linker_map = {
    "scopeV2.0.0": (
        [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACG"
    ),
    "scopeV2.0.1":(
        [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACGC"
    ),
    "scopeV2.1.0":(
        [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACG"
    ),
    "scopeV2.1.1":(
        [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACGC"
    ),
    "scopeV2.2.1":(
        [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACGC", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATC", "GATTGTCACTAACGCGATGCTGACTCCTAGTCC"
    ),
    "scopeV3.0.1":(
        [
            "AACGGACCT", "AGGACTCAC", "ACTGCCTAG", "GAACGCTAT", "GACTGGTTG",
            "CAGGACTTC", "TCGGTTCGT", "GTCTTGCGT", "TATCTCCGA", "CGCAACTAC",
            "CACTTCAGA", "GCTCTCACT", "ACGAAGCTC", "GTGTTAAGG", "GGCTCTCTA",
            "GCGTAGTAA", "AGCTCCTTG", "CACATTCAC", "TCCGTATCA", "CACCTGTAA",
            "ATGGTCTCA", "TGCATCAAG", "CGCCAATGA", "TCGACTGTC", "TGTGGACAC",
            "ACATCGGAC", "GCTTGAGGT", "GCCGTTATT", "AACACCGTT", "CCTAGTCTT",
            "AACACACAG", "CAACATCGG", "TACCTCTCC", "CCAATGTCT", "CGAGATAGT",
            "AGTTCAGAG", "ATCGACACG", "GACCTTAGC", "TAACCTACC", "TGGAGAACC",
            "GATGTTACG", "TGGCATGAG", "CTGGTACTT", "GCGAGTAAC", "GATCCATGC",
            "GGCTTCTCA", "CGTTAGCGT", "ACACAGGCT", "AACTGGCGA", "AGACGTTCA",
            "TGCGGTTCT", "CAGTCTTCG", "AGCTGAGTC", "GTGCATATC", "GATGGCTCA",
            "TGAGCGAAG", "CTTGACGTT", "GTTACTGGT", "GAGCAGCTT", "GTTGGAGTG",
            "TGCCTGATC", "GCAGATGTG", "CAGAGTACA", "TGTAGTGTG", "ACGAATGGA",
            "TCACTGGAA", "ATGACAGCA", "CTCAGAACT", "AAGCTTGCG", "TCGGACATG",
            "GAGGATTGA", "TCTTGGACA", "GTGCAAGGT", "CGATCGGTA", "GGATACCAC",
            "ATCGTTGGC", "TGGAACGTA", "GCCTACGAT", "CATGTAGGC", "CTCACGTTC",
            "TGTGAGTCA", "AATCGCCAC", "GTGCGACTA", "TCTGGCGAA", "ATACGCGGA",
            "CCTTGAATC", "CCACACATT", "CGGTGATTG", "AGGAGCAAT", "ACACACCAA",
            "GGATAGATG", "TACCGTCTG", "TCTTGCTTG", "CCAGCTAAC", "CTTACGCAG",
            "AGTAGGAGG",
        
            "TTGGTGACC", "CATCGGTTC", "CAATGCAAC", "AAGGTGGTA", "TCAGGTAGA",
            "CACTAGGCA", "GGCATGCAA", "ATGCGATAC", "AGACGAAGT", "CGATAAGGC",
            "CGAGTTGCA", "GGTGATCAG", "AGTGGTGGT", "TAGCGATGA", "ATGCCTAAG",
            "CCATAATCG", "CTAACCAGA", "TGATGTGCC", "TAGCTACAC", "TCCTGGCTT",
            "GTAGTTCCT", "GTATCCTTC", "GCCATAACC", "GCATGTTGG", "AATCCGGTG",
            "ACGTTACGA", "AGCATAGCG", "GGTCCGTAA", "GTTCTACCG", "CGCTGTAAG",
            "CGTCATACC", "TGGTAACCG", "ACAACAGGT", "TGACTTCCG", "CACGCAATA",
            "ACCGTACTC", "GATGTGTGT", "ACACCAACG", "TTGAGACAG", "CGGATCATC",
            "TAATGGCCG", "ACCTCGACT", "CGGCTAGAT", "TGGACTTGT", "TGATCCTCT",
            "ATAGCGTGT", "TGTCGGTGT", "AGCCACATA", "GAAGAAGCC", "CTGTGGTAT",
            "CAGCCGAAT", "ACCTGCTAC", "TCTTCTCAG", "GATCAGGAC", "TAGACCACT",
            "TGAGTAGTC", "TTCGAGGAT", "ATGTATCGG", "TCTGTCTGC", "GCTTACAGG",
            "ATACCAGTC", "AACGTCCAA", "CTCCTCAAT", "TTCCAATCG", "AACGCTAGT",
            "TTACACGAC", "CATACGACC", "CTAATCGCG", "CGTAATTGG", "GTAGTGTTG",
            "TATAGCGGT", "CGTACTGAA", "TCGATGTGG", "CCAGAAGAT", "AGGCTGTTG",
            "TGAGGCCTT", "CCAGTCCTA", "CAATTGCGC", "TTGCCGTCA", "GAGTTGACA",
            "TGGTCAGTT", "GGTAGTCCA", "ATCCTTCCA", "ATGTGCAGC", "GGAAGACTC",
            "GTAATGGAC", "TATCGTGCA", "CCTACAAGG", "TAGTCCGGA", "AGATACGCA",
            "ACTCATCGT", "GATAACCGC", "GCTGCGATA", "TACTCACCA", "ACCAGGTCA",
            "GCAACTTCA",
        
            "TTCGGTCAA", "GTACGGACT", "GGCCATGTT", "TATGACACC", "GTGTCTGAA",
            "TAAGCTTGG", "AAGATCTGC", "CTTGTGCCA", "ACTGGTTCC", "GATTGCGAG",
            "ACACGGTAG", "TCGCATACT", "CCGATCGAA", "GTAGCACGA", "TCAGCACTG",
            "TGCTATCGC", "GTAAGATCG", "TGCACGAGA", "GAATCGCAA", "CCTCGATCA",
            "GCGTCTAGT", "TCCAGTTAG", "CAGGAACCA", "GTGAACCAA", "CTCTAACAC",
            "GTCTGATGC", "TTAAGAGGC", "GTCAATGCA", "CTTCATGGA", "CAGACACTC",
            "AGCGTGTAT", "CCATCCTAA", "TAGTAGAGC", "GACGTAGAG", "GGTGGTATC",
            "ATGAGTCCT", "TCCTAAGCC", "GATCTCTTG", "GAGTGATCT", "CTTAGCGAC",
            "GAGAGGAGT", "ACGTACTGT", "TGTACGCCT", "CGAAGGCAT", "AGTCTACGC",
            "CTGCGTAGT", "GTCGCCTAT", "TACGCGTAC", "TGTGCAAGT", "TTACAGAGG",
            "ACTATCGCC", "CATTCCGCT", "CGTGGCTAT", "TTCTTCGTG", "CAGTAAGAG",
            "GAACCTGTA", "ACGTCCATA", "GAACGAATG", "GCCACATTG", "AGTGCTACA",
            "GAAGCCATT", "GAGTTCGTC", "AACCGTGAT", "TGACAGCTA", "CGTGTGATG",
            "ATTGAGAGC", "CAAGGTGTT", "TAGAAGGCG", "CGACTGAGA", "AGAATCCGG",
            "AGAGAATGG", "TTGTGTACG", "TGGCAATTC", "AGAATGACC", "CAAGAGTAG",
            "GATGCACAT", "ACGCCACTT", "TCGTCATGC", "GCACTATTC", "AACAGGAAC",
            "ATTAGCCTG", "ACCAATCCG", "GTAGGTAAG", "CTCTGCTTA", "CTATGACCA",
            "GTGGAATAG", "GACTAACGG", "ACGGTGAAG", "GACATGGCT", "TCCAGACGA",
            "GCATTAGCA", "AAGGAGTCT", "TTCAACAGC", "TGACCAGCA", "ATCGCTCTG",
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACGC", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATC", "GATTGTCACTAACGCGATGCTGACTCCTAGTCC"
    ),
    "flv_rna":(
          [
            "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
            "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
            "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
            "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
            "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
            "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
            "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
            "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
            "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
            "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
            "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
            "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
            "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
            "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
            "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
            "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
            "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
            "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
            "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
            "GCCAAGAC"
        ],
        "ATCCACGTGCTTGAGATCAGCATGCGGCTACGCTGTCT", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGACTGTCT",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATCTGTCT", "GATTGTCACTAACGCGATGCTGACTCCTAGTCCTGTCT"
    ),
    "flv":(
          [
              "AACGTGAT", "AAACATCG", "ATGCCTAA", "AGTGGTCA", "ACCACTGT",
              "ACATTGGC", "CAGATCTG", "CATCAAGT", "CGCTGATC", "ACAAGCTA",
              "CTGTAGCC", "AGTACAAG", "AACAACCA", "AACCGAGA", "AACGCTTA",
              "AAGACGGA", "AAGGTACA", "ACACAGAA", "ACAGCAGA", "ACCTCCAA",
              "ACGCTCGA", "ACGTATCA", "ACTATGCA", "AGAGTCAA", "AGATCGCA",
              "AGCAGGAA", "AGTCACTA", "ATCCTGTA", "ATTGAGGA", "CAACCACA",
              "GACTAGTA", "CAATGGAA", "CACTTCGA", "CAGCGTTA", "CATACCAA",
              "CCAGTTCA", "CCGAAGTA", "CCGTGAGA", "CCTCCTGA", "CGAACTTA",
              "CGACTGGA", "CGCATACA", "CTCAATGA", "CTGAGCCA", "CTGGCATA",
              "GAATCTGA", "CAAGACTA", "GAGCTGAA", "GATAGACA", "GCCACATA",
              "GCGAGTAA", "GCTAACGA", "GCTCGGTA", "GGAGAACA", "GGTGCGAA",
              "GTACGCAA", "GTCGTAGA", "GTCTGTCA", "GTGTTCTA", "TAGGATGA",
              "TATCAGCA", "TCCGTCTA", "TCTTCACA", "TGAAGAGA", "TGGAACAA",
              "TGGCTTCA", "TGGTGGTA", "TTCACGCA", "AACTCACC", "AAGAGATC",
              "AAGGACAC", "AATCCGTC", "AATGTTGC", "ACACGACC", "ACAGATTC",
              "AGATGTAC", "AGCACCTC", "AGCCATGC", "AGGCTAAC", "ATAGCGAC",
              "ATCATTCC", "ATTGGCTC", "CAAGGAGC", "CACCTTAC", "CCATCCTC",
              "CCGACAAC", "CCTAATCC", "CCTCTATC", "CGACACAC", "CGGATTGC",
              "CTAAGGTC", "GAACAGGC", "GACAGTGC", "GAGTTAGC", "GATGAATC",
              "GCCAAGAC"
          ],
          "CGTAGCCGCATGCTGATCTCAAGCACGTGGAT", "TCAGCTTCTGACTACGATATGGCTGTCACCGA",
          "ATGCTAATACGTAGTCGAGACCTACATGTTCG", "GACTAGGAGTCAGCATCGCGTTAGTGACAATC"
    ),
    "bulk_vdj":(
          [
              "CTCCAT", "ATCCTC", "ACGTCT", "TGCGAA", "TTCTCG",
              "CTGCTA", "CATGAT", "TCAACT", "AGTCCT", "GTTGAG",
              "TAGCTG", "TCGCCA", "GAACTC", "TATGGT", "CGCAAC",
              "TGGCAG", "ATGCAT", "GACTAT", "GTGATT", "CTCTTG",
              "AAGCGT", "CCAACA", "GTTGGT", "TCTAGT", "TATGTG",
              "TGTGGC", "GACCTG", "TTCCGT", "AAGGCA", "TAGGAT",
              "AACTCC", "TCGATG", "CTGCGT", "GTCGGA", "TCACAT",
              "ATAGGT", "CGTAAT", "GCATGT", "AACTGA", "GCACAA",
              "GCCATC", "CAACCG", "GTCTGG", "TCCATT", "CAGACC",
              "ACGGAG", "ACATCA", "TATCCG", "GGAGAG", "CCAATG",
              "TTCTGA", "GTGACG", "ATGGTG", "ACTTGT", "ATAGAC",
              "CCTATA", "TTAAGG", "GATCAC", "TAGCCT", "AGCGCT",
              "AGACGC", "CTAAGA", "TATCGA", "CGCACA", "CAAGTT",
              "GAACCA", "TACACA", "CATTGG", "TCATGC", "AGGTTA",
              "TCGAAT", "TCTTGG", "CTCTAC", "GAGGTC", "ACAACG",
              "CAGATA", "CAGGTA", "TCTTAC", "CCTGTG", "TCGAGC",
              "CTGAAT", "ATTGGC", "CATCTT", "TCTCTA", "GCGTCA",
              "GTTCAT", "AATCAG", "CGGTGT", "TCCGTC", "CTCACC",
              "TTGACT", "GCCGTA", "CGACTC", "ATCCAA", "TGCCAT",
              "ACGATA"
          ],
          "GTGGTATCAACGCAGAGT"
    ),
    "bulk_rna":(
          [
              "AACGGACCT", "AGGACTCAC", "ACTGCCTAG", "GAACGCTAT", "GACTGGTTG",
              "CAGGACTTC", "TCGGTTCGT", "GTCTTGCGT", "TATCTCCGA", "CGCAACTAC",
              "CACTTCAGA", "GCTCTCACT", "ACGAAGCTC", "GTGTTAAGG", "GGCTCTCTA",
              "GCGTAGTAA", "AGCTCCTTG", "CACATTCAC", "TCCGTATCA", "CACCTGTAA",
              "ATGGTCTCA", "TGCATCAAG", "CGCCAATGA", "TCGACTGTC", "TGTGGACAC",
              "ACATCGGAC", "GCTTGAGGT", "GCCGTTATT", "AACACCGTT", "CCTAGTCTT",
              "AACACACAG", "CAACATCGG", "TACCTCTCC", "CCAATGTCT", "CGAGATAGT",
              "AGTTCAGAG", "ATCGACACG", "GACCTTAGC", "TAACCTACC", "TGGAGAACC",
              "GATGTTACG", "TGGCATGAG", "CTGGTACTT", "GCGAGTAAC", "GATCCATGC",
              "GGCTTCTCA", "CGTTAGCGT", "ACACAGGCT", "AACTGGCGA", "AGACGTTCA",
              "TGCGGTTCT", "CAGTCTTCG", "AGCTGAGTC", "GTGCATATC", "GATGGCTCA",
              "TGAGCGAAG", "CTTGACGTT", "GTTACTGGT", "GAGCAGCTT", "GTTGGAGTG",
              "TGCCTGATC", "GCAGATGTG", "CAGAGTACA", "TGTAGTGTG", "ACGAATGGA",
              "TCACTGGAA", "ATGACAGCA", "CTCAGAACT", "AAGCTTGCG", "TCGGACATG",
              "GAGGATTGA", "TCTTGGACA", "GTGCAAGGT", "CGATCGGTA", "GGATACCAC",
              "ATCGTTGGC", "TGGAACGTA", "GCCTACGAT", "CATGTAGGC", "CTCACGTTC",
              "TGTGAGTCA", "AATCGCCAC", "GTGCGACTA", "TCTGGCGAA", "ATACGCGGA",
              "CCTTGAATC", "CCACACATT", "CGGTGATTG", "AGGAGCAAT", "ACACACCAA",
              "GGATAGATG", "TACCGTCTG", "TCTTGCTTG", "CCAGCTAAC", "CTTACGCAG",
              "AGTAGGAGG"
          ],
          []
    ),
    "scope_5p3p":(
          [
              "AACGGACCT", "AGGACTCAC", "ACTGCCTAG", "GAACGCTAT", "GACTGGTTG",
              "CAGGACTTC", "TCGGTTCGT", "GTCTTGCGT", "TATCTCCGA", "CGCAACTAC",
              "CACTTCAGA", "GCTCTCACT", "ACGAAGCTC", "GTGTTAAGG", "GGCTCTCTA",
              "GCGTAGTAA", "AGCTCCTTG", "CACATTCAC", "TCCGTATCA", "CACCTGTAA",
              "ATGGTCTCA", "TGCATCAAG", "CGCCAATGA", "TCGACTGTC", "TGTGGACAC",
              "ACATCGGAC", "GCTTGAGGT", "GCCGTTATT", "AACACCGTT", "CCTAGTCTT",
              "AACACACAG", "CAACATCGG", "TACCTCTCC", "CCAATGTCT", "CGAGATAGT",
              "AGTTCAGAG", "ATCGACACG", "GACCTTAGC", "TAACCTACC", "TGGAGAACC",
              "GATGTTACG", "TGGCATGAG", "CTGGTACTT", "GCGAGTAAC", "GATCCATGC",
              "GGCTTCTCA", "CGTTAGCGT", "ACACAGGCT", "AACTGGCGA", "AGACGTTCA",
              "TGCGGTTCT", "CAGTCTTCG", "AGCTGAGTC", "GTGCATATC", "GATGGCTCA",
              "TGAGCGAAG", "CTTGACGTT", "GTTACTGGT", "GAGCAGCTT", "GTTGGAGTG",
              "TGCCTGATC", "GCAGATGTG", "CAGAGTACA", "TGTAGTGTG", "ACGAATGGA",
              "TCACTGGAA", "ATGACAGCA", "CTCAGAACT", "AAGCTTGCG", "TCGGACATG",
              "GAGGATTGA", "TCTTGGACA", "GTGCAAGGT", "CGATCGGTA", "GGATACCAC",
              "ATCGTTGGC", "TGGAACGTA", "GCCTACGAT", "CATGTAGGC", "CTCACGTTC",
              "TGTGAGTCA", "AATCGCCAC", "GTGCGACTA", "TCTGGCGAA", "ATACGCGGA",
              "CCTTGAATC", "CCACACATT", "CGGTGATTG", "AGGAGCAAT", "ACACACCAA",
              "GGATAGATG", "TACCGTCTG", "TCTTGCTTG", "CCAGCTAAC", "CTTACGCAG",
              "AGTAGGAGG",
          
              "TTGGTGACC", "CATCGGTTC", "CAATGCAAC", "AAGGTGGTA", "TCAGGTAGA",
              "CACTAGGCA", "GGCATGCAA", "ATGCGATAC", "AGACGAAGT", "CGATAAGGC",
              "CGAGTTGCA", "GGTGATCAG", "AGTGGTGGT", "TAGCGATGA", "ATGCCTAAG",
              "CCATAATCG", "CTAACCAGA", "TGATGTGCC", "TAGCTACAC", "TCCTGGCTT",
              "GTAGTTCCT", "GTATCCTTC", "GCCATAACC", "GCATGTTGG", "AATCCGGTG",
              "ACGTTACGA", "AGCATAGCG", "GGTCCGTAA", "GTTCTACCG", "CGCTGTAAG",
              "CGTCATACC", "TGGTAACCG", "ACAACAGGT", "TGACTTCCG", "CACGCAATA",
              "ACCGTACTC", "GATGTGTGT", "ACACCAACG", "TTGAGACAG", "CGGATCATC",
              "TAATGGCCG", "ACCTCGACT", "CGGCTAGAT", "TGGACTTGT", "TGATCCTCT",
              "ATAGCGTGT", "TGTCGGTGT", "AGCCACATA", "GAAGAAGCC", "CTGTGGTAT",
              "CAGCCGAAT", "ACCTGCTAC", "TCTTCTCAG", "GATCAGGAC", "TAGACCACT",
              "TGAGTAGTC", "TTCGAGGAT", "ATGTATCGG", "TCTGTCTGC", "GCTTACAGG",
              "ATACCAGTC", "AACGTCCAA", "CTCCTCAAT", "TTCCAATCG", "AACGCTAGT",
              "TTACACGAC", "CATACGACC", "CTAATCGCG", "CGTAATTGG", "GTAGTGTTG",
              "TATAGCGGT", "CGTACTGAA", "TCGATGTGG", "CCAGAAGAT", "AGGCTGTTG",
              "TGAGGCCTT", "CCAGTCCTA", "CAATTGCGC", "TTGCCGTCA", "GAGTTGACA",
              "TGGTCAGTT", "GGTAGTCCA", "ATCCTTCCA", "ATGTGCAGC", "GGAAGACTC",
              "GTAATGGAC", "TATCGTGCA", "CCTACAAGG", "TAGTCCGGA", "AGATACGCA",
              "ACTCATCGT", "GATAACCGC", "GCTGCGATA", "TACTCACCA", "ACCAGGTCA",
              "GCAACTTCA",
          
              "TTCGGTCAA", "GTACGGACT", "GGCCATGTT", "TATGACACC", "GTGTCTGAA",
              "TAAGCTTGG", "AAGATCTGC", "CTTGTGCCA", "ACTGGTTCC", "GATTGCGAG",
              "ACACGGTAG", "TCGCATACT", "CCGATCGAA", "GTAGCACGA", "TCAGCACTG",
              "TGCTATCGC", "GTAAGATCG", "TGCACGAGA", "GAATCGCAA", "CCTCGATCA",
              "GCGTCTAGT", "TCCAGTTAG", "CAGGAACCA", "GTGAACCAA", "CTCTAACAC",
              "GTCTGATGC", "TTAAGAGGC", "GTCAATGCA", "CTTCATGGA", "CAGACACTC",
              "AGCGTGTAT", "CCATCCTAA", "TAGTAGAGC", "GACGTAGAG", "GGTGGTATC",
              "ATGAGTCCT", "TCCTAAGCC", "GATCTCTTG", "GAGTGATCT", "CTTAGCGAC",
              "GAGAGGAGT", "ACGTACTGT", "TGTACGCCT", "CGAAGGCAT", "AGTCTACGC",
              "CTGCGTAGT", "GTCGCCTAT", "TACGCGTAC", "TGTGCAAGT", "TTACAGAGG",
              "ACTATCGCC", "CATTCCGCT", "CGTGGCTAT", "TTCTTCGTG", "CAGTAAGAG",
              "GAACCTGTA", "ACGTCCATA", "GAACGAATG", "GCCACATTG", "AGTGCTACA",
              "GAAGCCATT", "GAGTTCGTC", "AACCGTGAT", "TGACAGCTA", "CGTGTGATG",
              "ATTGAGAGC", "CAAGGTGTT", "TAGAAGGCG", "CGACTGAGA", "AGAATCCGG",
              "AGAGAATGG", "TTGTGTACG", "TGGCAATTC", "AGAATGACC", "CAAGAGTAG",
              "GATGCACAT", "ACGCCACTT", "TCGTCATGC", "GCACTATTC", "AACAGGAAC",
              "ATTAGCCTG", "ACCAATCCG", "GTAGGTAAG", "CTCTGCTTA", "CTATGACCA",
              "GTGGAATAG", "GACTAACGG", "ACGGTGAAG", "GACATGGCT", "TCCAGACGA",
              "GCATTAGCA", "AAGGAGTCT", "TTCAACAGC", "TGACCAGCA", "ATCGCTCTG",
              "CCTTAGGTG"
          ],
          []
    )
    
}

pattern_regex = re.compile(r'([CLUT])(\d+)')

def parse_pattern(pattern):
    
    ranges = {'C': [], 'L': [], 'U': [], 'T': []}
    matches = pattern_regex.findall(pattern)
    for match in matches:
        char, length = match
        ranges[char].append(int(length))
    return ranges

def determine_pattern(r1_sequence):
    """
    根据 R1 序列判断其所属的模式。
    """
    for pattern_name, (bclist, linker) in bclist_linker_map.items():
        pattern = PATTERN_DICT.get(pattern_name)
        if not pattern:
            continue

        ranges = parse_pattern(pattern)
        pos = 0
        is_match = True

        for i in range(len(ranges['C'])):
            c_length = ranges['C'][i]
            c_segment = r1_sequence[pos:pos + c_length]
            if c_segment not in bclist:
                is_match = False
                break
            pos += c_length

            if i < len(ranges['L']) and linker:
                l_length = ranges['L'][i]
                l_segment = r1_sequence[pos:pos + l_length]
                if l_segment != linker:
                    is_match = False
                    break
                pos += l_length

        if is_match:
            return pattern_name

    # 若不匹配任何已知模式，返回 "auto"
    return "auto"

def split_and_modify_fastq(r1_filename, r2_filename, output_r1_filename, output_r2_filename):
    """
    根据模式切割和重组 R1 和 R2 文件。若为 auto 模式，则不处理原始数据。
    """
    with pysam.FastxFile(r1_filename) as r1_file, \
         pysam.FastxFile(r2_filename) as r2_file, \
         gzip.open(output_r1_filename, "wt") as output_r1_file, \
         gzip.open(output_r2_filename, "wt") as output_r2_file:

        for r1_read, r2_read in zip(r1_file, r2_file):
            r1_sequence = r1_read.sequence
            r2_sequence = r2_read.sequence

            # 判断 R1 的模式
            pattern_name = determine_pattern(r1_sequence)

            if pattern_name == "auto":
                # 如果是 auto 模式，不做任何处理，直接将原始的 R1 和 R2 写入输出文件
                SeqIO.write(SeqRecord(Seq(r1_sequence), id=r1_read.name, description=""), output_r1_file, "fastq")
                SeqIO.write(SeqRecord(Seq(r2_sequence), id=r2_read.name, description=""), output_r2_file, "fastq")
                continue

            # 若非 auto 模式，进行切割和重组
            pattern = PATTERN_DICT.get(pattern_name)
            if not pattern:
                click.echo(f"Pattern '{pattern_name}' has no defined pattern model. Skipping read {r1_read.name}.", err=True)
                continue

            ranges = parse_pattern(pattern)
            pos = 0
            positions = {'C': [], 'L': [], 'U': None, 'T': None}
            for length in ranges['C']:
                positions['C'].append((pos, pos + length))
                pos += length
            if ranges['L']:
                for length in ranges['L']:
                    positions['L'].append((pos, pos + length))
                    pos += length
            if ranges['U']:
                positions['U'] = (pos, pos + ranges['U'][0])
                pos += ranges['U'][0]
            if ranges['T']:
                positions['T'] = (pos, pos + ranges['T'][0])
                pos += ranges['T'][0]

            # 切割 R1
            barcode_segments = [r1_sequence[start:end] for start, end in positions['C']]
            linker_segments = [r1_sequence[start:end] for start, end in positions['L']]
            umi = r1_sequence[positions['U'][0]:positions['U'][1]] if positions['U'] else ""
            polyT = r1_sequence[positions['T'][0]:positions['T'][1]] if positions['T'] else ""
            remaining_sequence = r1_sequence[positions['T'][1]:] if positions['T'] else r1_sequence[positions['U'][1]:]

            # 重构 R1 序列
            new_r1_seq = ''.join(barcode_segments) + ''.join(linker_segments) + umi + polyT

            # 重构 R2 序列
            new_r2_seq = umi + remaining_sequence + r2_sequence

            # 创建新的 SeqRecord 对象
            new_r1_record = SeqRecord(Seq(new_r1_seq), id=r1_read.name, description="")
            new_r2_record = SeqRecord(Seq(new_r2_seq), id=r2_read.name, description="")

            # 写入到输出文件
            SeqIO.write(new_r1_record, output_r1_file, "fastq")
            SeqIO.write(new_r2_record, output_r2_file, "fastq")

@click.command()
@click.option('--r1-filename', required=True, type=click.Path(exists=True), help="Path to R1.fastq.gz")
@click.option('--r2-filename', required=True, type=click.Path(exists=True), help="Path to R2.fastq.gz")
@click.option('--output-r1-filename', required=True, type=click.Path(), help="Path to output R1 fastq file")
@click.option('--output-r2-filename', required=True, type=click.Path(), help="Path to output R2 fastq file")
def cli(r1_filename, r2_filename, output_r1_filename, output_r2_filename):
    
    split_and_modify_fastq(r1_filename, r2_filename, output_r1_filename, output_r2_filename)

if __name__ == "__main__":
    cli()