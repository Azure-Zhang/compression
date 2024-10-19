
# -*- coding: utf-8 -*-
import click
import pysam
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

PATTERN_DICT = {
    "auto": None,
    "scopeV1": {"pattern": "C12U8T18", "total_length": 12 + 8 + 18},  # 38 bases
    "scopeV2.0.0": {"pattern": "C8L16C8L16C8U8T18", "total_length": 8 + 16 + 8 + 16 + 8 + 8 + 18},  # 82 bases
    "scopeV2.0.1": {"pattern": "C8L16C8L16C8L1U8T18", "total_length": 8 + 16 + 8 + 16 + 8 + 1 + 8 + 18},
    "scopeV2.1.0": {"pattern": "C8L16C8L16C8U12T18", "total_length": 8 + 16 + 8 + 16 + 8 + 12 + 18},
    "scopeV2.1.1": {"pattern": "C8L16C8L16C8L1U12T18", "total_length": 8 + 16 + 8 + 16 + 8 + 1 + 12 + 18},
    "scopeV2.2.1": {"pattern": "C8L16C8L16C8L1U12T18", "total_length": 8 + 16 + 8 + 16 + 8 + 1 + 12 + 18},
    "scopeV3.0.1": {"pattern": "C9L16C9L16C9L1U12T18", "total_length": 9 + 16 + 9 + 16 + 9 + 1 + 12 + 18},
    "flv_rna": {"pattern": "C8L16C8L16C8U9L6", "total_length": 8 + 16 + 8 + 16 + 8 + 9 + 6},
    "flv": {"pattern": "U9C8L16C8L16C8", "total_length": 9 + 8 + 16 + 8 + 16 + 8},
    "bulk_vdj": {"pattern": "L18C6U16", "total_length": 18 + 6 + 16},
    "bulk_rna": {"pattern": "C9U12", "total_length": 9 + 12},
    "scope_5p3p": {"pattern": "C9C9C9U8", "total_length": 9 + 9 + 9 + 8}

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
        ["ATCCACGTGCTTGAGATCAGCATGCGGCTACGC", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATC", "GATTGTCACTAACGCGATGCTGACTCCTAGTCC"]
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
        ["ATCCACGTGCTTGAGATCAGCATGCGGCTACGC", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGAC",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATC", "GATTGTCACTAACGCGATGCTGACTCCTAGTCC"]
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
        ["ATCCACGTGCTTGAGATCAGCATGCGGCTACGCTGTCT", "TCGGTGACAGCCATATCGTAGTCAGAAGCTGACTGTCT",
        "CGAACATGTAGGTCTCGACTACGTATTAGCATCTGTCT", "GATTGTCACTAACGCGATGCTGACTCCTAGTCCTGTCT"]
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
          ["CGTAGCCGCATGCTGATCTCAAGCACGTGGAT", "TCAGCTTCTGACTACGATATGGCTGTCACCGA",
          "ATGCTAATACGTAGTCGAGACCTACATGTTCG", "GACTAGGAGTCAGCATCGCGTTAGTGACAATC"]
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
    
    #����ģʽ�ַ�������ȡC, L, U, T�ĳ��ȣ����������еķָ�㡣
    #����ÿ��Ƭ�εĳ��ȼ��ܳ��ȡ�
    
    if not isinstance(pattern, str):
        raise TypeError("Expected string or bytes-like object, got {}".format(type(pattern)))

    ranges = {'C': [], 'L': [], 'U': [], 'T': []}
    total_length = 0
    matches = pattern_regex.findall(pattern)
    for match in matches:
        char, length = match
        length = int(length)
        ranges[char].append(length)
        total_length += length  # �����ܳ���
    return ranges, total_length

def match_bclist_and_linker(r1_sequence, pattern_name, cut_positions):
    
    #��CƬ����bclistƥ�䣬��LƬ����linkerƥ�䣬ƥ��ɹ��򷵻�True�����򷵻�False��
    
    bclist, linker = bclist_linker_map.get(pattern_name, (None, None))
    
    if bclist is None:
        return False

    # ��r1_sequence����ȡC��LƬ�ν���ƥ��
    c_segments = [r1_sequence[start:end] for start, end in cut_positions['C']]
    l_segments = [r1_sequence[start:end] for start, end in cut_positions['L']] if cut_positions['L'] else []

    # ����ƥ��ÿ��CƬ���Ƿ���bclist��
    for c_segment in c_segments:
        if c_segment not in bclist:
            return False

    # ���linkerΪ�գ�������LƬ��ƥ��
    if linker is None:
        return True

    # ƥ��LƬ���Ƿ����linker
    for l_segment in l_segments:
        if l_segment != linker:
            return False

    return True

def determine_pattern_and_cut_point(r1_sequence):
    #���� R1 �����ж���������ģʽ�������ؼ���ķָ�㡣
    
    for pattern_name, pattern in PATTERN_DICT.items():
        if pattern is None or not isinstance(pattern, str):  # ���ģʽ�Ƿ�Ϊ None ����ַ��������Ϊ None ����
            continue

        ranges, cut_point = parse_pattern(pattern)  # ����ģʽ
        # ����C��L���и��
        cut_positions = {'C': [], 'L': [], 'U': None, 'T': None}
        pos = 0
        for length in ranges['C']:
            cut_positions['C'].append((pos, pos + length))
            pos += length
        for length in ranges['L']:
            cut_positions['L'].append((pos, pos + length))
            pos += length

        # ƥ��bclist��linker
        if match_bclist_and_linker(r1_sequence, pattern_name, cut_positions):
            return pattern_name, cut_point, cut_positions

    return "auto", None, None

def split_and_modify_fastq(r1_filename, r2_filename, output_r1_filename, output_r2_filename, reverse=False):
    
    #����ģʽ�и������ R1 �� R2 �ļ�����Ϊ auto ģʽ���򲻴���ԭʼ���ݡ�
    #���reverse=True�����������ļ��ָ�Ϊԭ�ļ���
    
    default_quality = 30  # ����Ĭ�ϵ�����ֵ

    with pysam.FastxFile(r1_filename) as r1_file, \
         pysam.FastxFile(r2_filename) as r2_file, \
         gzip.open(output_r1_filename, "wt") as output_r1_file, \
         gzip.open(output_r2_filename, "wt") as output_r2_file:

        for r1_read, r2_read in zip(r1_file, r2_file):
            r1_sequence = r1_read.sequence
            r2_sequence = r2_read.sequence
            
            # ������������Ƿ�Ϊ�գ����Ϊ��������Ĭ����������
            r1_quality = r1_read.quality if r1_read.quality is not None else [default_quality] * len(r1_sequence)
            r2_quality = r2_read.quality if r2_read.quality is not None else [default_quality] * len(r2_sequence)

            # �ж� R1 ��ģʽ�ͷָ��
            pattern_name, cut_point, cut_positions = determine_pattern_and_cut_point(r1_sequence)

            if pattern_name == "auto" or cut_point is None:
                # ����� auto ģʽ��δ�ҵ�ģʽ�������κδ���ֱ�ӽ�ԭʼ�� R1 �� R2 д������ļ�
                SeqIO.write(SeqRecord(Seq(r1_sequence), id=r1_read.name, description="", letter_annotations={"phred_quality": r1_quality}), output_r1_file, "fastq")
                SeqIO.write(SeqRecord(Seq(r2_sequence), id=r2_read.name, description="", letter_annotations={"phred_quality": r2_quality}), output_r2_file, "fastq")
                continue

            # ����ǻָ�ģʽ��reverse=True���ָ��ļ�����Ϊԭʼ״̬
            if reverse:
                new_r1_seq = r1_sequence + r2_sequence[:cut_point]  # �ָ�R1����
                new_r1_quality = r1_quality + r2_quality[:cut_point]  # �ָ�R1��������
                new_r2_seq = r2_sequence[cut_point:]  # ʣ�ಿ��ΪR2����
                new_r2_quality = r2_quality[cut_point:]  # ʣ�ಿ��ΪR2��������
            else:
                # ʹ�÷ָ������и�
                new_r1_seq = r1_sequence[:cut_point]
                new_r1_quality = r1_quality[:cut_point]
                remaining_sequence = r1_sequence[cut_point:]
                remaining_quality = r1_quality[cut_point:]
                new_r2_seq = remaining_sequence + r2_sequence
                new_r2_quality = remaining_quality + r2_quality

            # �����µ� SeqRecord ���󣬲�ȷ������������Ч
            new_r1_record = SeqRecord(Seq(new_r1_seq), id=r1_read.name, description="", letter_annotations={"phred_quality": new_r1_quality})
            new_r2_record = SeqRecord(Seq(new_r2_seq), id=r2_read.name, description="", letter_annotations={"phred_quality": new_r2_quality})

            # д�뵽����ļ�
            SeqIO.write(new_r1_record, output_r1_file, "fastq")
            SeqIO.write(new_r2_record, output_r2_file, "fastq")
            
@click.command()
@click.option('--r1-filename', required=True, type=click.Path(exists=True), help="Path to R1.fastq.gz")
@click.option('--r2-filename', required=True, type=click.Path(exists=True), help="Path to R2.fastq.gz")
@click.option('--output-r1-filename', required=True, type=click.Path(), help="Path to output R1 fastq file")
@click.option('--output-r2-filename', required=True, type=click.Path(), help="Path to output R2 fastq file")
@click.option('--reverse', is_flag=True, help="Set to reverse the process and restore the original file")
def cli(r1_filename, r2_filename, output_r1_filename, output_r2_filename, reverse):
    
    #ʹ�ø�����R1��R2�ļ��������ݲ�ͬģʽ���зָ������顣��Ϊautoģʽ���򲻴���ֱ�������
    
    split_and_modify_fastq(r1_filename, r2_filename, output_r1_filename, output_r2_filename, reverse)

if __name__ == "__main__":
    cli()