#include <iostream>

#include "gtest/gtest.h"
#include "matepair.h"

TEST(HammingMatch,PerfectMatches)
{

    string adapter1 = "CTGTCTCTTATACACATCT";
    string adapter2 = "AGATGTGTATAAGAGACAG";
    
    int min_overlap=12;
    float similarity=0.85;
    string target = "AGATGTGTATAAGAGACAGGCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGA";
    string query = "AGATGTGTATAAGAGACAG";
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),0);
    
    target = "GCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGACGCCAGATGTGTATAAGAGACAG";
    query = "AGATGTGTATAAGAGACAG";
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),target.size()-query.size());


    string s1 = "GCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCC";    
    string s2 = "TGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGG";
    query = adapter1+adapter2;
    target = s1 + query + s2;
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),s1.size());    
}


TEST(HammingMatch,NoisyMatches)
{

    string adapter1 = "CTGTCTCTTATACACATCT";
    string adapter2 = "AGATGTGTATAAGAGACAG";
    
    int min_overlap=12;
    float similarity=0.85;
    string target = "AGGTGTGCATAAGAGACAGGCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGA";
    string query = "AGATGTGTATAAGAGACAG";
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),0);
    
    target = "GCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGACGCCAGATGTGTATATCAGACAG";
    query = "AGATGTGTATAAGAGACAG";
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),target.size()-query.size());


    string s1 = "GCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCC";    
    string s2 = "TGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGG";
    query = "CTGTCTCTTATAGTCATCTAGATGTGTATAAGACACAG";
    target = s1 + query + s2;
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),s1.size());    
}

TEST(HammingMatch,PartialMatches)
{

    string adapter1 = "CTGTCTCTTATACACATCT";
    string adapter2 = "AGATGTGTATAAGAGACAG";
   
    int min_overlap=12;
    float similarity=0.85;
    string target =adapter2.substr(4) +  "GATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGA";
    string query = adapter2;
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),-4);
    
    target = "GCGATTTCCTTACATTGACGTTTTTATTACTCACTGTCCTGTTCCTGTTATCACATTATCTGCTGAACAATTACTGATAGGTTAAAGAGAACCAGGCCTGGGCATTGGCGATGCCGCCAGTAACCCGACG" + adapter1.substr(0,15);
    query = adapter1;
    ASSERT_EQ(hamming_match(target,query,min_overlap,similarity),target.size()-15);
}

