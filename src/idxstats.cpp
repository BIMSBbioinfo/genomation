/*  This is a portion of code taken from:
        
        bam_index.c -- index and idxstats subcommands.
   
   Copyright (C) 2008-2011, 2013, 2014 Genome Research Ltd.
   Portions copyright (C) 2010 Broad Institute.
   Portions copyright (C) 2013 Peter Cock, The James Hutton Institute.
  
   Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */


#include "hts.h"
#include "sam.h"
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>
#include <Rcpp.h>

using namespace Rcpp;

Rcpp::DataFrame bam_idxstats(std::string &bam_file)
{
    hts_idx_t *idx = NULL;
    bam_hdr_t *header = NULL;
    htsFile *fp = NULL ;

    if (bam_file == "" ) {
        fprintf(stderr, "Usage: samtools idxstats <in.bam>\n");
        return 1;
    }
    fp = sam_open(bam_file.c_str(),"r");
    if (fp == NULL) { fprintf(stderr, "[%s] fail to open BAM.\n", __func__); return 1; }
    header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "[%s] failed to read header for '%s'.\n",
                __func__, bam_file.c_str());
        return 1;
    }
    idx = sam_index_load(fp, bam_file.c_str());
    if (idx == NULL) { fprintf(stderr, "[%s] fail to load the index.\n", __func__); return 1; }

    int i;
    Rcpp::CharacterVector name(header->n_targets);
    Rcpp::IntegerVector length(header->n_targets),
                        mapped(header->n_targets),
                        unmapped(header->n_targets);
    for (i = 0; i < header->n_targets; ++i) {
        // Print out contig name and length
        // printf("%s\t%d", header->target_name[i], header->target_len[i]);
        // Store contig name and length in vector
        name[i]= header->target_name[i]; 
        length[i]= (int) header->target_len[i];
        // Now fetch info about it from the meta bin
        uint64_t u, v;
        hts_idx_get_stat(idx, i, &u, &v); // u/v are mapped/unmapped reads respectively
        //printf("\t%" PRIu64 "\t%" PRIu64 "\n", u, v);
        mapped[i]= (int) u;
        unmapped[i]= (int) v; 
    }
    // Dump information about unmapped reads
    // printf("*\t0\t0\t%" PRIu64 "\n", hts_idx_get_n_no_coor(idx));
    
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    sam_close(fp);
    // see http://gallery.rcpp.org/articles/modifying-a-data-frame/
    return DataFrame::create(_["name"]=name,
                             _["length"]=length,
                             _["mapped"]=mapped,
                             _["unmapped"]=unmapped);
}


//' Fast retrieval of indexstats
//' 
//' Retrieve and print stats in the index file corresponding 
//' to the input file. Before calling idxstats, the input BAM file
//' must be indexed, e.g. via samtools index.
//' 
//' @param bam_file path to the indexed BAM file 
//' @export 
//[[Rcpp::export]]
Rcpp::DataFrame idxStats(std::string bam_file)
{
  return(bam_idxstats(bam_file));
}