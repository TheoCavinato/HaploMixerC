/*******************************************************************************
 * Copyright (C) 2020 Théo Cavinato, University of Lausanne
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <io/genotype_reader_writer.h>

genotype_reader_writer::genotype_reader_writer(){}

genotype_reader_writer::~genotype_reader_writer(){}

void genotype_reader_writer::readAndWriteGenotypes(string fvcfin, string fvcfout, string region, vector<vector<int>> &parents_haplotypes, int lower_limit, int upper_limit){
	//-----------INITIALISE VCF TO READ---------------//
	//Initialize HTSlib VCF/BCF reader
	htsFile * fin= hts_open(fvcfin.c_str(),"r");
	bcf_hdr_t * hdr_in = bcf_hdr_read(fin);
	//To read VCF records
	unsigned int nset = 0;
	int ngt, *gt_arr = NULL, ngt_arr = 0;
	bcf1_t * rec_in = bcf_init();

	int n_samples = bcf_hdr_nsamples(hdr_in);

	//-----------INITIALISE VCF TO WRITE---------------//
	string fname = fvcfout;
	string file_format="w";
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; }
	bcf_hdr_t * hdr_out = bcf_hdr_init("w");
	htsFile * fout = hts_open(fname.c_str(),file_format.c_str());

	//Create VCF header
	bcf_hdr_append(hdr_out, string("##fileDate="+tac.date()).c_str());
	bcf_hdr_append(hdr_out, "##source=shapeit4.1.3");
	bcf_hdr_append(hdr_out, string("##contig=<ID="+region+">").c_str());
	bcf_hdr_append(hdr_out, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr_out, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr_out, "##INFO=<ID=CM,Number=A,Type=Float,Description=\"Interpolated cM position\">");
	bcf_hdr_append(hdr_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");
	for(int i = 0; i < n_samples; i++){
		bcf_hdr_add_sample(hdr_out, hdr_in->samples[i]);
	}
	bcf_hdr_add_sample(hdr_out, NULL);      // to update internal structures
	if (bcf_hdr_write(fout, hdr_out) < 0) vrb.error("Failing to write VCF/header");

	//To write VCF record
	bcf1_t *rec_out = bcf_init();
	int * genotypes_out = (int*)malloc(n_samples*2*sizeof(int));

	//-----------READ fvcfin AND WRITE fvcfout---------------//
	int n_variants = 1;
	while((nset=bcf_read(fin, hdr_in, rec_in))==0) {
		bcf_clear(rec_out);
		bcf_unpack(rec_in, BCF_UN_STR);
        if (rec_in->pos > lower_limit && rec_in->pos < upper_limit ){ //to only read positions that are between in the gmap limits
			//std::cout << rec_in->pos << std::endl;
			bool a0, a1;
			ngt = bcf_get_genotypes(hdr_in, rec_in, &gt_arr, &ngt_arr);
			int count_alt = 0;

			for(int i = 0; i < n_samples*2; i+=2){
				int pos_h0 = parents_haplotypes[n_variants][i];
				int pos_h1 = parents_haplotypes[n_variants][i+1];
				a0 = (bcf_gt_allele(gt_arr[pos_h0])==1);
				a1 = (bcf_gt_allele(gt_arr[pos_h1])==1);
				genotypes_out[i] = bcf_gt_phased(a0);
				genotypes_out[i+1] = bcf_gt_phased(a1);
				count_alt += a0+a1;
			}

			//Add all the 9 first column to genotypes_out
			rec_out->pos = rec_in->pos;
			rec_in->rid = rec_in->rid;

			string ref = string(rec_in->d.allele[0]);
			string alt = string(rec_in->d.allele[1]);
			string alleles = ref + "," + alt;
			bcf_update_alleles_str(hdr_out, rec_out, alleles.c_str());

			bcf_update_info_int32(hdr_out, rec_out, "AC", &count_alt, 1);
			float freq_alt = count_alt * 1.0 / (2 * n_samples);
			bcf_update_info_float(hdr_out, rec_out, "AF", &freq_alt, 1);

			bcf_update_id(hdr_out, rec_out, rec_in->d.id);

			//Write genotypes_out to the fout
			bcf_update_genotypes(hdr_out, rec_out, genotypes_out, n_samples*2);
			if (bcf_write(fout, hdr_out, rec_out) < 0) vrb.error("Failing to write VCF/record");

			if(n_variants % 100000 == 0){
				vrb.bullet("[ " + stb.str(n_variants) + " SNPs written (" +  stb.str(tac.rel_time()*1.0/1000, 2) + "s) ]");
				tac.clock();
			}

			n_variants++;
		}
	}

	//Free allocated memory
	free(genotypes_out);
	bcf_hdr_destroy(hdr_out);
	bcf_destroy(rec_out);

	free(gt_arr);

	if (hts_close(fin)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	if (hts_close(fout)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

};