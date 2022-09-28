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

#include <containers/recombination.h>
#include <random>
#include <algorithm>

recombination::recombination() {
}

void recombination::readBcfPos(string fvcf, int lower_limit, int upper_limit){
    tac.clock();
    htsFile * f = hts_open(fvcf.c_str(), "r");
    bcf1_t  * line_vcf = bcf_init();
    bcf_hdr_t * hdr = bcf_hdr_read(f);
	unsigned int nset = 0;

	n_samples = bcf_hdr_nsamples(hdr);

	while((nset=bcf_read(f, hdr, line_vcf))==0) {
		bcf_unpack(line_vcf, BCF_UN_STR);
        if (line_vcf->pos > lower_limit && line_vcf->pos < upper_limit ){ //to only read positions that are comprised in the gmap
            bcfPosBp.push_back(line_vcf->pos);
        }
    }

    bcf_hdr_destroy(hdr);
    bcf_destroy(line_vcf);

	if (hts_close(f)) vrb.error("Non zero status when closing VCF/BCF file descriptor");

	vrb.bullet("Retrieve bp from bcf (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void recombination::bpToRecRate(vector<int> &gmap_bp_pos, vector <double> &gmap_cm_pos){
    tac.clock();
    //step 0: find the position of each bcf_pos in gmap_bp_pos
    vector <int> bcf_to_gmap_idx_down;
    vector <int> bcf_to_gmap_idx_up;
    int bcf_itr = 0;
    int gmap_bp_itr = 0;
    while (bcf_itr < bcfPosBp.size()){
        if (bcfPosBp[bcf_itr] <= gmap_bp_pos[gmap_bp_itr+1] && bcfPosBp[bcf_itr] >= gmap_bp_pos[gmap_bp_itr]){
            bcf_to_gmap_idx_down.push_back(gmap_bp_itr);
            bcf_to_gmap_idx_up.push_back(gmap_bp_itr+1);
            bcf_itr++;
        }
        else{
            gmap_bp_itr+=1;
        }
    }

    //step 1: convet bp to cM by linear conversion using the positions found previously
    vector <double> bcf_pos_cm;
    double bcf_pos, g_bp_down, g_bp_up, g_cm_down, g_cm_up;
    for (int i = 0; i< bcf_to_gmap_idx_down.size(); i++){
        bcf_pos = bcfPosBp[i];
        g_bp_down = gmap_bp_pos[bcf_to_gmap_idx_down[i]];
        g_bp_up = gmap_bp_pos[bcf_to_gmap_idx_up[i]];
        g_cm_down = gmap_cm_pos[bcf_to_gmap_idx_down[i]];
        g_cm_up = gmap_cm_pos[bcf_to_gmap_idx_up[i]];
        bcf_pos_cm.push_back(((bcf_pos-g_bp_down)/(g_bp_up-g_bp_down))*(g_cm_up-g_cm_down)+g_cm_down);
    }

    //step 2: calculate recombination rate based on the obtained values
    bcfRecRateBtwPos.push_back(bcf_pos_cm[0]/100);
    for (int i = 1; i< bcf_pos_cm.size(); i++){
        bcfRecRateBtwPos.push_back((bcf_pos_cm[i]-bcf_pos_cm[i-1])/100);
    }

	vrb.bullet("Calculate rec rate (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void recombination::simulateRecombination(string recvalid){

    tac.clock();
    vector <int> shuffled_haplotypes (n_samples*2);
    for(int i = 0; i<n_samples*2; i++) shuffled_haplotypes[i] = i;
    std::random_shuffle(shuffled_haplotypes.begin(), shuffled_haplotypes.end());

    haploToSelect.push_back(shuffled_haplotypes);

    double random_double;

    vector<int> save_rec_sites;

    for (int i =0; i < bcfRecRateBtwPos.size(); i++){
        for (int j=0; j<n_samples*2; j+=2){
            random_double = rng.getDouble();
            // check if we change haplotype
            if (random_double < bcfRecRateBtwPos[i]){
                int h0 = shuffled_haplotypes[j], h1 = shuffled_haplotypes[j+1];
                shuffled_haplotypes[j] = h1;
                shuffled_haplotypes[j+1] = h0;
                // write recombination sites if --recvalid option activated
                if (recvalid!="None") save_rec_sites.push_back(bcfPosBp[i]);
            }
        }
        haploToSelect.push_back(shuffled_haplotypes);
    }

    // if --recvalid option is on
    if (recvalid!="None") {
        output_file fd(recvalid);
        for(int i : save_rec_sites) fd << i << endl;
        fd.close();
    }

	vrb.bullet("Simulate recombination sites (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

}

recombination::~recombination() {
}