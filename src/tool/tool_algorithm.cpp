#include <tool/tool_header.h>

#define MODE_MAF	0
#define MODE_MAC	1
#define MODE_AF		2
#define MODE_AC		3

void tool::runMainTask() {
	vrb.title("Compute main TASK");

	//simulate recombination sites
	REC_SITES.readBcfPos(options["vcf"].as < string > (), GMAP.pos_bp[0], GMAP.pos_bp.back());
	REC_SITES.bpToRecRate(GMAP.pos_bp, GMAP.pos_cm);

	string valid_path="None";
	if (options.count("recvalid")) valid_path = options["recvalid"].as < string > ();
	REC_SITES.simulateRecombination(valid_path);

	//read bcf and write simulated individuals
	GEN.readAndWriteGenotypes(options["vcf"].as < string > (), options["output"].as <string> (), options["region"].as < string > (), REC_SITES.haploToSelect, GMAP.pos_bp[0], GMAP.pos_bp.back());

}
