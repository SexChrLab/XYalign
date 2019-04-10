from __future__ import print_function
from xyalign import variants as xyv
import numpy as np
import sys

out_file = sys.argv[1]
path_prefix = sys.argv[2]

if path_prefix[-1] == "/":
	path_prefix = path_prefix[:-1]

amp = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.ampliconic.vcf.gz".format(path_prefix))
het = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.heterochromatic.vcf.gz".format(path_prefix))
oth = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.other.vcf.gz".format(path_prefix))
par = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.par.vcf.gz".format(path_prefix))
xde = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.xdegen.vcf.gz".format(path_prefix))
xtr = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.xtrans.vcf.gz".format(path_prefix))

amp_parse = amp.parse_platypus_VCF(30, 30, 4, "chrY")
het_parse = het.parse_platypus_VCF(30, 30, 4, "chrY")
oth_parse = oth.parse_platypus_VCF(30, 30, 4, "chrY")
par_parse = par.parse_platypus_VCF(30, 30, 4, "chrY")
xde_parse = xde.parse_platypus_VCF(30, 30, 4, "chrY")
xtr_parse = xtr.parse_platypus_VCF(30, 30, 4, "chrY")

xyv.hist_read_balance(
	"chrY", amp_parse[2], "Ampliconic", False, "{}/Ampliconic_chrY".format(
		path_prefix))
xyv.hist_read_balance(
	"chrY", het_parse[2], "Heterochromatic", False, "{}/Heterochromatic_chrY".format(
		path_prefix))
xyv.hist_read_balance(
	"chrY", oth_parse[2], "Other", False, "{}/Other_chrY".format(
		path_prefix))
xyv.hist_read_balance(
	"chrY", par_parse[2], "PAR", False, "{}/PAR_chrY".format(
		path_prefix))
xyv.hist_read_balance(
	"chrY", xde_parse[2], "X-degenerate", False, "{}/X_degen_chrY".format(
		path_prefix))
xyv.hist_read_balance(
	"chrY", xtr_parse[2], "X-transposed", False, "{}/X_trans_chrY".format(
		path_prefix))

# with fixed
xyv.hist_read_balance(
	"chrY", amp_parse[2], "Ampliconic", False, "{}/Ampliconic_chrY_with_fixed".format(
		path_prefix), include_fixed=True)
xyv.hist_read_balance(
	"chrY", het_parse[2], "Heterochromatic", False, "{}/Heterochromatic_chrY_with_fixed".format(
		path_prefix), include_fixed=True)
xyv.hist_read_balance(
	"chrY", oth_parse[2], "Other", False, "{}/Other_chrY_with_fixed".format(
		path_prefix), include_fixed=True)
xyv.hist_read_balance(
	"chrY", par_parse[2], "PAR", False, "{}/PAR_chrY_with_fixed".format(
		path_prefix), include_fixed=True)
xyv.hist_read_balance(
	"chrY", xde_parse[2], "X-degenerate", False, "{}/X_degen_chrY_with_fixed".format(
		path_prefix), include_fixed=True)
xyv.hist_read_balance(
	"chrY", xtr_parse[2], "X-transposed", False, "{}/X_trans_chrY_with_fixed".format(
		path_prefix), include_fixed=True)

parse_list = [amp_parse, het_parse, oth_parse, par_parse, xde_parse, xtr_parse]
parse_list_regions = [
	"ampliconic", "heterochromatic", "other", "par", "xdegen", "xtr"]

with open(out_file, "w") as f:
	f.write("region\tmean_read_balance\tnum_sites\n")
	for idx, i in enumerate(parse_list):
		f.write("{}\t{}\t{}\n".format(
			parse_list_regions[idx], np.mean(i[2]), len(i[2])))
