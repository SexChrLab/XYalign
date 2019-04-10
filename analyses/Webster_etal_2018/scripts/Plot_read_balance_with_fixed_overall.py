from __future__ import print_function
from xyalign import variants as xyv
import numpy as np
import sys

out_file = sys.argv[1]
path_prefix = sys.argv[2]

if path_prefix[-1] == "/":
	path_prefix = path_prefix[:-1]

xx = xyv.VCFFile("{}/HG00513_wgs_hg19.noprocessing.vcf.gz".format(path_prefix))
xy = xyv.VCFFile("{}/HG00512_wgs_hg19.noprocessing.vcf.gz".format(path_prefix))

xx_19 = xx.parse_platypus_VCF(30, 30, 4, "chr19")
xx_x = xx.parse_platypus_VCF(30, 30, 4, "chrX")
xy_19 = xy.parse_platypus_VCF(30, 30, 4, "chr19")
xy_x = xy.parse_platypus_VCF(30, 30, 4, "chrX")
xy_y = xy.parse_platypus_VCF(30, 30, 4, "chrY")

xyv.hist_read_balance(
	"chr19", xx_19[2], "HG000513_chr19", False, "{}/HG000513_chr19_with_fixed".format(
		path_prefix), include_fixed=True)

xyv.hist_read_balance(
	"chrX", xx_x[2], "HG000513_chrX", False, "{}/HG000513_chrX_with_fixed".format(
		path_prefix), include_fixed=True)

xyv.hist_read_balance(
	"chr19", xy_19[2], "HG000512_chr19", False, "{}/HG000512_chr19_with_fixed".format(
		path_prefix), include_fixed=True)

xyv.hist_read_balance(
	"chrX", xy_x[2], "HG000512_chrX", False, "{}/HG000512_chrX_with_fixed".format(
		path_prefix), include_fixed=True)

xyv.hist_read_balance(
	"chrY", xy_y[2], "HG000512_chrY", False, "{}/HG000512_chrY_with_fixed".format(
		path_prefix), include_fixed=True)

parse_list = [xx_19, xx_x, xy_19, xy_x, xy_y]
parse_list_names = ["XX_chr19", "XX_chrX", "XY_chr19", "XY_chrX", "XY_chrY"]

with open(out_file, "w") as f:
	f.write("sample_chrom\tmean_read_balance\tnum_sites\n")
	for idx, i in enumerate(parse_list):
		f.write("{}\t{}\t{}\n".format(
			parse_list_names[idx], np.mean(i[2]), len(i[2])))
