from __future__ import print_function
import argparse
import csv


def parse_args():
	""" Parse command line arguments """
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--bed", required=True,
		help="Input bed file. Fourth column must contain CCDS ID.")

	parser.add_argument(
		"--gene_status", required=True,
		help="Input gene status file. Sixth column contains CCDS ID")

	parser.add_argument(
		"--outfile", required=True,
		help="Output file.")

	parser.add_argument(
		"--missing_ids_outfile", required=True,
		help="Output file for IDs present in input bed but missing from gene"
		"status.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()
	input_ccds_list = []
	with open(args.bed, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			input_ccds_list.append(i[3])

	print("Present in bed file: ", input_ccds_list)
	found_list = []
	missing_in_bed_list = []
	out_list = []
	with open(args.gene_status, "rU") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			if i[5] in input_ccds_list:
				out_list.append(i)
				found_list.append(i[5])
			else:
				missing_in_bed_list.append(i[5])

	print("Found in bed file and gene status file:", found_list)
	print("Not present in bed file:", missing_in_bed_list)

	missing_in_gene_status = [x for x in input_ccds_list if x not in found_list]
	missing_out_list = []
	with open(args.bed, "r") as f:
		reader = csv.reader(f, delimiter="\t")
		for i in reader:
			if i[3] in missing_in_gene_status:
				missing_out_list.append(i)

	with open(args.outfile, "w") as o:
		writer = csv.writer(o, delimiter="\t")
		writer.writerows(out_list)

	with open(args.missing_ids_outfile, "w") as o:
		writer = csv.writer(o, delimiter="\t")
		writer.writerows(missing_out_list)


if __name__ == "__main__":
	main()
