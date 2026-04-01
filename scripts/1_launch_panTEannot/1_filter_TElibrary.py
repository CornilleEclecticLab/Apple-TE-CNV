from Bio import SeqIO
import statistics


f = "../../data/TElibrary/Final_V2_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150.fasta" #"Final_V1_DUAN_DACCORD_ATIP_DATABASE_non-redundant_150sup.fasta"
f_list = "../../data/TElists/250923_TE_structure.txt"
f_list2 = "../../data/TElists/250923_TE_residuals.txt"



def get_lengths( records ):

	d_seq_lengths = {record.id: len(record.seq) for record in records}

	print(f"min = {min(d_seq_lengths.values())} ; mean = {statistics.mean(d_seq_lengths.values())} ; "
	f"median = {statistics.median(d_seq_lengths.values())} ; max = {max(d_seq_lengths.values())} bp\n")



def select_candidate( f_candidate_TEs ):

	ids_to_include = []
	with open(f_candidate_TEs, "r") as f_candidate_TEs:
		for line in f_candidate_TEs:
			ID = line.strip().replace(">", "")
			ids_to_include.append(ID)

	records = [record for record in SeqIO.parse(f, "fasta") if record.id in ids_to_include]

	print("\nSelected TEs:")
	get_lengths(records)

	return records



def filter_TElibrary( records, f_list, extension ):

	output_file = f + "." + extension
	with open(output_file, "w") as out_f:
		SeqIO.write(records, out_f, "fasta")

	print(f"Filtered file saved as: {output_file}")



if __name__ == "__main__":

	## all
	records_all = [record for record in SeqIO.parse(f, "fasta")]
	print("\nAll sequences:")
	get_lengths(records_all)


	## filtered
	filtered_records = select_candidate( f_list )
	filter_TElibrary( filtered_records, f_list, 'TE_structure' )

	filtered_records2 = select_candidate( f_list2 )
	filter_TElibrary( filtered_records2, f_list2, 'TE_residuals' )
