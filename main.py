"""
# From someone in biostars.org
currentCid = ''
buffer = []

for record in SeqIO.parse(open(sys.argv[1]),"fasta"):
	cid = str(record.description).split('.')[0][1:]

	if currentCid == '':
		currentCid = cid
	else:
		if cid != currentCid:
			buffer.sort(key = lambda x : len(x[1]))
			print('>' + buffer[-1][0])
			print(buffer[-1][1])
			currentCid = cid
			buffer = [(str(record.description),str(record.seq))]
		else:
			buffer.append((str(record.description),str(record.seq)))

buffer.sort(key = lambda x : len(x[1]))
"""
from collections import Counter
from operator import itemgetter
from Bio import SeqIO
import re

file1 = 'dna.example.fasta'
file2 = 'NC_005816.fna'


def SeqIO_parse(file):
    the_record = SeqIO.parse(file, "fasta")
    return the_record


# (1) How many records are in the file?
def num_of_records(file):
    the_record = SeqIO_parse(file)
    print("There are %s records in file '%s'\n" % (len(list(the_record)), file))


# num_of_records('dna.example.fasta')


# (2) What are the lengths of the sequences in the file?
def len_of_seq(file):
    the_record = SeqIO_parse(file)
    print("Here are the lengths of the sequences in file '%s':" % file)
    for record in the_record:
        print(len(record.seq))
    print('\n')


# len_of_seq(file1)


#   What is the longest sequence and what is the shortest sequence?
def longest_n_shortest_seq(file):
    the_record_list = list(SeqIO_parse(file))
    sorted_record = sorted(the_record_list, key=lambda record: len(record.seq))
    print(u"The shortest sequence in the file %s is:%s, with the length of %s\n"
          % (file, sorted_record[0].id, len(sorted_record[0].seq)))
    print("The longest sequence on the other hand is:%s, with the length of %s"
          % (sorted_record[-1].id, len(sorted_record[-1].seq)))


# identify all ORFs present in each sequence of the FASTA file
def identify_all_ORFs(file):
    """
	Finds all ORFs in the given file and returns as a list.
	:param file: File
	:return: List
	"""
    the_record = SeqIO_parse(file)
    ORF_list = []
    for record in the_record:
        ORFs = re.findall(r"ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)", str(
            record.seq))  # I just wish that the website told me that re.findall method returns one group for each CAPTURING GROUP. I should've just fucking looked at the module description. gg
        ORF_list += ORFs
    return ORF_list


# from stackoverflow: max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', record.seq), key = len)


# The longest ORF from identify_ORFs
def longest_ORF(ORF_list):
    """
	Takes a ORF list and finds the longest ORF.
	:param ORF_list: list
	:return: string
	"""
    ORF = max(ORF_list, key=len)
    return ORF


# What is the identifier of the sequence containing the longest ORF?
def ID_of_ORF(ORF, file):
    the_record = SeqIO_parse(file)
    for record in the_record:
        if ORF in record.seq:
            return record.id


# what is the longest ORF contained in the sequence represented by that identifier?

def ORFs_for_ID(ID, file):
    """
	Returns a list of ORFs for a given ID.
	:param ID: The ID to be searched in the file.
	:param file: FILE!
	:return: A list of ORFs for the ID.
	"""
    the_record = SeqIO_parse(file)
    for record in the_record:
        if ID == record.id:
            return re.findall(r"ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)", str(record.seq))


# What is the starting position of the longest ORF in the sequence that contains it?
def pos_ORF(ORF, file):
    """

    :param ORF:
    :type ORF: str
    :param file:
    :type file: file
    :return: Start and end position of the given ORF
    :rtype: print
    """
    the_record = SeqIO_parse(file)
    for record in the_record:
        if ORF in record.seq:
            i = re.search(ORF, str(record.seq))
            print("The ORF %s starts from %s and ends at %s" % (ORF, i.start(), i.end()))


#
# i = identify_all_ORFs(file2)
# i1 = longest_ORF(i)
# i2 = ID_of_ORF(i1, file2)
# print("Longest ORF is ", i1, "\nand its ID is ", i2)
#
# i3 = ORFs_for_ID(i2, file2)
# print("For ID %s, there are %s ORFs and the longest ORF is %s" % (i2, len(i3), longest_ORF(i3)))



# Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. Your program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.

def identify_all_repeats(file, repeat_len) -> list:
    """
    Identifies all overlapping repeats with a given length and returns them as a list.
    :param file:
    :param repeat_len: desired repeat length
    :return: list
    """
    the_record = SeqIO_parse(file)
    repeats = []
    for record in the_record:
        for repeat in range(len(record.seq)):
            if record.seq[repeat:(repeat + repeat_len)] == record.seq[(repeat + repeat_len - 1):(repeat+ 2 * repeat_len - 1)]:
                repeats.append(str(record.seq[repeat:(repeat + repeat_len)]))
    return repeats


def repeats_occurrences(file, repeat_len):
    repeats = identify_all_repeats(file, repeat_len)
    occur = []
    for repeat in set(repeats):
        occur.append([repeat, repeats.count(repeat)])
    print("From the most frequent repeat sequence to the least, in the format [Seq, Freq]:\n", sorted(occur, key= itemgetter(1), reverse=True))
    # return Counter(repeats)-> This is better than the looping but it's 2 ez.

repeats_occurrences(file2,3)