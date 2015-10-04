from operator import itemgetter
import re

from Bio import SeqIO

file1 = 'dna.example.fasta'
file2 = 'NC_005816.fna'


def SeqIO_parse(file):
    """
    Turns a fasta file into an iterator.
    :param file: single fasta file
    :type file: str
    :return: iterator
    :rtype: SeqRecords
    """
    the_record = SeqIO.parse(file, "fasta")
    return the_record


def num_of_records(file):
    """
    :param file: single fasta file
    :type file: str
    :return: Prints the # of records (ID and Seq)
    :rtype:
    """
    the_record = SeqIO_parse(file)
    print("There are %s records in file '%s'\n" % (len(list(the_record)), file))


def len_of_seq(file):
    """
    :param file: single fasta file
    :type file: str
    :return: Prints the length of each sequence.
    :rtype:
    """
    the_record = SeqIO_parse(file)
    print("Here are the lengths of the sequences in file '%s':" % file)
    for record in the_record:
        print(len(record.seq))
    print('\n')


def longest_n_shortest_seq(file):
    """
    :param file: single fasta file
    :type file: str
    :return: Prints the shortest and longest sequence in the file.
    :rtype:
    """
    the_record_list = list(SeqIO_parse(file))
    sorted_record = sorted(the_record_list, key=lambda record: len(record.seq))
    print(u"The shortest sequence in the file %s is:%s, with the length of %s\n"
          % (file, sorted_record[0].id, len(sorted_record[0].seq)))
    print("The longest sequence on the other hand is:%s, with the length of %s"
          % (sorted_record[-1].id, len(sorted_record[-1].seq)))


def identify_all_ORFs(file):
    """
	:param file: single fasta file
	:return: a list of all ORFs
	"""
    the_record = SeqIO_parse(file)
    ORF_list = []
    for record in the_record:
        ORFs = re.findall(r"ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)", str(
            record.seq))
        ORF_list += ORFs
    return ORF_list


def longest_ORF(ORF_list):
    """
	:param ORF_list: ORF list
	:return: the longest ORF
	"""
    ORF = max(ORF_list, key=len)
    return ORF


def ID_of_ORF(ORF, file):
    """
    :param ORF: Single ORF
    :type ORF: str
    :param file: single fasta file
    :type file: str
    :return: ID of the ORF
    :rtype: str
    """
    the_record = SeqIO_parse(file)
    for record in the_record:
        if ORF in record.seq:
            return record.id


def ORFs_for_ID(ID, file):
    """
	:param ID: The ID to be searched in the file.
	:param file: single fasta file
	:return: A list of ORFs for the ID.
	"""
    the_record = SeqIO_parse(file)
    for record in the_record:
        if ID == record.id:
            return re.findall(r"ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)", str(record.seq))


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

def identify_all_repeats(file, repeat_len) -> list:
    """
    Identifies all overlapping repeats with a given length and returns them as a list.
    :param file: single fasta file
    :param repeat_len: desired repeat length
    :return: a list of all repeats
    """
    the_record = SeqIO_parse(file)
    repeats = []
    for record in the_record:
        for repeat in range(len(record.seq)):
            if record.seq[repeat:(repeat + repeat_len)] == record.seq[
                                                           (repeat + repeat_len - 1):(repeat + 2 * repeat_len - 1)]:
                repeats.append(str(record.seq[repeat:(repeat + repeat_len)]))
    return repeats


def repeats_occurrences(file, repeat_len):
    """
    :param file:  single fasta file
    :type file: str
    :param repeat_len: desired repeat length
    :type repeat_len: int
    :return: prints a sorted list of repeats
    :rtype:
    """
    repeats = identify_all_repeats(file, repeat_len)
    occur = []
    for repeat in set(repeats):
        occur.append([repeat, repeats.count(repeat)])
    print("From the most frequent repeat sequence to the least, in the format [Seq, Freq]:\n",
          sorted(occur, key=itemgetter(1), reverse=True))
    # return Counter(repeats)-> This is better than the looping but it's 2 ez.


repeats_occurrences(file2, 3)
