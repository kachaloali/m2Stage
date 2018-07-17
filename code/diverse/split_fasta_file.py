# -*- coding: utf-8 -*-
"""
Created on Tuesday, March 13, 2018, 16:06:49

@author: Kachalo Ali
"""

from Bio import SeqIO
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
   
if __name__ == '__main__':
	import subprocess, commands, os
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	def get_output_name(filename):
		"""
		Returns the file name without its extension.
		"""
		infile_extension = filename.split(".")[-1]
		return filename[:len(filename) - len(infile_extension) - 1]
		
	# Checking arguments
	parser = ArgumentParser(description=__doc__,
		                formatter_class=ArgumentDefaultsHelpFormatter)            	    
	parser.add_argument("-f", "--fasta", dest = "file", type = str, 
						help = "The input fasta file name")
	parser.add_argument("-n", dest = "number", type = int, 
					help = "The number of sequences from which to cut")
	args = parser.parse_args()
	if args.file is None or args.number is None: 
			parser.print_help()
			exit(0)
	
	record_iter = SeqIO.parse(open(args.file),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, args.number)):
		filename = get_output_name(args.file) + "_%s.fasta" % (i + 1)
		with open(filename, "w") as handle:
		    count = SeqIO.write(batch, handle, "fasta")
		print("Wrote %i records to %s" % (count, filename))
