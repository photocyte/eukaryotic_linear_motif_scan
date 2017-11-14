#! /usr/bin/python 
import Bio
import Bio.SeqIO
import argparse
import re
import urllib2
import StringIO
import os
import csv
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_input',nargs=1,type=str,help='The path to the transdecoder peptide FASTA file')
    parser.add_argument('--elm_string',type=str,help='ELM classes to include.  Default is "TRG" See http://elm.eu.org/downloads.html E.g. "TRG_PTS" for peroxisomal targetting sequences')
    parser.add_argument('--regular_expression',type=str,help='Provide a whole regular expression to search with, which is a much more powerful way to look for specific sequences. Searches this motif instead of downloading motifs from the ELM. See https://en.wikipedia.org/wiki/Regular_expression for information and http://pythex.org to design regular expressions E.g. "FYY" - FYY anywhere in the sequence. Equivalent to "find" in a text editor.  Careful as "*" is a special character with respect to regular expressions, use "\*" instead. E.g. "((FYY\*$)|(FYY$))" - FYY at the sequence C-terminal end with or without an explicit stop codon')
    parser.add_argument('--turn_off_stopcodon_autodetect',action='store_true')
    parser.add_argument('-q',default=True,action='store_false',help='Turn off informational messages on stderr')
    args = parser.parse_args()

def regular_expression_search(fasta_file,elm_string=None,regular_expression=None,verbosity=False,web_object=None):

    header_row=[]
    #Loading FASTA records
    if verbosity: sys.stderr.write("Loading "+fasta_file+"\n")

    if type(fasta_file) is file or isinstance(fasta_file,StringIO.StringIO): ##if the variables type is file
        cached_records = list(Bio.SeqIO.parse(fasta_file,"fasta"))
    elif type(fasta_file) is str and os.path.isfile(fasta_file): ##if there is a physical file
        fasta_file_handle = open(fasta_file,"r")
        cached_records = list(Bio.SeqIO.parse(fasta_file_handle,"fasta"))
    else:
        sys.stderr.write("Can't find file "+fasta_file+"\n")
        return 'error'

    if verbosity: sys.stderr.write("Loaded "+str(len(cached_records))+" fasta records\n")

    if regular_expression == None:
        ##Loading eukaryotic linear motif records, which are just regular expressions.
        if verbosity: sys.stderr.write("Loading eukaryotic linear motif records from http://elm.eu.org...\n")
        if elm_string == None:
            elm_string = "TRG"

        elm_url='http://elm.eu.org/elms/elms_index.tsv?q='+elm_string
        if verbosity: sys.stderr.write("ELM string is "+elm_string+"\n")
        if verbosity: sys.stderr.write("ELM URL is "+elm_url+"\n")
        response = urllib2.urlopen(elm_url)

        html_data = response.read()
        f = StringIO.StringIO(html_data)
        elm_classes = list(csv.reader(f, delimiter='\t'))
        f.close()
        if verbosity: sys.stderr.write("Loaded "+str(len(elm_classes))+" ELM classes\n")
	if verbosity: sys.stderr.write("Going through ELM classes regex with explicit stop codon workaround\n")

	##In case of stop codons being explicit, have to change the regex
	for i in range(0,len(elm_classes)):
            if len(elm_classes[i]) > 2 and not args.turn_off_stopcodon_autodetect:
                elm_classes[i][4] = elm_classes[i][4].replace("$","[*]$")
    else:
        elm_classes = [] 
        header_row=['Accession', 'ELMIdentifier', 'FunctionalSiteName', 'Description', 'Regex', 'Probability', '#Instances', '#Instances_in_PDB']
        if regular_expression != None:
            elm_classes.append(['null','null','null','null',regular_expression,'null','null','null'])


    ##assign header row if not assigned
    if len(header_row) < 4: 
        for row in elm_classes:
            if len(row) > 4:
                header_row = row
                break

    if verbosity: sys.stderr.write("Iterating through records and matching to regular expressions...\n")
    i=0
    j=0

    results_buffer = ''
    if type(web_object) is not None:
        #web_object.write('\t'.join(header_row)+'\n'.replace('\n','<br />')) ##Write directly to the browser!
        #web_object.flush()
	pass
    else: 
        results_buffer += '\t'.join(header_row)+'\n'
	print('\t'.join(header_row)+'\n')
    for record in cached_records:
        i+=1	
        for motif_row in elm_classes:
            if len(motif_row) > 4:
                match = re.search(motif_row[4],str(record.seq))
                if match != None:
                    j+=1
                    output_row = [record.id]+motif_row
                    print('\t'.join(output_row))
                    if type(web_object) is not None:
                        #web_object.write('\t'.join(output_row)+'\n'.replace('\n','<br />')) ##Write directly to the browser!
                        #if j% 5 == 0:
                        #    web_object.flush() ##Flush the buffer w/ every 5th match.
                        pass
		    else:
                        results_buffer += '\t'.join(output_row)+'\n' ##Output as TSV
    if verbosity: sys.stderr.write("Checked "+str(i)+" records and found "+str(j)+" matches\n")
    return results_buffer

if __name__ == "__main__":
    results = regular_expression_search(args.fasta_input[0],elm_string=args.elm_string,regular_expression=args.regular_expression,verbosity=args.q)
    print results
