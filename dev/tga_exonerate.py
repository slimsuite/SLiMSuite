# -*- coding: utf-8 -*-

"""
Make a python method that will:

Take input attributes qryfas, genome, model
Run exonerate with multiple output formats keeping the top 10 hits:

exonerate eogA_NOTSC__TS10XV2PRIEOG0907027J.nt.fa /srv/slimsuite/data/genomes/apollo/tigersnake.v1.0.fasta --model est2genome --showtargetgff --showcigar --bestn 10

exonerate eogA_NOTSC__TS10XV2PRIEOG0907027J.aa.fa /srv/slimsuite/data/genomes/apollo/tigersnake.v1.0.fasta --model protein2genome --showtargetgff --showcigar --bestn 10

» Do this as os.popen(‘exonerate %s %s --model %s --showtargetgff --showcigar --bestn 10’ % (qryfas, genome, model))

Extract the different outputs into a dictionary, removing line prefixes as required: { query: {‘gff’:[outputlines], ‘cigar’:[outputlines], ‘alignment’:[outputlines], ‘vulgar’:[outputlines] }

Get this working with a single query first and then multiple queries.

Ultimately, I will want to ape GABLAM output tables as much as possible. The next thing will be to convert the vulgar format into a table with headers and convert each line into a dictionary of {header:value} with the header list as the first element of the vulgar list, e.g.

{ ‘vulgar’:[[headerlist], {linedict1}, {linedict2}, …] }

Note: os.popen() is "Deprecated since version 2.6: This function is obsolete. Use the subprocess module. Check especially the Replacing Older Functions with the subprocess Module section." (https://docs.python.org/2/library/os.html#file-object-creation)
Note: subprocess has now been superceeded by subprocess32.  "See also POSIX users (Linux, BSD, etc.) are strongly encouraged to install and use the much more recent subprocess32 module instead of the version included with python 2.7. It is a drop in replacement with better behavior in many situations." (https://docs.python.org/2/library/subprocess.html#module-subprocess)
"""
import string
from subprocess import Popen, PIPE

def exonerate(qryfas, genome, model,exonerate='exonerate',bestn=0,callobj=None):
    #{ query: {‘gff’:[outputlines], ‘cigar’:[outputlines], ‘alignment’:[outputlines], ‘vulgar’:[[headerlist], {header:value}, {header:value}, …] }
    query_dic = {}
    header_list = ['query_id', 'query_start', 'query_end', 'query_strand', 'target_id', 'target_start', 'target_end', 'target_strand', 'score', '<label, query_length, target_length> triplets']
    excmd = [exonerate, qryfas, genome, '--model', model, '--showtargetgff', '--showcigar']
    if bestn: excmd += ['--bestn', '%d' % bestn]
    if callobj: callobj.printLog('#RUN',string.join(excmd))
    process = Popen(excmd, stdout=PIPE)
    output_format = ''
    while True:
        line = process.stdout.readline().rstrip()
        if line:
            if line.startswith('         Query:'):
                query = line.split(':', 1)[1].split(' ')[1]
            if line == 'C4 Alignment:':
                output_format = 'alignment'
            elif line == '# --- START OF GFF DUMP ---':
                output_format = 'gff'
            elif line.startswith('vulgar:'):
                output_format = 'vulgar'
                fields = line.split(' ', 10)[1:]
                if output_format in query_dic[query]:
                    query_dic[query][output_format].append({})
                else:
                    query_dic[query][output_format] = [header_list, {}]
                for header, field in zip(header_list, fields):
                    query_dic[query][output_format][-1][header] = field
            elif line.startswith('cigar:'):
                output_format = 'cigar'
                if output_format in query_dic[query]:
                    query_dic[query][output_format].append(line.replace('cigar: ', ''))
                else:
                    query_dic[query][output_format] = [line.replace('cigar: ', '')]
            elif line == '------------' or line.startswith('Command line:') or line.startswith('Hostname:') or line == '# --- END OF GFF DUMP ---' or line == '#':
                pass
            elif output_format:
                if query in query_dic:
                    if output_format in query_dic[query]:
                        query_dic[query][output_format].append(line)
                    else:
                        query_dic[query][output_format] = [line]
                else:
                    query_dic[query] = {output_format:[line]}
        elif process.poll() is not None:
            break
        elif output_format == 'alignment':
            try: query_dic[query][output_format].append(line)
            except: pass
        if callobj: callobj.bugPrint(line)
    return query_dic

def main(qryfas, genome, model):
    exonerate(qryfas, genome, model)

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "usage: python exonerate.py [options] query.fa genome.fa model"
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error(usage)
    main(args[0], args[1], args[2])
