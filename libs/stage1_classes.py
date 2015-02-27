#!/usr/bin/env python
# -*- coding: utf-8 -*-

import bio_file_parsers as parser
import stage1_funcs as s1f

#
# Stage 1 classes
#

class Segment:
    # Class to hold info on a single V,D or J segment
    def __init__(self, hit, e, query_start, query_end, hit_start, hit_end,
                 hit_length, identity, align_len):
                     
        self.hit = hit
        self.e = e
        self.q_start = query_start
        self.q_end = query_end
        self.del5 = min(hit_start, hit_end) - 1
        self.del3 = hit_length - max(hit_start, hit_end)
        self.identity = identity
        self.align_len = align_len

class Cluster:
    # Class to hold a single cluster record
    
    def __init__(self, cluster_id, size, seq):
        self.id = cluster_id
        self.size = size
        self.seq = str(seq)
        self.v = []
        self.d = []
        self.j = []
    
    def add_segment(self, blast_record, segment, e_cutoff, top_num):
        """Will send the blast record to the correct function determined by
           the segment provided.
        """
        if segment == 'J':
            self.add_j(blast_record, e_cutoff, top_num)
        elif segment == 'D':
            self.add_d(blast_record, e_cutoff, top_num)
        elif segment == 'V':
            self.add_v(blast_record, e_cutoff, top_num)
    
    def add_vj(self, seg, blast_record, e_cutoff, top_num):
        """Given a BLAST result record for J it will add top_num J segments
           that are below e_cutoff.
        """
        if seg == "J":
            seg_record = self.j
        elif seg == "V":
            seg_record = self.v
        
        # Iterate through hits
        for i in range(len(blast_record['hits'])):
            # Stop if added top_num hits already
            if len(seg_record) == top_num:
                return 0 
            # If e is less then cutoff add hit to record else stop
            e = blast_record['hits'][i]['e']
            if e <= e_cutoff:
                seg_record.append(Segment(
                        blast_record['hits'][i]['hit_id'],
                        blast_record['hits'][i]['e'],
                        blast_record['hits'][i]['qstart'],
                        blast_record['hits'][i]['qend'],
                        blast_record['hits'][i]['sstart'],
                        blast_record['hits'][i]['send'],
                        blast_record['hits'][i]['hit_len'],
                        blast_record['hits'][i]['identity'],
                        blast_record['hits'][i]['align_len']))
            else:
                return 0
                
    def add_d(self, blast_record, e_cutoff, top_num):
        """J segment must already exist, checks that any D hits have 
           E lower than the cutoff and not overlapping the J segment
        """
        if len(self.j) > 0 and len(self.v) > 0:
            # Iterate through hits
            for i in range(len(blast_record['hits'])):
                # Stop if added top_num hits already
                if len(self.d) == top_num:
                    return 0
                # If e is less then cutoff add hit to record and break
                e = blast_record['hits'][i]['e']
                query_start = blast_record['hits'][i]['qstart']
                if e <= e_cutoff:
                    # Check that start of D is beyond end of J and before V
                    if (query_start > self.j[0].q_end and
                        query_start < self.v[0].q_start): 
                        self.d.append(Segment(
                                blast_record['hits'][i]['hit_id'],
                                blast_record['hits'][i]['e'],
                                blast_record['hits'][i]['qstart'],
                                blast_record['hits'][i]['qend'],
                                blast_record['hits'][i]['sstart'],
                                blast_record['hits'][i]['send'],
                                blast_record['hits'][i]['hit_len'],
                                blast_record['hits'][i]['identity'],
                                blast_record['hits'][i]['align_len']))
                else:
                    return 0
        else:
            #~ print 'Warning: no J segment was found.'
            return 1
    
    def get_n_seq(self):
        """If D and J exist then it will calculate the calc the dist 
           between them and save as attribute n
        """
        if len(self.d) > 0 and len(self.j) > 0 and len(self.v) > 0:
            # n1_rev is rev comp of joining region between V and D
            n1_rev = self.seq[self.d[0].q_end:self.v[0].q_start - 1]
            n1 = s1f.reverse_complement_DNA(n1_rev)
            # n2_rev is rev comp of joining region between D and J
            n2_rev = self.seq[self.j[0].q_end:self.d[0].q_start - 1]
            n2 = s1f.reverse_complement_DNA(n2_rev)
            return [str(n1), str(n2)]
        elif len(self.d) == 0 and len(self.j) > 0 and len(self.v) > 0:
            # n1_rev is rev comp of joining region between V and J
            n1_rev = self.seq[self.j[0].q_end:self.v[0].q_start - 1]
            n1 = s1f.reverse_complement_DNA(n1_rev)
            return [str(n1)]
        else:
            return None
            
    def get_top_segment(self, segment):
        """Input is 'D' or 'J' and function returns object containing segment.
        """
        if segment == 'D':
            return str(s1f.first_entry_name(self.d))
        elif segment == 'J':
            return str(s1f.first_entry_name(self.j))
        elif segment == 'V':
            return str(s1f.first_entry_name(self.v))

class Records:
    # Class to hold all cluster records and perform methods on the whole set
    
    def __init__(self):
        self.title = None # Title of the records collection
        self.record_list = [] # To store a list of records ordered by size
        self.id_refs = {} # Needed to associate BLAST results with records
        self.total_clusters_size = 0

    def add_new_record(self, cluster_id, size, seq, sort=True):
        """ Add a new record to the collection, if sort == True then
            records will be added in correct order, if false the list
            should be sorted later.
        """
        new_record = Cluster(cluster_id, size, seq)
        # Add reference to self.id_refs
        self.id_refs[cluster_id] = new_record
        # Add to record_list in place if sort = True, else just append
        if sort:
            added = False
            for i in range(len(self.record_list)):
                if new_record.size > self.record_list[i].size:
                    self.record_list.insert(i, new_record)
                    added = True
                    break
            if added == False:
                self.record_list.append(new_record)
        else:
            self.record_list.append(new_record)
        # Add the size of the record to the total size of all records
        self.total_clusters_size += new_record.size
        return 0
    
    def add_records_from_clustering(self, consensus_list):
        """Will take consensus list produced by clustering and add all records
           to the class.
        """
        for record in consensus_list:
            self.add_new_record(record['id'],
                                record['size'],
                                record['seq'],
                                sort=False)
        # Sort the record list
        self.record_list.sort(key=lambda x: x.size, reverse=True)
        return 0
        
    def add_VDJ_BLAST_xml(self, xml_file, e_cutoff, top_num, segment):
        """Takes BLAST xml output file and directs each blast_record to the
           correct cluster_record so that J then D hits can be added.
        """
        with open(xml_file, 'r') as xml_handle:
            for blast_record in parser.blast_xml_parser(xml_handle):
                # Fetch correct cluster record
                cluster_id = blast_record['query_id']
                cluster_record = self.id_refs[cluster_id]
                # Add J and D hits
                if segment == 'J' or segment == 'V':
                    cluster_record.add_vj(segment, blast_record, e_cutoff,
                                          top_num)
                elif segment == 'D':
                    cluster_record.add_d(blast_record, e_cutoff, top_num)
                
    def write_top_to_fasta(self, top, out_file, size_out=False):
        """Will write the top N sequences (by size) to a fasta file.
           If top=None then all records will be written to fasta.
           If size_out=True then the size of the cluster will be included in
           then fasta header line.
        """
        with open(out_file, 'w') as out_handle:
            for record in self.record_list[:top]:
                # Add size_out to header if True
                if size_out:
                    header = 'centroid={0};size={1}'.format(record.id,
                                                            record.size)
                else:
                    header = record.id
                # Write the seq record to fasta file
                parser.write_fasta(out_handle,
                                   header,
                                   record.seq)
                
        return 0
        
    def write_tabulated(self, out_file, top):
        """Will write a .tab summary file showing the major clones"""
        with open(out_file, 'w') as out_handle:
            # Write header
            header = '\t'.join(['#', 'ID', 'Cluster size', 'Proportion',
                                'V region', 'D region', 'J region',
                                "3' sequence"])
            out_handle.write(header + '\n')
            # Write top N records
            cluster_num = 0
            for record in self.record_list[:top]:
                cluster_num += 1
                proportion = self.get_proportion(record)
                entry = '\t'.join([str(cluster_num),
                                   record.id,
                                   str(record.size),
                                   '{0:.3f}'.format(proportion),
                                   record.get_top_segment('V'),
                                   record.get_top_segment('D'),
                                   record.get_top_segment('J'),
                                   "3'...{0}...5'".format(record.seq)])
                out_handle.write(entry + '\n')
        return 0

    def write_detail_output(self, out_file, top):
        """Will output a detailed summary showing alternative BLAST results,
           and the rearrangement short-hand notation"""
        with open(out_file, 'w') as out_handle:
            cluster_num = 0
            for record in self.record_list[:top]:
                cluster_num += 1
                proportion = self.get_proportion(record)
                # Write header lines
                header = '\n'.join([
                         '> {0} {1}'.format(cluster_num, record.id),
                         '- Reads: {0} Proportion: {1:.3f}'.format(record.size,
                                                                  proportion),
                         ''])
                out_handle.write(header)
                # Write V, D and J hits
                for seg in ['V', 'D', 'J']:
                    out_handle.write('- {0} gene top hits:\n'.format(seg))
                    if seg == 'V':
                        hit_list = record.v
                    elif seg == 'D':
                        hit_list = record.d
                    elif seg == 'J':
                        hit_list = record.j
                    if len(hit_list) == 0:
                        out_handle.write("\tNone\n")
                    hit_num = 0
                    for entry in hit_list:
                        hit_num += 1
                        out_handle.write(s1f.detailed_string(entry, hit_num))
                # Write the shorthand rearrangement
                n_seq = record.get_n_seq()
                # If a D exists (and so there are two N regions)
                if n_seq == None:
                    pass
                elif len(n_seq) == 2:
                    out_handle.write('\n'.join([
                        '- Shorthand rearrangement:',
                        "\t5'~{0}(-{1})...{2}...(-{3}){4}(-{5})...{6}...(-{7}){8}~3'".format(
                                record.v[0].hit,
                                record.v[0].del3,
                                len(n_seq[0]),
                                record.d[0].del5,
                                record.d[0].hit,
                                record.d[0].del3,
                                len(n_seq[1]),
                                record.j[0].del5,
                                record.j[0].hit),
                        "\tN1 sequence: 5'~{0}~3'".format(n_seq[0]),
                        "\tN2 sequence: 5'~{0}~3'\n".format(n_seq[1]),
                        ]))
                # If no D exists (and so there is 1 N region)
                elif len(n_seq) == 1:
                    out_handle.write('\n'.join([
                        '- Shorthand rearrangement:',
                        "\t5'~{0}(-{1})...{2}...(-{3}){4}~3'".format(
                                record.v[0].hit,
                                record.v[0].del3,
                                len(n_seq[0]),
                                record.j[0].del5,
                                record.j[0].hit),
                        "\tN sequence: 5'~{0}~3'\n".format(n_seq[0])
                        ]))
                # Write seq with marker underneath
                out_handle.write('-Sequence domains\n')
                rev_seq = s1f.reverse_complement_DNA(record.seq)
                out_handle.write("\t5'-" + rev_seq + "-3'\n")
                if n_seq != None:
                    if len(n_seq) == 2:
                        marker = 'V' * record.v[0].align_len
                        marker += 'n' * len(n_seq[0])
                        marker += 'D' * record.d[0].align_len
                        marker += 'n' * len(n_seq[1])
                        marker += 'J' * record.j[0].align_len
                    elif len(n_seq) == 1:
                        marker = 'V' * record.v[0].align_len
                        marker += 'n' * len(n_seq[0])
                        marker += 'J' * record.j[0].align_len
                    # Sometimes the marker is shorter than seq as V isnt
                    # aligned at 5' end. Therefore add buffer
                    len_diff = len(rev_seq) - len(marker)
                    if len_diff > 0:
                        marker = ' ' * len_diff + marker
                    out_handle.write("\t5'-" + marker + "-3'\n")
                        
                out_handle.write('\n')

    def get_proportion(self, record):
        """Given a record it will return the total proportion of the reads
           that are contained within that cluster.
        """
        return float(record.size)/self.total_clusters_size

