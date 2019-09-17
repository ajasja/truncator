"""Taken from Ryab Kibler"""

#import json
from collections import Counter
from Bio import SeqIO

from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.Alphabet import IUPAC

import os
import copy


def replace_sequence_in_record(record, location, new_seq):
    #print(record, location, new_seq)
    #print(dir(location))
    #print(location.extract(record.seq))
    #print(record.seq[location.start:location.end])

    if location.strand >= 0:
        adjusted_seq = record.seq[:location.start] + new_seq.seq + record.seq[location.end:]
    else:
        adjusted_seq = record.seq[:location.start] + new_seq.reverse_complement().seq + record.seq[location.end:]

    #exit(adjusted_seq)
    record.seq = adjusted_seq

    #print(help(location))
    #exit(dir(location))

    seq_diff = len(new_seq) - len(location)
    orig_start = location.start
    orig_end = location.end

    processed_features = []

    #print("-=-=-=-=-=-=-=-=-=-==---=-=-=-=-=-=")
    #print(location)
    #print("diff: %d" % seq_diff)

    #adjust all features
    for feat in record.features:
        #print("----------------------------------------")
        #print(feat.qualifiers['label'][0])
        #print(feat.location)

        f_loc = feat.location

        loc_list = []

        for subloc in f_loc.parts:

            assert(subloc.start <= subloc.end)

            #type 1: where the start and end are contained within the original location
            #-> do not add it to the processed_features list
            if subloc.start > location.start and subloc.start < location.end and subloc.end > location.start and subloc.end < location.end:
                #print("start: %d and end: %d are contained within %s" % (subloc.start, subloc.end, location))
                #print("omit")
                continue

            #type 1b: where the start and end are the same which will happen a lot for storing constraints and objectives
            elif subloc.start == location.start and subloc.end == location.end:
                new_loc = FeatureLocation(location.start, location.end + seq_diff, strand=subloc.strand)



            #type 2: where they start or end inside the location
            #-> chop off. don't forget to add on approprate amount
            ##THINK! does strand even matter? How is start and end defined? I'm assuming that for strand -1 things are flipped but that's probably not how it's implemented. Also, consider strand = 0 (no direction). There is probably an easier way.
            elif subloc.start >= location.start and subloc.start <= location.end:
                #print("start: %d is in %s" % (subloc.start, location))
                new_loc = FeatureLocation(location.end + seq_diff, subloc.end + seq_diff, strand=subloc.strand)

            elif subloc.end >= location.start and subloc.end <= location.end:
                #print("end: %d is in %s" % (subloc.end, location))
                new_loc = FeatureLocation(subloc.start, location.start, strand=subloc.strand)


            #type 3: where they span the location
            #-> keep the leftmost point same and add diff to rightmost. do not split
            elif location.start >= subloc.start and location.start <= subloc.end and location.end >= subloc.start and location.end <= subloc.end:
                #print("loc spans insert. keep start and add diff to end")
                new_loc = FeatureLocation(subloc.start, subloc.end + seq_diff, strand=subloc.strand)

            #type 4: where they start and end before location
            #-> add it to list unchanged
            elif subloc.start <= location.start and subloc.end <= location.start:
                #print("loc is before insert location so just keep")
                new_loc = subloc

            #type 5: where they start and end after location
            #-> add diff to whole location
            elif subloc.start >= location.end and subloc.end >= location.end:
                #print("loc is after insert location so apply offset and keep")
                new_loc = subloc + seq_diff

            loc_list.append(new_loc)
            #print("new loc:")
            #print(new_loc)

        #if the list is empty, it means that all the sublocs were contained within the insert
        if loc_list:
            feat.location = sum(loc_list)
            processed_features.append(feat)


    record.features = processed_features


    return record

def find_annotation(record, label):
    for feat in record.features:
        if label == feat.qualifiers['label'][0]:
            #I will be replacing it so remove it:
            #vector.features.remove(feat)
            return feat
            #return None
    raise ValueError("label not found: " + label)

def insert_into_vector(vector, destination_label, new_seq_record):
    """Inserts sequence record into a vector at a feature labled with destination_label """
    vector = copy.deepcopy(vector)

    destination_annotation = find_annotation(vector, destination_label)
    #print(destination_annotation)

    location = destination_annotation.location

    #print(vector)
    #print(vector.features)
    #print(dir(vector.features))
    vector = replace_sequence_in_record(vector, location, new_seq_record)

    #re-annotate the thing
    insert_loc = FeatureLocation(location.start, location.start + len(new_seq_record), strand=location.strand)
    destination_annotation.location = insert_loc
    destination_annotation.qualifiers['label'][0] = new_seq_record.name
    vector.features.append(destination_annotation)

    return vector, insert_loc
