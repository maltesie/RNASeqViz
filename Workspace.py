#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 11:26:50 2021

@author: abc
"""

import numpy as np

def make_fasta(names, sequences):
    my_fasta = ""
    max_len = 0
    for sequence in sequences:
        print(sequence)
    for name, sequence in zip(names, sequences):
        my_fasta += ">{}\n".format(name)
        joined_seq = ''.join(sequence)
        max_len = max(max_len, len(joined_seq))
        my_fasta += "{}\n".format(joined_seq)
    return my_fasta, max_len

def get_inserts(sequence):
    sequence_counter = 0
    insert_counter = 0
    insertions = []
    for nuti in sequence:
        if nuti != '-': 
            if insert_counter != 0:
                insertions.append((sequence_counter, insert_counter))
                insert_counter = 0
            sequence_counter += 1
        else:
            insert_counter += 1
    if insert_counter != 0: insertions.append((sequence_counter, insert_counter))
    return insertions

ins = get_inserts("asdasdsad----asd--asd")
print(ins)

def get_sequence_with_inserts(sequence, inserts):
    my_seq = ""
    if not inserts: return sequence
    my_inserts = inserts[::-1]
    current_insert = my_inserts.pop()
    current_sequence_position = 0
    slice_position = 0
    for i,n in enumerate(sequence):
        if current_sequence_position == current_insert[0]:
            my_seq += sequence[slice_position:i] + '-' * current_insert[1]
            slice_position = i
            if my_inserts: current_insert = my_inserts.pop()
            else: return my_seq + sequence[i:]
        if n != '-': current_sequence_position += 1
    my_seq += sequence[slice_position:]
    if current_insert[0] == current_sequence_position: my_seq += '-' * current_insert[1]
    return my_seq

seq = "asdasdsad----asd--asd"
inserts = [(1,2),(9,3),(15,1)]
new_seq = get_sequence_with_inserts(seq, inserts)
print(new_seq)

def get_sequence_with_absolute_inserts(sequence, inserts):
    my_seq = str(sequence)
    if not inserts: return sequence
    insert_counter = 0
    for insert in inserts:
        if insert[0]+insert_counter > len(my_seq) : return my_seq
        my_seq = my_seq[:insert[0]+insert_counter] + '-' * insert[1] + my_seq[insert[0]+insert_counter:]
        insert_counter += insert[1]
    return my_seq

seq = "asdasdsad----asd--asd"
inserts = [(1,2),(9,3),(15,1)]
new_seq = get_sequence_with_absolute_inserts(seq, inserts)
print(new_seq)

def get_absolute_inserts(sequence, inserts):
    my_inserts = np.array(inserts).reshape(-1,2)
    seq_ins = get_inserts(sequence)
    for seq_in in seq_ins: my_inserts[:,0][my_inserts[:,0] >= seq_in[0]] += seq_in[1]
    return [tuple(seq_in) for seq_in in my_inserts]

seq = "asdasdsad----asd--asd"
inserts = get_inserts(seq)
abs_inserts = get_absolute_inserts(seq, inserts)
print(inserts, abs_inserts)

        
def multialigned(seqname, alignments):
    if not alignments: return make_fasta(['no_data'], [''])
    if len(alignments) == 1: return make_fasta([seqname, alignments[0][0]], alignments[0][1:3])
    spec, seqaln, refaln, score = alignments[0]
    species = [spec]
    refalns = [refaln]
    for current_species, current_seqaln, current_refaln, _ in alignments[1:]:
        seq_inserts = get_inserts(seqaln)
        current_seq_inserts = get_inserts(current_seqaln)
        abs_current_seqaln_inserts = get_absolute_inserts(current_seqaln, seq_inserts)
        abs_new_seqaln_inserts = get_absolute_inserts(seqaln, current_seq_inserts)
        new_refaln = get_sequence_with_absolute_inserts(current_refaln, abs_current_seqaln_inserts)
        seqaln = get_sequence_with_inserts(seqaln, current_seq_inserts)
        for i, aligned_refaln in enumerate(refalns): refalns[i] = get_sequence_with_absolute_inserts(aligned_refaln, abs_new_seqaln_inserts)
        refalns.append(new_refaln)
        species.append(current_species)
    return make_fasta([seqname] + species, [seqaln] + refalns)
    
alignments = [
    ['test1', 'ABCDE-FG', 'ABCDEFFG', 1],
    ['test2', 'ABCDEFG', 'ABC-EFG', 1],
    ['test1', 'ABC-DEFG', 'ABCCDEFG', 1],
    ]

name = 'supertest'

test = multialigned(name, alignments)

print(test)

