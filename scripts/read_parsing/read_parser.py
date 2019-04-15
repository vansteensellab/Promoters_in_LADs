#!/usr/bin/env python

# Parse reads (version 1.6)
## Tom van Schaik, 160823

# - Functionality:
# "Hack" cutadapt and use its functions to remove a bunch of adapters and keep
# the information I'm interested in. All other cutadapt functionalies are lost
# (different file input, quality trimming, ...), so these are best done BEFORE
# this read parsing (as this only looks at the sequence, not qualities).
# Based on a read structure file, the read file(s) is/are searched for one
# element at the time, to be written as new read file in the output directory.
# Additional files can be returned, based on the input and options given.

# - Input (See "main" below):
# * Required input:
#  * Read file (.fastq(.gz))
#  * Read structure file (for more details, see this file and the description
#    given there)
#  * Output directory
#
# * Optional input:
#  * Paired read file if paired end (default None)
#  * Basename to be used for ALL output files (default read file basename)
#  * Quiet - do not give execution information (default False)
#  * All - return an "all file", with all elements found for each read (default
#    None)
#  * Number - number of reads written to all file (default 1000)
#  * Barcode sequence - if barcodes are given, also write original (complete)
#    sequences to barcode file
#  * Max error rate - to be used for adapter searches
#  * Min overlap - to be used for adapter searches
#  * Min length - minimum length of reads after trimming

# - Output:
# * Standard output:
#  * Trimmed read files (= output_dir/basename.fastq.gz)
#  * Trimming statistics - information on elements found and final reads written
#    (tab delimeted) (= output_dir/basename.statistics.txt)
#
# * Optional output:
#  * All file - file with all elements found for each read (tab delimeted)
#    (= output_dir/basename.all.txt.gz)
#  * Barcode file - file with read IDs linked to barcode sequences. Optionally
#    include original read sequences (tab delimeted)
#    (= output_dir/basename.barcode.txt.gz)
#  * DNA file - file with fastq sequences for additional DNA elements (fastq
#    format) (= output_dir/basename_[ID].fastq.gz)
#  * [Write the trimmed reads to stdout]
#
# * Variations:
#  * Paired - if paired, modify returned read name to
#    output_dir/basename.1(/2).fastq.gz (note the ".")
#  * Index - if index is given, separate reads on this and put them in files
#    output_dir/basename_[index].fastq.gz

# - Method (briefly):
# Note that this is basically a "hacked" version of cutadapt. This means that
# it uses all great tools written by them. However, this also means that it
# has to follow their logic, which can be quite unlogical. In summary:
# 1) Using the read matrix, input and output will be determined and adapter
#    modifiers will be created
# 2) For each read, each adapter will be removed and information saved
# 3) + All required files will be appended
# 4) Finally, get all important statistics and report these

# - Versions (+ changes):
# * 1.0 - initial version (with cutadapt version 1.11)
# * 1.1 - all-file now contains "-" for not found, and "NA" for skipped,
#         added max_error_rate option
# * 1.2 - added min_overlap option, added "group"-requirement possibility
#         (filter for one of multiple features), extended read_structure.txt
#         file)
# * 1.3 - added min_length option. Partially implemented a better indexing
#         strategy
# * 1.4 - changed barcode name, added min_length statistic and bugs
# * 1.5 - fixed reading indices from file, fix n_tooshort statistic
# * 1.6 - added index to barcode file

# - To do:
# * Allow stdout output for PE data(?)
#   (Not very important in my opinion, as the interleaved output is not ideal
#   anyway)
# * Build more checks, as whether the read_file exists before continuing too far
#   (Not very important in my opinion, although very useful)
# * Incorporate index in statistics
# * Merge more functions between single and paired end read processing, instead
#   of this copy/pasting + some modifications
#   (Not very important in my opinion, although it would make things more
#    structured)
# * Remove unnecessary "filters" from the code
#   (Not very important in my opinion)
#(* This could be organized better and commented more clearly)
#   (Not very important in my opinion, although useful)
# * Extract variable part of const_xxx with sequence alignment instead of 2x
#   cutadapt
#   (Not very important in my opinion, although it would prevent the const_xxx
#    problem!)

from __future__ import print_function, division, absolute_import

import sys, os, optparse, itertools
import logging
import re
# import new
from Bio import Seq
import editdistance

###############################################################################
##### Cutadapt main script starts here - heavily tweaked ######################
###############################################################################

# Print a helpful error message if the extension modules cannot be imported.
from cutadapt import check_importability
check_importability()

import time
import errno
import functools
import platform
import textwrap

from cutadapt import seqio, __version__
from cutadapt.xopen import xopen
from cutadapt.adapters import AdapterParser, Adapter
from cutadapt.modifiers import (LengthTagModifier, SuffixRemover, PrefixSuffixAdder,
	DoubleEncoder, ZeroCapper, PrimerTrimmer, QualityTrimmer, UnconditionalCutter,
	NEndTrimmer, AdapterCutter, NextseqQualityTrimmer)
from cutadapt.filters import (NoFilter, PairedNoFilter, Redirector, PairedRedirector,
	LegacyPairedRedirector, TooShortReadFilter, TooLongReadFilter,
	Demultiplexer, NContentFilter, DiscardUntrimmedFilter, DiscardTrimmedFilter)
from cutadapt.report import Statistics, print_report, redirect_standard_output
from cutadapt.compat import next

def process_single_reads(reader, modifiers, filters, all_file, all_number,
                         barcode_file, has_index=False, barcode_seq=False,
                         DNAs=None, bc_named=False, min_length=None):
   """
   Cutadapt process of single reads, modified.

   Description:
   This is a modified version of the cutadapt process_single_reads. Basically,
   reads are parsed into cutadapt "Sequence"-classes, after which each "modifier"
   (adapter/...) is removed (if found) and this information stored. Finally,
   all required files are written.

   Input:
   - reader: object (cutadapt)
   - modifiers: list (instance = modifier (cutadapt), ID, type, "req", "second",
                      "keep_bases", "const_bar")
   - filters: list (cutadapt, index sequence or None) - just the "write to file"
              filter here
   - all_file: all file to write to, or None
   - all_number: if all_file, number of reads
   - barcode_file: barcode file to write to, or None
   - has_index: index present
   - barcode_seq: report sequences in barcode_file

   Output:
   - Statistics object (cutadapt)
   + fastq file (+ barcode / all_file)

   Method:
   Loop over reads, find adapters and trim reads, apply modifiers and
   output modified reads. If a cutting requirement is not met, stop the process
   but do write the all_file
   """
   n = 0             # no. of processed reads
   n_written = 0     # no. of written reads
   total_bp = 0
   n_tooshort = 0

   barcode_IDs = []

   for read in reader:
      all_info = []
      barcodes = []
      DNA_hits = {}
      req_groups = {}

      n += 1
      total_bp += len(read.sequence)
      keep = True
      seq = read.sequence

      for modifier, _ID, _type, _req, _second, keep_bases, const_var in modifiers:

         # Remember barcodes IDs
         if n == 1 and _type in ("barcode", "const_bar") and bc_named:
            barcode_IDs.append(_ID)

         # If a "group-requirement" is given (one of the adapters of a group is
         # present), add this to the group dictionary
         try:
            if _req.startswith("group:") and _req not in req_groups.keys():
               req_groups[_req] = False
         except AttributeError:
            pass

         # Only continue when reads are still interesting
         if keep:
            read = modifier(read)

            # this is the found match
            match = read.match_info[0][5] if read.match_info else "-"

            if all_file:
               all_info.append(match)

            if _type in ("barcode", "const_bar"):
               if _type == "barcode":
                  barcodes.append(match)
               else:
                  # For the const_bar method, find the barcode within the match
                  barcode = find_var(read, const_var)
                  if barcode == "" and _req == True:
                     keep = False
                  barcodes.append(barcode)
                  if all_file:
                     all_info.append(barcode)

            if match and match != "-":
               try:
                  if _req == False:
                     keep = False
                  elif _req.startswith("group:"):
                     # If a match is found, set the group dictionary to True
                     req_groups[_req] = True
               except AttributeError:
                  pass

               # Let's only write fastq entries when there is something to write
               if _type in ("DNA", "const_DNA"):
                  if _type == "DNA":
                     qual = read.match_info[0][9]
                     dna = seqio.Sequence(read.name, match, qual)
                  else:
                     dna = find_var(read, const_var, "sequence")
                     if dna.sequence == "" and _req == True:
                        keep = False
                  # Write DNA sequences to file using a "just write" filter
                  DNA_hits[_ID] = dna
                  if all_file:
                     all_info.append(dna.sequence)

               if _type == "index":
                  # This is the sequence matched to, not the match itself
                  # Note: there can only be one index to prevent too complex output files
                  # (for now?)
                  index = modifier.adapter

               if keep_bases:
                  # Keep n bases from the match removed
                  if keep_bases[1] == "_5":
                     read.sequence = match[-keep_bases[0]:] + read.sequence
                     read.qualities = read.match_info[0][9][-keep_bases[0]:] + read.qualities
                  else:
                     read.sequence = read.sequence + match[:keep_bases[0]]
                     read.qualities = read.qualities + read.match_info[0][9][:keep_bases[0]]

            else:
               if _req == True:
                  keep = False
               # If req is "-", do nothing and just continue

         else:
            all_info.append("NA")

      # Read trimming done, now append findings to files
      if all_file:
         if n <= all_number:
            print("\t".join([read.name, seq] + all_info), file=all_file)

      # Filter reads where not all group requirements are True
      if not all(req_groups.values()):
         keep = False

      # Filter reads on minimum length
      if min_length:
         if len(read) < min_length:
            n_tooshort += 1
            keep = False

      if keep:
         # Only report barcodes of good reads
         if barcodes:
            if has_index == False:
               if barcode_seq:
                  print("\t".join([read.name, read.sequence] + barcodes), file=barcode_file)
               else:
                  print("\t".join([read.name] + barcodes), file=barcode_file)
            else:
               if barcode_seq:
                  print("\t".join([read.name, read.sequence] + barcodes + [index]), file=barcode_file)
               else:
                  print("\t".join([read.name] + barcodes + [index]), file=barcode_file)

            if bc_named:
               for i, _ID in enumerate(barcode_IDs):
                  # read.name += "|" + _ID + ":" + barcodes[i]
                  read.name = barcode_in_readname(read.name, _ID, barcodes[i])

         # Only report DNA hits of "good" reads
         if DNA_hits:
            for _ID in DNA_hits:
               DNAs[_ID](DNA_hits[_ID])

         # Old cutadapt style - first remove all bits from the read, then put reads
         # through filters. But redundant, as I've removed all the extra features
         # for which this would be useful
         for filter, _seq in filters:
            if has_index == False:
               # This is the cutadapt-way; only write when the filter is passed
               # and then "break"
               if filter(read):
                  break
            else:
               if index == _seq:
                  if filter(read):
                     break

         n_written += 1

   return (Statistics(n=n, total_bp1=total_bp, total_bp2=None), n_written, n_tooshort)

def process_paired_reads(paired_reader, modifiers, filters, all_file, all_number,
                         barcode_file, has_index=False, barcode_seq=False,
                         DNAs=None, bc_named=False, min_length=None):
   """
   Cutadapt process of paired reads, modified.

   Description:
   This is a modified version of the cutadapt process_paired_reads. Basically,
   reads are parsed into cutadapt "Sequence"-classes, after which each "modifier"
   (adapter/...) is removed (if found) and this information stored. Finally,
   all required files are written.


   Input:
   - paired_reader: object (cutadapt)
   - modifiers: list (instance = modifier (cutadapt), ID, type, "req", "second",
                      "keep_bases", "const_bar")
   - filters: list (cutadapt) - just the "write to file" filter
   - all_file: all file to write to, or None
   - all_number: if all_file, number of reads
   - barcode_file: barcode file to write to, or None
   - has_index: index present
   - barcode_seq: report sequences in barcode_file

   Output:
   - Statistics object (cutadapt)
   + two fastq files (+ barcode / all_file)

   Method:
   Loop over reads, find adapters and trim reads, apply modifiers and
   output modified reads. If a cutting requirement is not met, stop the process
   but do write the all_file
   """

   n = 0             # no. of processed reads
   n_written = 0     # no. of written reads
   total1_bp = 0
   total2_bp = 0
   n_tooshort = 0

   # Keep track of barcode IDs for later (this can probably be done when creating
   # adapters as well)
   barcode_IDs = []

   # Keep track of the barcode starts and lengths for "const_bar_comp"
   bc_starts = {}
   bc_len = {}

   for read1, read2 in paired_reader:
      all_info = []
      barcodes = []
      DNA_hits = {}
      req_groups = {}

      n += 1
      total1_bp += len(read1.sequence)
      total2_bp += len(read2.sequence)

      keep = True
      seq1 = read1.sequence
      seq2 = read2.sequence

      for modifier, _ID, _type, _req, _second, keep_bases, const_var in modifiers:

         if n == 1:
            # Remember barcodes IDs
            if _type in ("barcode", "const_bar") and bc_named:
               barcode_IDs.append(_ID)

         # If a "group-requirement" is given (one of the adapters of a group is
         # present), add this to the group dictionary
         try:
            if _req.startswith("group:") and _req not in req_groups.keys():
               req_groups[_req] = False
         except AttributeError:
            pass

         # Only continue with good reads
         if keep:

            # "Cheesy" way to select the correct read for each modifier
            # In this way, you keep the modifier order and stop when a requirement
            # (in either of the reads) is not met.
            if _second:
               read = read2
            else:
               read = read1

            # Create a read-specific complementary barcode-sequence Adapter
            if _type == 'const_bar_comp':

               # Get adapter and "current" modifier sequence, using the last barcode found
               bar_comp = str(Seq.Seq(barcodes[-1]).reverse_complement())
               old_seq = modifier.adapters[0].sequence

               # Replace barcode in modifier sequence with new sequence
               if _ID not in bc_starts:
                  # At first, replace "barcode" with start + 4 ([BC])
                  bc_starts[_ID] = re.search('\[BC\]', old_seq).start()
                  new_seq = old_seq.replace(old_seq[bc_starts[_ID]:(bc_starts[_ID] + 4)], bar_comp)
               else:
                  # Then, replace "barcode" with start + length previous barcode
                  new_seq = old_seq.replace(old_seq[bc_starts[_ID]:(bc_starts[_ID] + bc_len[_ID])], bar_comp)
               bc_len[_ID] = len(bar_comp)

               # Replace Adapter with new Adapter
               modifier.adapters[0] = Adapter(new_seq, where = modifier.adapters[0].where)

            read = modifier(read)
            match = read.match_info[0][5] if read.match_info else "-"

            if _type in ("barcode", "const_bar"):
               if _type == "barcode":
                  barcodes.append(match)
               else:
                  barcode = find_var(read, const_var)
                  barcodes.append(barcode)

            if all_file:
               all_info.append(match)

            if match and match != "-":
               try:
                  if _req == False:
                     keep = False
                  elif _req.startswith("group:"):
                     # If a match is found, set the group dictionary to True
                     req_groups[_req] = True
               except AttributeError:
                  pass

               # Only write DNA fastq entries of matched reads
               if _type in ("DNA", "const_DNA"):
                  if _type == "DNA":
                     qual = read.match_info[0][9]
                     dna = seqio.Sequence(read.name, match, qual)
                  else:
                     dna = find_var(read, const_var, "sequence")
                  DNA_hits[_ID] = dna

               if _type == "index":
                  # This is the sequence matched to, not the match itself
                  # Note: there can only be one index to prevent too complex output files
                  # (for now?)
                  index = modifier.adapter

               if keep_bases:
                  # Keep n bases from the match removed
                  if keep_bases[1] == "_5":
                     read.sequence = match[-keep_bases[0]:] + read.sequence
                     read.qualities = read.match_info[0][9][-keep_bases[0]:] + read.qualities
                  else:
                     read.sequence = read.sequence + match[:keep_bases[0]]
                     read.qualities = read.qualities + read.match_info[0][9][:keep_bases[0]]

            else:
               if _req == True:
                  keep = False
               # If req is "-", do nothing and just continue

            # Cheesy paired end method end, ready for the next modifier
            if _second:
               read2 = read
            else:
               read1 = read

         else:
            all_info.append("NA")

      # Read trimming done, now append findings to file
      if all_file:
         if n <= all_number:
            print("\t".join([read1.name, seq1, seq2] + all_info), file=all_file)

      # Filter reads where not all group requirements are True
      if not all(req_groups.values()):
         keep = False

      # Filter reads on minimum length
      if min_length:
         if len(read1) < min_length or len(read2) < min_length:
            n_tooshort += 1
            keep = False

      if keep:
         if barcodes:
            if has_index == False:
               if barcode_seq:
                  print("\t".join([read.name, read1.sequence, read2.sequence] + barcodes),
                        file=barcode_file)
               else:
                  print("\t".join([read.name] + barcodes), file=barcode_file)
            else:
               if barcode_seq:
                  print("\t".join([read.name, read1.sequence, read2.sequence] + barcodes + [index]),
                        file=barcode_file)
               else:
                  print("\t".join([read.name] + barcodes + [index]), file=barcode_file)

            if bc_named:
               for i, _ID in enumerate(barcode_IDs):
                  # read1.name += "|" + _ID + ":" + barcodes[i]
                  # read2.name += "|" + _ID + ":" + barcodes[i]
                  read1.name = barcode_in_readname(read1.name, _ID, barcodes[i])
                  read2.name = barcode_in_readname(read2.name, _ID, barcodes[i])


         # Only report DNA hits of "good" reads
         if DNA_hits:
            for _ID in DNA_hits:
               DNAs[_ID](DNA_hits[_ID])

         # Old cutadapt style - first remove all bits from the read, then put reads
         # through filters. But redundant, as I've removed all the extra features
         # for which this would be useful
         for filter, _seq in filters:
            if has_index == False:
               if filter(read1, read2):
                  break
            else:
               if index == _seq:
                  if filter(read1, read2):
                     break

         n_written += 1

   return (Statistics(n=n, total_bp1=total1_bp, total_bp2=total2_bp), n_written, n_tooshort)

def run(reads, read_matrix, output_base,
        all_file=False,
        all_number=1000,
        barcode_seq=False,
        stdout=False,
        bc_named=False,
        file_format="fastq",
        max_error_rate=0.1,
   	  min_overlap=3,
        min_length=None,
   	  match_read_wildcards=False,
   	  match_adapter_wildcards=True,
   	  indels=True):
   """
   Main function of cutadapt, heavily modified

   Description:
   This is a modified version of the cutadapt "run". This encompasses creating
   all required objects, starting the read-trimming and writing the statistics

   Input:
   - reads: read file name (.fastq / .fastq.gz)
   - read_matrix: read structure matrix as described elsewhere
   - output_base: basename of all written files
   - all_file: all-file for all found adapters / sequences
   - all_number: number of reads in the all-file
   - barcode_seq: if a barcode is present, return the read sequence(s) in this
     file
   - stdout: write trimmed reads to stdout (only single end)
   - bc-named: add barcode sequence to read name
   - ... (argument names make sense)

   Output:
   Written files:
   1) fastq files at output_base + "fastq.gz" (this could be stdout)
   2) (Optional:) all_file output_base + ".all_file.txt.gz"
   3) (Optional:) barcode_file output_base + ".barcode.txt.gz"
   4) statistics file at output_base + ".statistics.txt"

   Method:
   - Setup input / output file names
   - Parse adapters
   - Setup output fastq file name (using index information)
   - Process reads
   - Report statistics
   """

   logging.debug("Trimming with (a modified) cutadapt - version: %s\n" % __version__)

   ############################################################################
   ### Part 1: create AdapterCutters and organize input / output files ########
   ############################################################################

   if reads[1]:
      logging.debug("Paired end detected")
      paired = True
   else:
      paired = False

   if stdout:
      assert paired == False, "No stdout output for paired end"

   # Input file setup
   input_filename = reads[0]
   if paired:
      input_paired_filename = reads[1]
   else:
      input_paired_filename = None

   # Cutadapt methods for reading and writing files
   try:
   	reader = seqio.open(input_filename, file2=input_paired_filename,
   			qualfile=None, colorspace=False,
   			fileformat=file_format, interleaved=False)
   except (seqio.UnknownFileType, IOError) as e:
   	raise

   # Open writer(s)
   open_writer = functools.partial(seqio.open, mode='w',
   	qualities=reader.delivers_qualities, colorspace=False,
   	interleaved=False)

   # Set-up "all-file"
   if all_file is not None:
   	all_file = xopen(output_base + ".all.txt.gz", 'w')

   # Set-up adapter parser(s) and create adapters from read_structure
   adapter_parser = AdapterParser(
   	colorspace=False,
   	max_error_rate=max_error_rate,
   	min_overlap=min_overlap,
   	read_wildcards=match_read_wildcards,
   	adapter_wildcards=match_adapter_wildcards,
   	indels=indels)

   # Tweak a method in the AdapterCutter class to save the adapter
   # "reference" sequence found
   AdapterCutter._best_match = best_match_index

   # Information objects for adapter parsing:
   modifiers = []       # List of AdapterCutter objects
   all_IDs = []         # List of IDs found for the all_file header
   barcode_file = None  # If barcodes present, file to write to
   barcode_IDs = []     # List of barcodes for the barcode_file header
   DNAs = {}            # Dictionary of DNA elements ID to writers
   has_index = False    # Value to store whether an index is present (max 1)

   # Create adapters from the read matrix
   features = read_features(read_matrix)
   adapters_list = create_adapters(features, adapter_parser)

   # Loop over multiple adapter instances and save them as "modifiers"
   for adapters, _ID, _type, _req, _second, keep_bases, const_bar, index_names in adapters_list:

      adapter_cutter = AdapterCutter(adapters=adapters,
                                     times=1,
                                     wildcard_file=None,
                                     info_file=None,
                                     rest_writer=None,
                                     action="trim")
      # "Manually" set to store match info
      adapter_cutter.keep_match_info = True

      modifiers.append((adapter_cutter, _ID, _type, _req, _second, keep_bases, const_bar))

      if all_file:
         all_IDs.append(_ID)

      if _type == "index":
         # Store that an index has been found and determine output_bases for each
         # index
         assert has_index == False, "One index allowed for now"
         has_index = True

         output_index = []
         for i, adapter in enumerate(adapters):
            seq = adapter.sequence
            # If index names are present, replace the output basename with this
            # index basename. Note that the statistisc file and everything are
            # still written with the "default" basename!
            if index_names:
               index_base = os.path.join(os.path.dirname(output_base), index_names[i])
               output_index.append((index_base, seq))
            else:
               output_index.append((output_base + "_" + seq, seq))

      if _type in ("barcode", "const_bar"):
         # If a barcode, set-up the barcode file
         if barcode_file == None:
            barcode_file = xopen(output_base + ".barcode.txt.gz", "w")
         barcode_IDs.append(_ID)
         if _type == "const_bar" and all_file:
            all_IDs.append(_ID + "_barcode")

      if _type in ("DNA", "const_DNA"):
         # If a DNA sequence, set-up the DNA file
         wr = open_writer("%s_%s.fastq.gz" % (output_base, _ID),
                          None)
         DNAs[_ID] = NoFilter(wr)
         if _type == "const_DNA" and all_file:
            all_IDs.append(_ID + "_DNA")

      if _type == "const_bar_comp":
         # If a complementary barcode, make sure a barcode is already present!
         if not barcode_IDs:
            raise BaseException("Cannot have a const_bar_comp feature without a barcode first!")

   # Write headers for the all_file and barcode_file (if present)
   if all_file:
      if paired:
         print("\t".join(["ID", "sequence_1", "sequence_2"] + all_IDs), file=all_file)
      else:
         print("\t".join(["ID", "sequence"] + all_IDs), file=all_file)

   # Set-up the header of the barcode file. Note that I add an "index" when
   # there is an index.
   if barcode_file:
      if has_index == False:
         if barcode_seq:
            if paired:
               print("\t".join(["ID", "sequence_1", "sequence_2"] + barcode_IDs),
                     file=barcode_file)
            else:
               print("\t".join(["ID", "sequence"] + barcode_IDs),
                     file=barcode_file)
         else:
            print("\t".join(["ID"] + barcode_IDs), file=barcode_file)
      else:
         if barcode_seq:
            if paired:
               print("\t".join(["ID", "sequence_1", "sequence_2"] + barcode_IDs + ["index"]),
                     file=barcode_file)
            else:
               print("\t".join(["ID", "sequence"] + barcode_IDs + ["index"]),
                     file=barcode_file)
         else:
            print("\t".join(["ID"] + barcode_IDs + ["index"]), file=barcode_file)

   # Set-up filters
   # Note: filters are "old"-cutadapt style. Basically, a list of filters is
   # created for each writer-object. Whenever a read passes a filter, it will
   # be written using the writer-object, and using a "break" the writing process
   # will be stopped. Instead of rewriting this part, I removed all filters
   # except the "NoFilter" filter (which is basically no filter, but just
   # writing). This feels a bit confusing, and can be rewritten
   filters = []

   if has_index == False:
      if paired:
         output_file = output_base + ".1.fastq.gz"
         output_paired_file = output_base + ".2.fastq.gz"
      else:
         output_file = output_base + ".fastq.gz"
         output_paired_file = None

      if stdout:
         writer = open_writer(sys.stdout)
      else:
         writer = open_writer(output_file, output_paired_file)

      # Append a tuple with the filter and the "index" sequence (in this case
      # always None)
      if not paired:
      	filters.append((NoFilter(writer), None))
      else:
      	filters.append((PairedNoFilter(writer), None))

   else:
      if stdout:
         logging.debug("Sorry, no stdout with indices involved")

      for index, _seq in output_index:
         if paired:
            output_file = index + ".1.fastq.gz"
            output_paired_file = index + ".2.fastq.gz"
         else:
            output_file = index + ".fastq.gz"
            output_paired_file = None

         writer = open_writer(output_file, output_paired_file)

         # Append a tuple with the filter and the index sequence. This will
         # be later used to match against
         if not paired:
         	filters.append((NoFilter(writer), _seq))
         else:
         	filters.append((PairedNoFilter(writer), _seq))

   ############################################################################
   ### Part 2: process each read and write output #############################
   ############################################################################

   # Main read processing (as done in cutadapt)
   start_time = time.clock()
   try:
   	if paired:
         stats = process_paired_reads(reader, modifiers, filters, all_file,
                                      all_number, barcode_file, has_index,
                                      barcode_seq, DNAs, bc_named, min_length)
   	else:
         stats = process_single_reads(reader, modifiers, filters, all_file,
                                      all_number, barcode_file, has_index,
                                      barcode_seq, DNAs, bc_named, min_length)
   except KeyboardInterrupt as e:
   	print("Interrupted", file=sys.stderr)
   	sys.exit(130)
   except IOError as e:
   	if e.errno == errno.EPIPE:
   		sys.exit(1)
   	raise
   except (seqio.FormatError, EOFError) as e:
   	sys.exit("cutadapt: error: {0}".format(e))

   ############################################################################
   ### Part 3: write statistics ###############################################
   ############################################################################

   # Prepare and report statistics
   elapsed_time = time.clock() - start_time

   statistics_file = xopen(output_base + ".statistics.txt", "w")

   # Note that I have to select the first element in each list to get the
   # actual adapter / modifier / filter
   stats[0].collect([x[0] for x in adapters_list], elapsed_time,
                    [x[0] for x in modifiers], [], [x[0] for x in filters])

   # And use a custom function to write statistics to
   stats = gather_statistics(stats, modifiers, statistics_file)

   ############################################################################
   ### Finally: close open files ##############################################
   ############################################################################

   # close open files
   for f in [writer, all_file, barcode_file, statistics_file]:
      if f is not None and f is not sys.stdin and f is not sys.stdout:
         f.close()

   for f in DNAs.values():
      f.writer.close()

   logging.debug("Closed all open files")
   logging.debug("Done.\n")

###############################################################################
##### Additional (supporting) functions #######################################
###############################################################################

def setup_logging(quiet, filename=False):
   """
   Small function for logging
   """

   if filename:
      logging.basicConfig(filename=filename, level=logging.INFO, format='%(message)s')
   else:
      logging.basicConfig(level=logging.INFO, format='%(message)s')

   if not quiet:
      logging.getLogger().setLevel(logging.DEBUG)

def nonblank_lines(f):
   """
   Filter out empty lines from file
   """
   for l in f:
      line = l.rstrip()
      if line:
         yield line

def read_features(structure):
   """
   Loop through the structure file, saving all variables on the go
   Discard the header
   """
   features = []

   logging.debug("-------\nFeatures:")
   with open(structure, "r") as input:
      for line in nonblank_lines(input):
         if line:
            if line.startswith("#") == True:
               continue
            line = line.split()
            logging.debug("\t".join(line))
            features.append(line)
   logging.debug("-------\n")

   # The first feature is the "description-line", which shouldn't be returned
   return features[1:]

def process_index(indices, max_error_rate, override = True):
   """
   This function checks indices for in-between distances and separates a possible
   name from the given sequence
   """

   sequences = []
   index_names = []

   # Separate sequence from names, using the semicolon
   i = 0
   for index in indices:
      index = index.split(":")
      if len(index) == 2:
         index_names.append(index[0])
         sequences.append(index[1])
      else:
         sequences.append(index[0])
      i += 1
   if len(index_names) != 0:
      assert len(sequences) == len(index_names) == i, \
         "Not all indices can be separated by name, are all ':'-separated?"

   # Check whether the maximum error rate allows for proper separation,
   # implemented using the simple method:
   # - Determine the minimum levenshtein distance
   # - Determine the shortest index
   # - Determine whether these are separable, and report this
   shortest = len(min(sequences))
   min_distance = len(max(sequences))

   for i in sequences:
      for j in sequences:
         if i == j:
            pass
         else:
            distance = editdistance.eval(i, j)
            if distance < min_distance:
               min_distance = distance

   # Report a warning if the max_error_rate and shortest length can lead to
   # "unknown" calls (two indices could be present)
   errors = (min_distance-1)/2
   if int(shortest * max_error_rate) > errors:
      logging.debug("With this max_error_rate, indices with errors can overlap!")

   # If override is true, set the maximum_error_rate "as high as possible"
   # Note: for now, only decrease the maximum error rate, not increase
   max_error_rate_new = max_error_rate
   if override == True:
#      if int(shortest * max_error_rate) != errors:
      if int(shortest * max_error_rate) > errors:
         max_error_rate_new = float(errors) / shortest + 0.01
         logging.debug("Overriding max_error_rate to %f, which is %d error(s)" % (
            max_error_rate_new, max_error_rate_new * shortest))

   return sequences, index_names, max_error_rate_new

def process_end(read_end, which, _ID, _type, _pos,
                max_error_rate, adapter_parser, const_var = None):
   """
   Modify given read ends to suitable adapter_parser input
   """

   index_names = []

   if read_end == "-":
      read_end = []
   else:
      # If a txt file is given, read this for the sequences
      if read_end.endswith(".txt"):
         assert os.path.isfile(read_end), "%s is not a file" % read_end
         with open(read_end, "r") as f:
            read_end = f.read().replace("\n", ",")

      # Separate (possible) multiple sequences
      read_end = read_end.strip(",")
      read_end = read_end.split(",")

      # If integer is given, convert this to N-sequence
      try:
         read_end = ["N"*int(x) for x in read_end]
      except ValueError:
         pass

      # If the element is a const_xxx, remove brackets and create linkedadapter
      # cutter
      # + Update max_error_rate depending on size constant region / total region
      if _type in ("const_bar", "const_DNA"):

         # Let's keep it at one of these elements
         assert len(read_end) == 1, "only one of this const_xxx things please"

         # Extract barcode location (and replace with Ns if int given)
         end = re.split(r'{|}', read_end[0])
         try:
            end[1] = "N" * int(end[1])
         except ValueError:
            pass

         # Create a second (linked) parser to extract the barcode
         # Note: this is more lenient with errors as this is only executed when
         # the complete const_bar type is found, so the barcode MUST be there
         adapter_parser.constructor_args["max_error_rate"] = 0.2
         const_var = {5 : None, 3 : None}
         if end[0]:
            seq = ["^" + end[0]]
            adpt = adapter_parser.parse_multi([], [], seq)
            mod = AdapterCutter(adpt)
            mod.keep_match_info = True
            const_var[5] = mod
         if end[2]:
            seq = [end[2] + "$"]
            adpt = adapter_parser.parse_multi(seq, [], [])
            mod = AdapterCutter(adpt)
            mod.keep_match_info = True
            const_var[3] = mod

         # Adjust max error rate for adapter parser based on the ratio known
         # (constant) bases versus unknown (barcode) bases. In this way, N matches
         # "do not count" towards the max error rate
         tot_len = sum(len(s) for s in end)
         mer_new = max_error_rate * (tot_len - len(end[1])) / tot_len
         adapter_parser.constructor_args["max_error_rate"] = mer_new

         max_errs = int(tot_len * mer_new)
         logging.debug("Changed max_error_rate for %s to %f (max %d error(s))" % (
            _ID, mer_new, max_errs))

         # Update the "read_end"
         read_end = ["".join(end)]

      if _type == "index":
         sequences, index_names, max_error_rate_new = process_index(read_end,
            max_error_rate, override = True)

         # Change max_error_rate if a new one is given
         adapter_parser.constructor_args["max_error_rate"] = max_error_rate_new

         read_end = sequences

      # Fix the read position if requested using cutadapt special characters
      if _pos in ("fixed", "-"):
         if which == "_5":
            read_end = ["^" + x for x in read_end]
         else:
            read_end = [x + "$" for x in read_end]

   # Return the "processed" element and an optional linked parser
   return read_end, const_var, index_names

def create_adapters(features, adapter_parser, max_error_rate=0.1):
   """
   Create Adapter classes for each element in the read matrix, and return these
   adapters as lists. For information about the input, see "read_structure.txt"
   """

   adapters = []

   # try:
   for _ID, _5, _3, _type, _req, _second, _pos, keep_bases in features:

      # Assert that everything is looking okay
      assert _ID not in ("-", ""), "ID cannot be empty"
      assert _5 != "-" or _3 != "-", "Both sequences cannot be empty"
      assert _type in ("const", "barcode", "index", "const_bar", "DNA", "const_DNA",
                       "const_bar_comp"), "unknown type-argument"
      assert _req in ("-", "present", "absent") or _req.startswith("group:"), "unknown req-argument"
      assert _second in ("True", "False", "-"), "unknown second-argument"
      assert _pos in ("-", "fixed", "var"), "unknown pos-argument"

      # For simplicity:
      if _type in ("index", "barcode", "const_bar", "DNA", "const_DNA", "const_bar_comp"):
         assert _5 == "-" or _3 == "-", "Both ends doesn't make sense for this type"

      # Index must be present - what's the output otherwise?
      if _type == "index":
         _req = True
      elif _req == "absent":
         _req = False
      elif _req == "present":
         _req = True

      if _second == "True":
         _second = True
      else:
         _second = False

      if keep_bases != "-":
         try:
            keep_bases = int(keep_bases)
         except ValueError:
            raise "keep-bases should be an integer"
         assert _5 == "-" or _3 == "-", "Both ends doesn't make sense when keeping bases"
         if _5:
            keep_bases = (keep_bases, "_5")
         else:
            keep_bases = (keep_bases, "_3")
      else:
         keep_bases = None

      # Process sequences for 5' and 3' adapters
      # "const_var" is an (optional) second linked adapter, or None
      _5, const_var, index_names1 = process_end(_5, "_5", _ID, _type, _pos,
                                                max_error_rate, adapter_parser)
      _3, const_var, index_names2 = process_end(_3, "_3", _ID, _type, _pos,
                                                max_error_rate, adapter_parser, const_var)
      index_names = index_names1 + index_names2

      # Create adapter class, stored in a tuple with other important information
      adpt = (adapter_parser.parse_multi(_3, [], _5), _ID, _type, _req,
              _second, keep_bases, const_var, index_names)

      adapters.append(adpt)

      # Return max_error_rate if swapped for class "const_bar"
      adapter_parser.constructor_args["max_error_rate"] = max_error_rate

#   except ValueError as e:
#       raise BaseError("Error in read parsing")

   logging.debug("%d adapters listed\n" % len(adapters))

   return adapters

def create_output_name(read_file, output_dir, basename):
   """
   - Create output_base from given output directory and the basename
   - Create the directory if not present and report this
   """

   # Set basename for output files
   if not basename:
      basename = os.path.basename(read_file).split(".")[0]
   output_base = os.path.join(output_dir, basename)
   logging.debug("Output base: %s" % output_base)

   return output_base

def best_match_index(self, read):
   """
   Cutadapt adapter function:
      Find the best matching adapter in the given read.
      Return either a Match instance or None if there are no matches.

   Modified to store the adapter matched to ("self.adapter")
   """
   best = None
   self.adapter = None
   for adapter in self.adapters:
      match = adapter.match_to(read)
      if match is None:
         continue

      # the no. of matches determines which adapter fits best
      if best is None or match.matches > best.matches:
         best = match
         try:
            self.adapter = adapter.sequence
         except AttributeError:
            self.adapter = None
   return best

def find_var(read, const_var, give = "dna"):
   """
   Small function to extract the barcode sequence from the found match, using
   a very lenient linked adapter

   give can be "dna" (just the DNA stretch) or "sequence" (a sequence object)
   """
   if read.match_info:
      # Extract sequence and create new "Sequence" object
      seq = read.match_info[0][5]
      qual = read.match_info[0][9]
      match = seqio.Sequence(read.name, seq, qual)

      # Find barcode
      for i in const_var:
         if const_var[i]:
            match = const_var[i](match)
            if seq == match.sequence:
               #logging.debug("Pattern matched but no barcode found for read %s - %s" % (
               #              match.name, match.sequence))
               return ""
            seq = match.sequence

      if give == "dna":
         barcode = match.sequence
         return barcode
      else:
         return match

   return ""

def barcode_in_readname(_name, _ID, barcode):
   """
   Add the barcode found to the read name, before the first space.

   Example:
      read.name = barcode_in_readname(read.name, _ID, barcodes[i])
   """

   _name = _name.split(" ")
   _name[0] += "|" + _ID + ":" + barcode
   _name = " ".join(_name)

   return _name

def gather_statistics(stats, modifiers, statistics_file):
   """
   Create a small text file with the trimmed adapter counts
   """
   stats, n_written, n_tooshort = stats
   n, bp_total, bp_written = stats.n, stats.total_bp, stats.total_written_bp

   logging.debug("\n%d reads processed, %d reads written" % (n, n_written))

   IDs = []
   counts = []

   for modifier, _ID, _type, _req, _second, keep, const_bar in modifiers:
      IDs.append(_ID)
      counts.append(modifier.with_adapters)

   print("\t".join(["reads", "bp", "reads_written", "bp_written", "n_tooshort"] + IDs),
         file=statistics_file)
   print("\t".join([str(x) for x in [n, bp_total, n_written, bp_written, n_tooshort] + counts]),
         file=statistics_file)

   # Add versions
   print("## Versions ##", file=statistics_file)
   print("# read_parser: %s" % version, file=statistics_file)
   print("# cutadapt: %s" % __version__, file=statistics_file)

###############################################################################
##### Main ####################################################################
###############################################################################

version = "1.6"

def main():

   optParser = optparse.OptionParser(

      usage = "%prog [options] read-file read-structure output-dir",

      description=
         "This script / module takes a read file (or paired end) and parses " +
         "this based on a read-structure file given, writing to the output " +
         "directory + basename.",

      epilog =
         "Written as part of the Module-system - " +
         "160717 - Tom van Schaik" )

   optParser.add_option("-p", "--paired", type = "string", dest = "paired",
      default = None,
      help = "Paired end read-file")

   optParser.add_option("-b", "--basename", type = "string", dest = "basename",
      default = None,
      help = "basename of <output files> (default: input basename)")

   optParser.add_option("-q", "--quiet", action = "store_true", dest = "quiet",
      help = "suppress progress report")

   optParser.add_option("-a", "--all", action = "store_true", dest = "all_file",
      help = "report complete read structure table")

   optParser.add_option("-n", "--number", type = "int", dest = "all_number",
      default = 1000,
      help = "if all-file is given, the number of reads to be put here " +
             "(Note: more than number of reads present = all) (default: 1000)")

   optParser.add_option("-s", "--seq", action = "store_true", dest = "barcode_seq",
      help = "Include read sequence(s) in barcode file")

   optParser.add_option("-o", "--stdout", action = "store_true", dest = "stdout",
      help = "Write filtered reads to stdout (only for single end and no index!)")

   optParser.add_option("-r", "--bc-named", action = "store_true", dest = "bc_named",
      help = "Add barcode sequence to read name in output")

   optParser.add_option("-e", "--max_error_rate", type = "float", dest = "max_error_rate",
      default = 0.1, help = "Maximum error rate (default: 0.1)")

   optParser.add_option("-m", "--min_overlap", type = "int", dest = "min_overlap",
      default = 3, help = "Minimum overlap (default: 3)")

   optParser.add_option("-M", "--min_length", type = "int", dest = "min_length",
      default = None, help = "Minimum length (default: None)")

   optParser.add_option("-l", "--log_file", type = "string", dest = "log_file",
      default = False, help = "Output debug messages to log file instead of stderr")

   optParser.add_option("-x", "--read_wildcards", action = "store_true",
      default = False, dest = "read_wildcards",
      help = "Match read wildcards (default: False)")

   if len( sys.argv ) == 1:
      optParser.print_help()
      sys.exit(1)

   (opts, args) = optParser.parse_args()

   if len( args ) != 3:
      sys.stderr.write( sys.argv[0] + ": Error: Please provide three arguments.\n" )
      sys.stderr.write( "  Call with '-h' to get usage information.\n" )
      sys.exit( 1 )

   try:
      os.makedirs(args[2])
   except OSError:
      pass

   # Setup logger
   setup_logging(opts.quiet, opts.log_file)

   logging.debug("---- read_parser: version %s -----" % version)

   # Set-up output basename
   output_base = create_output_name(args[0], args[2], opts.basename)

   # Main
   reads = [args[0], opts.paired]
   read_matrix = args[1]

   run(reads, read_matrix, output_base,
       all_file = opts.all_file,
       all_number = opts.all_number,
       barcode_seq = opts.barcode_seq,
       stdout = opts.stdout,
       bc_named = opts.bc_named,
       max_error_rate = opts.max_error_rate,
       min_overlap = opts.min_overlap,
       min_length = opts.min_length,
       match_read_wildcards = opts.read_wildcards)

if __name__ == "__main__":
   main()
