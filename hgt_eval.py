#!/usr/bin/env python

import pysam
import numpy as np
import sys
import csv
import argparse
import operator
import random
import time
import string

__author__ = "Kathrin Trappe"



# parse candidate file and keep entries with acc-don or don-acc translocations
def get_candidates(candfile, acceptor, donor):
    with open(candfile) as input:
        for line in input:
            line = line.replace('\n', '')
            if (line.split('\t')[0] == acceptor and line.split('\t')[3] == donor):
                a_pos = int(line.split('\t')[1])
                a_base = line.split('\t')[2]
                d_pos = int(line.split('\t')[4])
                d_base = line.split('\t')[5]
                split_support = int(line.split('\t')[6])
                yield(CAND(a_pos, a_base, d_pos, d_base, split_support))
            if (line.split('\t')[3] == acceptor and line.split('\t')[0] == donor):
                a_pos = int(line.split('\t')[4])
                a_base = line.split('\t')[5]
                d_pos = int(line.split('\t')[1])
                d_base = line.split('\t')[2]
                split_support = int(line.split('\t')[6])
                yield(CAND(a_pos, a_base, d_pos, d_base, split_support))

def read_primary_pairs(inputfile):
    first = None
    second = None
    name = None
    for aln in inputfile:
        if (aln.qname != name) and (name is not None):
            if (first is not None) and (second is not None):
                yield (first, second)
            first = None
            second = None
        name = aln.qname
        if aln.is_secondary:
            continue
        if aln.is_read1:
            first = aln
        if aln.is_read2:
            second = aln
    if (first is not None) and (second is not None):
        yield (first, second)

def read_pairs(inputfile):
    first = None
    second = None
    name = None
    for aln in inputfile:
        if (first is not None) and (second is not None):
            yield (first, second)
        if (aln.qname is not name) and (name is not None):
            first = None
            second = None
        name = aln.qname
        #if aln.is_secondary:
            #continue
        if aln.is_read1:
            first = aln
        if aln.is_read2:
            second = aln
    if (first is not None) and (second is not None):
        yield (first, second)

# Get all read names from phage sam file
def read_phage_pairs(inputfile):
    name = None
    for aln in inputfile:
        if (aln.qname != name):
            name = aln.qname
            yield (name)

def get_random_regions(hgt_size, tabu_start, tabu_end, gen_size, num_regions):
    randrange = random.randrange
    return (get_random_region(hgt_size, tabu_start, tabu_end, gen_size, randrange) for x in range(num_regions))
    
# Generate random regions outside of candidate HGT region ("tabu region")
def get_random_region(hgt_size, tabu_start, tabu_end, gen_size, randrange):
    start = 0
    end = gen_size + 1
    while (end > gen_size):
        start = randrange(0, gen_size - hgt_size)
        end = start + hgt_size
        # random region is overlapping with tabu region
        if (start == tabu_start) or (end == tabu_end) or (start < tabu_start and end > tabu_start) or (start < tabu_end and end > tabu_end):
            end = gen_size + 1
    return (start, end)

# Add donor pair if matching, i.e. both mates lie within donor HGT region
def adpim(don_start, don_end, aln1, aln2):
    if (aln1.aend < don_start) or (aln1.pos > don_end): return 0
    if (aln2.aend < don_start) or (aln2.pos > don_end): return 0
    return 1

# Add pair if matching, i.e. first mate is on acceptor close to HGT boundary, second within donor HGT region
def apim(acc_start, acc_end, don_start, don_end, don_size, aln1, aln2):
    if (aln1.aend < acc_end) and (aln1.pos > acc_start): return 0
    if (aln1.aend < (acc_start - don_size/2)) or (aln1.pos > (acc_end + don_size/2)): return 0
    if (aln2.aend < don_start) or (aln2.pos > don_end): return 0
    return 1

# Add phage hit if read pair is in phage_list
def phage_hit(qname, phage_list):
    if qname in phage_list:
        return 1
    return 0

# CAND class holding information of single boundaries
class CAND:
    def __init__(self, acc_pos, acc_base, don_pos, don_base, split_support):
        self.acc_pos = acc_pos
        self.acc_base = acc_base
        self.don_pos = don_pos
        self.don_base = don_base
        self.split_support = split_support

# HGT class holding all HGT related infos including sampling values
# acc_length/don_length: length of acceptor/donor genome (required for sampling)
# add_pair_if_matching functions evaluate two given alignments (reads) if they fall into the HGT boundaries
class HGT:
    def __init__(self, acc_start, acc_end, acc_start_base, acc_end_base, don_start, don_end, split_support, options, acc_covs, don_covs):
        self.acc_start = acc_start
        self.acc_end = acc_end
        self.acc_start_base = acc_start_base
        self.acc_end_base = acc_end_base
        self.don_start = don_start
        self.don_end = don_end
        self.don_size = don_end - don_start
        self.acc_mean = np.mean(acc_covs[acc_start:acc_end],)
        self.don_mean = np.mean(don_covs[don_start:don_end],)
        self.don_pair_support = 0
        self.pair_support = 0
        self.phage_hits = 0
        self.split_support = split_support
        # Create sampling regions
        self.settup_bootstrap(options, acc_covs, don_covs)
    def add_pair_if_matching(self, aln1, aln2, phage_list):
        self.add_pair_if_matching_bootstrap(aln1, aln2)
        if (apim(self.acc_start, self.acc_end, self.don_start, self.don_end, self.don_size, aln1, aln2) != 0):
            self.pair_support += 1
            self.phage_hits += phage_hit(aln2.qname, phage_list)
    def add_don_pair_if_matching(self, aln1, aln2, phage_list):
        self.add_don_pair_if_matching_bootstrap(aln1, aln2)
        if (aln1.aend < self.don_start) or (aln1.pos > self.don_end): return
        if (aln2.aend < self.don_start) or (aln2.pos > self.don_end): return
        self.don_pair_support += 1
        self.phage_hits += phage_hit(aln1.qname, phage_list)
    def settup_bootstrap(self, options, acc_covs, don_covs):
        # Bootstrapping regions
        it = list(get_random_regions((self.acc_end - self.acc_start), self.acc_start, self.acc_end, options.acc_length, options.boot_num))
        l1, l2 = zip(*it)
        self.boot_acc_start_list = list(l1)
        self.boot_acc_end_list = list(l2)
        self.boot_acc_coverage_list = list(map(lambda x: np.mean(acc_covs[x[0]:x[1]]), it))
        it = list(get_random_regions((self.don_end - self.don_start), self.don_start, self.don_end, options.don_length, options.boot_num))
        l1, l2 = zip(*it)
        self.boot_don_start_list = list(l1)
        self.boot_don_end_list = list(l2)
        self.boot_don_coverage_list = list(map(lambda x: np.mean(don_covs[x[0]:x[1]]), it))
        if (options.paired_reads):
            self.boot_crossing_pairs_list = [0] * options.boot_num
            self.boot_don_pairs_list = [0] * options.boot_num
    def add_pair_if_matching_bootstrap(self, aln1, aln2):
        hits = map(lambda x: apim(x[0], x[1], self.don_start, self.don_end, self.don_size, aln1, aln2), zip(self.boot_acc_start_list, self.boot_acc_end_list))
        self.boot_crossing_pairs_list = [sum(x) for x in zip(self.boot_crossing_pairs_list, hits)]
    def add_don_pair_if_matching_bootstrap(self, aln1, aln2):
        hits = map(lambda x: adpim(x[0], x[1], aln1, aln2), zip(self.boot_don_start_list, self.boot_don_end_list))
        self.boot_don_pairs_list = [sum(x) for x in zip(self.boot_don_pairs_list, hits)]

# Check for duplicate HGT entry, i.e. a HGT with all 4 boundaries within tolerance, resp.
def duplicate_entry(hgt_list, tolerance, acc_start, acc_end, don_start, don_end, split_support):
    for i in range(0, len(list(hgt_list))):
        last_hgt = hgt_list[i]
        if (abs(last_hgt.acc_start - acc_start) < tolerance and
            abs(last_hgt.acc_end - acc_end) < tolerance and
            abs(last_hgt.don_start - don_start) < tolerance and
            abs(last_hgt.don_end - don_end) < tolerance):
            last_hgt.split_support += split_support
            return 1
    return 0

def main():
    # Parsing options
    parser = argparse.ArgumentParser(description='Hgt candidate evaluation.')
    parser.add_argument('bamfile', help="BAM alignments file for both acceptor and donor, sorted via qname")
    parser.add_argument('candfile', help="File with inter-chr translocation breakpoints")
    parser.add_argument('acceptor', help="Name of acceptor reference (gi as in fasta)")
    parser.add_argument('donor', help="Name of donor reference (gi as in fasta)")

    # optional arguments (has default values)
    parser.add_argument('--phagefile', default=None, help="SAM alignments file for phage database")
    parser.add_argument('-o','--outfile', default="hgt_eval.vcf", help="Output VCF file with filtered, evaluated candidates")
    parser.add_argument('-min', dest='min_size', default=100, type=int, help="minimal HGT size, default 100")
    parser.add_argument('-max', dest='max_size', default=50000, type=int, help="maximal HGT size, default 50000")
    parser.add_argument('--tolerance', default=20, type=int, help="Position range to remove duplicates, default 20")
    parser.add_argument('--pair-support', dest='paired_reads', default=True, action='store_false', help="Turn on/off paired reads support, default TRUE")
    parser.add_argument('--num-boot-regions', dest='boot_num', default=100, type=int, help="Number of sampled regions in acceptor for bootstrapping expected number of pairs and coverage for filtering, default 100")
    parser.add_argument('--boot-sens', dest='boot_sens', default=95, type=int, choices=range(0, 100), help="Percent of cases for the candidate region to exceed values of sampled regions, default 95 percent")

    options = parser.parse_args()

    # Extracting parameters
    hgt_minLength = options.min_size
    hgt_maxLength = options.max_size
    bf = pysam.Samfile(options.bamfile,'rb')
    if (options.phagefile is not None):
        psf = pysam.Samfile(options.phagefile,'r')
    acc_tid = bf.gettid(options.acceptor)
    don_tid = bf.gettid(options.donor)
    options.out_all = options.outfile.strip('vcf')+'tsv'

    print ('acc_tid', acc_tid, options.acceptor)
    print ('don_tid', don_tid, options.donor)
    print ('hgt_minLength', hgt_minLength)
    print ('hgt_maxLength', hgt_maxLength)
    print ('paired_reads', options.paired_reads)
    print ('num_boot_regions', options.boot_num)
    print ('boot_sensitivity', options.boot_sens)

    # Getting acceptor genome length
    options.acc_length = bf.lengths[acc_tid]
    options.don_length = bf.lengths[don_tid]

    # Extracting coverages
    covs = [np.zeros((l,)) for l in bf.lengths]
    num_matches = 0 # number of read matches, unused by now
    for read in bf:
        if not read.is_unmapped:
            r_start = read.pos # start position
            r_end = read.pos + read.qlen # end
            covs[read.tid][r_start:r_end] += 1
            num_matches += 1
    bf.reset()

    acc_cov = covs[acc_tid]
    don_cov = covs[don_tid]

    # Extracting phage read IDs, if any
    phage_list = []
    if (options.phagefile != None):
        phage_list = [qname for qname in read_phage_pairs(psf)]

    # Extracting candidate positions from candifate file
    cand_list = [cand for cand in get_candidates(options.candfile, options.acceptor, options.donor)]
    
    #indices = range(len(acc_pos))
    #indices.sort(key=acc_pos.__getitem__)
    cand_list.sort(key=operator.attrgetter('acc_pos'))

    tstart = time.clock()

    # Pair up single boundary candidates to HGT candidates conforming to size constraints
    hgt_list = []
    for ps in range(0, len(cand_list)):
        cand_start = cand_list[ps]
        acc_start = cand_start.acc_pos
        don_start_temp = cand_start.don_pos
        for pe in range(ps, len(cand_list)):
            cand_end = cand_list[pe]
            acc_end = cand_end.acc_pos + 1
            don_end_temp = cand_end.don_pos + 1
            # Exchange don_start/end so that don_start < don_end (we don't care about the orientation for now)
            if (don_end_temp < don_start_temp):
                don_start = don_end_temp
                don_end = don_start_temp
            else:
                don_start = don_start_temp
                don_end = don_end_temp
            # Avoid similar hgt entries within tolerance
            if (duplicate_entry(hgt_list, options.tolerance, acc_start, acc_end, don_start, don_end, cand_end.split_support)):
                continue
            # Continue if acceptor HGT region exceeds max length
            if (abs(acc_end - acc_start) > hgt_maxLength):
                break
            if (abs(don_end - don_start) > hgt_minLength and abs(don_end - don_start) < hgt_maxLength):
                split_support = (cand_start.split_support + cand_end.split_support)
                hgt_list.append(HGT(acc_start, acc_end, cand_start.acc_base, cand_end.acc_base, don_start, don_end, split_support, options, acc_cov, don_cov))

    print (time.clock() - tstart)

    # Get (all/primary) read pairs which map to both donor and acceptor and add them to hgt attributes
    if (options.paired_reads): 
        #for aln1, aln2 in read_pairs(bf):
        for aln1, aln2 in read_primary_pairs(bf):
            if aln1.is_unmapped:
                continue
            if aln2.is_unmapped:
                continue
            if (aln1.tid == don_tid) and (aln2.tid == don_tid):
                for hgt in hgt_list:
                    hgt.add_don_pair_if_matching(aln1, aln2, phage_list)
            if (aln1.tid != acc_tid) or (aln2.tid != don_tid):
                aln1, aln2 = aln2, aln1 # Ensures that aln1 belongs to acceptor and aln2 to donor
            if (aln1.tid != acc_tid) or (aln2.tid != don_tid):
                continue
            for hgt in hgt_list:
                hgt.add_pair_if_matching(aln1, aln2, phage_list)
    print ("total sample time ", time.clock() - tstart)
    
    # Output results       
    print ("writing output")
    with open(options.outfile, 'w') as vcf_out, open(options.out_all, 'w') as output:
        # VCF header
        writevcf = csv.writer(vcf_out, delimiter = '\t')
        writevcf.writerow(['##fileformat=VCFv4.2'])
        writevcf.writerow(['##source=DAISY'])
        writevcf.writerow(['##INFO=<ID=EVENT,Number=1,Type=String,Description=\"Event identifier for breakends.\">'])
        writevcf.writerow(['##contig=<ID='+options.acceptor+'>'])
        writevcf.writerow(['##contig=<ID='+options.donor+'>'])
        writevcf.writerow(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

        # TSV header
        writetsv = csv.writer(output, delimiter = '\t')
        writetsv.writerow(['#AS: Acceptor start position'])
        writetsv.writerow(['#AE: Acceptor end position'])
        writetsv.writerow(['#DS: Donor start position'])
        writetsv.writerow(['#DE: Donor end position'])
        writetsv.writerow(['#MC: Mean coverage in region'])
        writetsv.writerow(['#Split: Total number split-reads per region (including duplicates!)'])
        writetsv.writerow(['#PS-S: Pairs spanning HGT boundaries'])
        writetsv.writerow(['#PS-W: Pairs within HGT boundaries'])
        writetsv.writerow(['#Phage: PS-S and PS-W reads mapping to phage database'])
        writetsv.writerow(['#BS:MC/PS-S/PS-W: Percent of bootstrapped random regions with MC/PS-S/PS-W smaller than candidate'])
        if (options.paired_reads):
            writetsv.writerow(['AS', 'AE', 'MC', 'BS:MC', 'DS', 'DE', 'MC', 'Split', 'PS-S', 'PS-W', 'Phage', 'BS:MC', 'BS:PS-S', 'BS:PS-W'])
        else:
            writetsv.writerow(['AS', 'AE', 'MC', 'BS:MC', 'DS', 'DE', 'MC', 'Split', 'BS:MC'])

        # write candidates
        hgt_vcf_counter = 1
        # Define sensitivity values for candidates to be reported in VCF
        sens = float(options.boot_sens)/float(100)
        sens_acc = 1.0 - sens
        thresh = sens * float(options.boot_num)
        for hgt in hgt_list:
            # evaluate bootstrap
            acc_cov_test = sum(1 for i in range(len(list(hgt.boot_acc_coverage_list))) if hgt.acc_mean > hgt.boot_acc_coverage_list[i])
            don_cov_test = sum(1 for i in range(len(list(hgt.boot_don_coverage_list))) if hgt.don_mean > hgt.boot_don_coverage_list[i])
            if (options.paired_reads):
                # write all candidates to tsv file
                cross_pair_test = sum(1 for i in range(len(list(hgt.boot_crossing_pairs_list))) if hgt.pair_support > hgt.boot_crossing_pairs_list[i])
                don_pair_test = sum(1 for i in range(len(list(hgt.boot_don_pairs_list))) if hgt.don_pair_support > hgt.boot_don_pairs_list[i])
                phage_support = 0
                if ((hgt.pair_support + hgt.don_pair_support) > 0):
                    phage_support = float(hgt.phage_hits)/float((hgt.pair_support + hgt.don_pair_support))
                writetsv.writerow([hgt.acc_start, hgt.acc_end, "%.2f" % hgt.acc_mean, acc_cov_test, hgt.don_start, hgt.don_end, "%.2f" % hgt.don_mean, hgt.split_support,  hgt.pair_support, hgt.don_pair_support, "%.4f" % phage_support, don_cov_test, cross_pair_test, don_pair_test])
            else:
                writetsv.writerow([hgt.acc_start, hgt.acc_end, "%.2f" % hgt.acc_mean, acc_cov_test, hgt.don_start, hgt.don_end, "%.2f" % hgt.don_mean, hgt.split_support, don_cov_test])

            # Write only canidates to VCF file that passed the filter (boot_sens)
            if (options.boot_sens > 0):
                if (acc_cov_test < thresh) and (acc_cov_test > sens_acc * options.boot_num): continue
                if (options.paired_reads):
                    if (cross_pair_test < thresh) or (don_pair_test < thresh): continue
                if (don_cov_test < thresh): continue

            # write filtered candidates to VCF file
            writevcf.writerow([options.acceptor, hgt.acc_start, 'BND_'+str(hgt_vcf_counter)+'_1', hgt.acc_start_base, hgt.acc_start_base+'['+options.donor+':'+str(hgt.don_start)+'[', 'PASS', 'SVTYPE=BND;EVENT=HGT'+str(hgt_vcf_counter), '.', '1'])
            writevcf.writerow([options.acceptor, hgt.acc_end, 'BND_'+str(hgt_vcf_counter)+'_2', hgt.acc_end_base, ']'+options.donor+':'+str(hgt.don_end)+']'+hgt.acc_end_base, 'PASS', 'SVTYPE=BND;EVENT=HGT'+str(hgt_vcf_counter), '.', '1'])
            hgt_vcf_counter += 1

        if (hgt_vcf_counter == 1):
            print ('No canidates written to VCF, try to rerun with lower sampling sensitivity (--boot_sens)')


if __name__ == '__main__':
        #profile.run('re.compile("main")')
        main()

