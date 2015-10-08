#! /usr/bin/env python

import sys
import os
import argparse
import subprocess32 as sp
import select
import logging
import time
import hgt_eval

# https://gist.github.com/hangtwenty/6390750
def call(popenargs, logger, stdout_log_level=logging.DEBUG, stderr_log_level=logging.ERROR, **kwargs):
    """
    Variant of subprocess.call that accepts a logger instead of stdout/stderr,
    and logs stdout messages via logger.debug and stderr messages via
    logger.error.
    """

    try:
        child = sp.Popen(popenargs, stdout=sp.PIPE, stderr=sp.PIPE, **kwargs)

        log_level = {child.stdout: stdout_log_level, child.stderr: stderr_log_level}

        def check_io():
            ready_to_read = select.select([child.stdout, child.stderr], [], [], 1000)[0]
            for io in ready_to_read:
                line = io.readline()
                if line:
                    logger.log(log_level[io], line[:-1])

        # keep checking stdout/stderr until the child exits
        while child.poll() is None:
            check_io()
        check_io()  # check again to catch anything after the process exits 
    except Exception as e:
        raise(e)
    
    return child.wait()

def printAndWrite(pmessage, lmessage, logger, mode):
    print(pmessage)
    if mode=='warning':
        logger.warning(lmessage)
    elif mode=='info':
        logger.info(lmessage)
    elif mode=='exception':
        logger.exception(lmessage)
    else:
        logger.error(lmessage)

# check if all/one file exists
def checkExistence(logger, job, *filename):
    b = True
    for f in filename:
        b = b and os.path.isfile(f)
    if b:
        logger.info('Skipping {}. File(s): {} already exist(s).'.format(job, ' '.join(filename)))
    return b

def validInt(s):
    try:
        int(s)
        return(str(s))
    except:
        raise argparse.ArgumentTypeError('{} is not an int.'.format(s))

# parse donor file with GIs
def get_GIs(donorfile):
    with open(donorfile) as input:
        for line in input:
            line = line.replace('\n', '')
            name = line.split('|')[3]
            yield(line,name)

def producedict():
    argsD = dict()
    argsD['acceptor'] = '--acceptor'
    argsD['accref'] = '--accref'
    argsD['b_bed'] = '--bedTools'
    argsD['b_cat'] = '--concatRef'
    argsD['b_eval'] = '--hgt_eval'
    argsD['b_gmj'] = '--gustafMateJoin'
    argsD['b_can'] = '--get_candidates'
    argsD['b_gustaf'] = '--gustaf'
    argsD['b_indexSam'] = '--indexSam'
    argsD['b_new'] = '--overwrite'
    argsD['b_phage'] = '--phage_mapping'
    argsD['b_sak'] = '--sak'
    argsD['b_sortBamQ'] = '--sortBam'
    argsD['b_sortMapping'] = '--sortMapped'
    argsD['b_stellar'] = '--stellar'
    argsD['b_yaraI'] = '--yaraIndex'
    argsD['b_yaraM'] = '--yaraMapper'
    argsD['donor'] = '--donor'
    argsD['donor2'] = '--donor2'
    argsD['donref'] = '--donref'
    argsD['donref2'] = '--donref2'
    argsD['g_bth'] = '--gustaf_bth'
    argsD['g_cbp'] = '--gustaf_cbp'
    argsD['g_gth'] = '--gustaf_gth'
    argsD['g_ith'] = '--gustaf_ith'
    argsD['g_le'] = '--gustaf_le'
    argsD['g_ll'] = '--gustaf_ll'
    argsD['g_nth'] = '--gustaf_nth'
    argsD['g_st'] = '--gustaf_st'
    argsD['h_min'] = '--hgteval_minsize'
    argsD['h_max'] = '--hgteval_maxsize'
    argsD['h_tol'] = '--hgteval_tolerance'
    argsD['h_pairs'] = '--hgteval_pairsupport'
    argsD['h_bootnum'] = '--hgteval_bootnum'
    argsD['h_bootsens'] = '--hgteval_bootsens'
    argsD['outdir'] = '--outdir'
    argsD['phage_dir'] = '--phage_dir'
    argsD['phage_ref'] = '--phage_ref'
    argsD['read1fasta'] = '--read1'
    argsD['read2fasta'] = '--read2'
    argsD['rev'] = '--reverse'
    argsD['sam_F'] = '--sam_F'
    argsD['stellar_l'] = '--stellar_l'
    argsD['task'] = '--task'
    argsD['yara_e'] = '--yara_e'
    argsD['yara_le'] = '--yara_le'
    argsD['yara_ll'] = '--yara_ll'
    argsD['yara_t'] = '--yara_t'
    return argsD

def parser():
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',description='DAISY',add_help=False)
    required = parser.add_argument_group('required arguments')
    mutrequired = parser.add_mutually_exclusive_group(required=True)
    optional = parser.add_argument_group('optional arguments')
    yara = parser.add_argument_group('Yara arguments')
    sam = parser.add_argument_group('Samtools arguments')
    stellar  = parser.add_argument_group('Stellar arguments')
    gustaf = parser.add_argument_group('Gustaf arguments')
    hgt = parser.add_argument_group('HGT evaluation arguments')
    optional.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    optional.add_argument('-od', '--outdir', dest='outdir', nargs='?', type=str, default='./', help='Output directory for all output files. (Default: %(default)s)')
    optional.add_argument('-pod', '--phage_dir', dest='phage_dir', nargs='?', type=str, default='./', help='Output/input directory phage index. (Default: %(default)s)')
    required.add_argument('-r1', '--read1', dest='read1fasta', nargs='?', type=str, required=True, default=None, help='Single-end reads/Paired-end reads 1')
    optional.add_argument('-r2', '--read2', dest='read2fasta', nargs='?', type=str, default=None, help='Paired-end reads 2')
    required.add_argument('-ar', '--acceptor_ref', dest='accref', nargs='?', type=str, required=True, default=None, help='Reference file of acceptor')
    required.add_argument('-dr', '--donor_ref', dest='donref', nargs='?', type=str, required=True, default=None, help='Reference file of donor')
    optional.add_argument('-dr2', '--donor_ref2', dest='donref2', nargs='?', type=str, default=None, help='Second reference file of donor')
    optional.add_argument('-pr', '--phage_ref', dest='phage_ref', nargs='?', type=str, default=None, help='Phage database reference file')
    required.add_argument('-a', '--acceptor', dest='acceptor', nargs='?', type=str, required=True, default=None, help='Name of acceptor (gi from reference file)')
    mutrequired.add_argument('-d', '--donor', dest='donor', nargs='?', type=str, default=None, help='Name(s) of donor (gi from reference file)')
    mutrequired.add_argument('-d2', '--donor2', dest='donor2', nargs='?', type=str, default=None, help='File of names of donor (gi from reference file); one gi per line')
    optional.add_argument('-t', '--task', dest='task', nargs='?', type=str, default='', help='Define optional taskname. Will be used as prefix for created files.')
    optional.add_argument('-new', '--overwrite', dest='b_new', action='store_true', default=False, help='Overwrite existing files. (Default: %(default)s)')
    optional.add_argument('-rev', '--reverse', dest='rev', action='store_false', default=True, help='Define which modules are not to be run instead of which are to be run. (Default: %(default)s)')
    optional.add_argument('-f', '--saveArgs', dest='argF', nargs='?', type=argparse.FileType('w'), help='Store arguments to file <ARGF>. Use ./daisy.py @<ARGF> to use stored arguments.')
    optional.add_argument('-cr', '--concatRef', dest='b_cat', action='store_true', default=False, help='Run join candidate genomes.')
    optional.add_argument('-yi', '--yaraIndex', dest='b_yaraI', action='store_true', default=False, help='Run Yara index creation.')
    optional.add_argument('-ym', '--yaraMapper', dest='b_yaraM', action='store_true', default=False, help='Run Yara mapping.')
    optional.add_argument('-sm', '--sortMapped', dest='b_sortMapping', action='store_true', default=False, help='Run sort Mapped Sam.')
    optional.add_argument('-bed', '--bedTools', dest='b_bed', action='store_true', default=False, help='Run unmapped FASTQ creation.')
    optional.add_argument('-is', '--indexSam', dest='b_indexSam', action='store_true', default=False, help='Run index Sam.')
    optional.add_argument('-sb', '--sortBam', dest='b_sortBamQ', action='store_true', default=False, help='Run sort Bam by qname.')
    optional.add_argument('-sak', '--sak', dest='b_sak', action='store_true', default=False, help='Run sak.')
    optional.add_argument('-gmj', '--gustafMateJoin', dest='b_gmj', action='store_true', default=False, help='Run Gustaf mate joining.')
    optional.add_argument('-ste', '--stellar', dest='b_stellar', action='store_true', default=False, help='Run Stellar.')
    optional.add_argument('-gus', '--gustaf', dest='b_gustaf', action='store_true', default=False, help='Run Gustaf.')
    optional.add_argument('-can', '--get_candidates', dest='b_cand', action='store_true', default=False, help='Run extracting hgt candidates.')
    optional.add_argument('-pm', '--phage_mapping', dest='b_phage', action='store_true', default=False, help='Run Phage database mapping.')
    optional.add_argument('-eva', '--hgt_eval', dest='b_eval', action='store_true', default=False, help='Run HGT evaluation.')
    yara.add_argument('-yll', '--yara_ll', dest='yara_ll', nargs='?', type=validInt, default='200', help='Yara library length parameter (Default: %(default)s)')
    yara.add_argument('-yle', '--yara_le', dest='yara_le', nargs='?', type=validInt, default='20', help='Yara library error parameter (Default: %(default)s)')
    yara.add_argument('-ye', '--yara_e', dest='yara_e', nargs='?', type=validInt, default='10', help='Yara error rate parameter (Default: %(default)s)')
    yara.add_argument('-yt', '--yara_t', dest='yara_t', nargs='?', type=validInt, default='5', help='Yara thread parameter (Default: %(default)s)')
    sam.add_argument('-F', '--sam_F', dest='sam_F', nargs='?', type=validInt, default='2', help='F parameter for sorting unmapped Bam (Default: %(default)s)')
    stellar.add_argument('-l', '--stellar_l', dest='stellar_l', nargs='?', type=validInt, default='30', help='Stellar l parameter (Default: %(default)s)')
    gustaf.add_argument('-st', '--gustaf_st', dest='g_st', nargs='?', type=validInt, default='3', help='Gustaf st parameter (Default: %(default)s)')
    gustaf.add_argument('-bth', '--gustaf_bth', dest='g_bth', nargs='?', type=validInt, default='70', help='Gustaf bth parameter (Default: %(default)s)')
    gustaf.add_argument('-gth', '--gustaf_gth', dest='g_gth', nargs='?', type=validInt, default='100', help='Gustaf gth parameter (Default: %(default)s)')
    gustaf.add_argument('-ith', '--gustaf_ith', dest='g_ith', nargs='?', type=validInt, default='50', help='Gustaf ith parameter (Default: %(default)s)')
    gustaf.add_argument('-gll', '--gustaf_ll', dest='g_ll', nargs='?', type=validInt, default='500', help='Gustaf ll parameter (Default: %(default)s)')
    gustaf.add_argument('-gle', '--gustaf_le', dest='g_le', nargs='?', type=validInt, default='50', help='Gustaf le parameter (Default: %(default)s)')
    gustaf.add_argument('-cbp', '--gustaf_cbp', dest='g_cbp', action='store_true', default=False, help='Gustaf cbp parameter (Default: %(default)s)')
    gustaf.add_argument('-nth', '--gustaf_nth', dest='g_nth', nargs='?', type=validInt, default='5', help='Gustaf nth parameter (Default: %(default)s)')
    hgt.add_argument('-min', '--hgteval_minsize', dest='h_min', nargs='?', type=validInt, default='100', help='HGT evaluation min_size (Default: %(default)s)')
    hgt.add_argument('-max', '--hgteval_maxsize', dest='h_max', nargs='?', type=validInt, default='55000', help='HGT evaluation max_size (Default: %(default)s)')
    hgt.add_argument('-tol', '--hgteval_tolerance', dest='h_tol', nargs='?', type=validInt, default='20', help='HGT evaluation tolerance for HGT duplicate removal (Default: %(default)s)')
    hgt.add_argument('-num', '--hgteval_bootnum', dest='h_bootnum', nargs='?', type=validInt, default='100', help='HGT evaluation number of sample regions (Default: %(default)s)')
    hgt.add_argument('-sen', '--hgteval_bootsens', dest='h_bootsens', nargs='?', type=validInt, default='95', help='HGT evaluation percent of sample region sensitivity (Default: %(default)s)')
    return parser

def pipeline(args):

    tstart = time.time()

    rundir = os.path.dirname(sys.argv[0])+'/'

    if (args.read2fasta is not None):
        args.h_pairs = True
    else:
        args.h_pairs = False

    if args.rev:
        args.b_cat = not args.b_cat
        args.b_yaraI = not args.b_yaraI
        args.b_yaraM = not args.b_yaraM
        args.b_sortMapping = not args.b_sortMapping
        args.b_bed = not args.b_bed
        args.b_indexSam = not args.b_indexSam
        args.b_sortBamQ = not args.b_sortBamQ
        args.b_sak = not args.b_sak
        args.b_gmj = not args.b_gmj
        args.b_stellar = not args.b_stellar
        args.b_gustaf = not args.b_gustaf
        args.b_cand = not args.b_cand
        args.b_phage = not args.b_phage
        args.b_eval = not args.b_eval
    
    if args.argF:
        argsD = producedict()
        for arg in sorted(vars(args)):
            if arg != 'argF':
                args.argF.write('{}={}\n'.format(argsD[arg],getattr(args,arg)))
        args.argF.close()
    
    readname = os.path.basename(args.read1fasta).split('.')[0]
    accref = os.path.basename(args.accref).split('.')[0]
    donref = os.path.basename(args.donref).split('.')[0]
    if (args.donref2 is not None):
        donref2 = os.path.basename(args.donref2).split('.')[0]
    if (args.task != ''):
        args.task += '_'

    # Create a new logger if none exists.
    try:
        logger = args.logger
    except AttributeError:
        logging.basicConfig(format='%(levelname)s:%(message)s', filename='{}{}{}_{}_{}.log.txt'.format(args.outdir, args.task, readname, accref, donref), filemode='w',  level=logging.DEBUG)
        logger = logging.getLogger('log_subprocesses')

    # Joining candidate genomes
    # cat acc.fa don.fa > acc_don.fa

    candidateRef = '{}{}_{}.fa'.format(args.outdir, accref, donref)
    if (args.donref2 is not None):
        candidateRef = '{}{}_{}_{}.fa'.format(args.outdir, accref, donref, donref2)
    if args.b_cat and (args.b_new or not checkExistence(logger,'Joining Candidate Genomes',candidateRef)):
        cmd = 'cat {} {} > {} && sed -i "/^[^>]/ s/[^ACTG]/N/g" {}'.format(args.accref, args.donref, candidateRef, candidateRef)
        if (args.donref2 is not None):
            cmd = 'cat {} {} {} > {} && sed -i "/^[^>]/ s/[^ACTG]/N/g" {}'.format(args.accref, args.donref, args.donref2, candidateRef, candidateRef)
        printAndWrite('Start Joining Candidate Genomes', 'Start Joining Candidate Genomes', logger, 'info')
        logger.debug(cmd)
        try:
            call(cmd, logger, shell=True)
        except sp.CalledProcessError:
            printAndWrite('Error in cat. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        printAndWrite('End Joining Candidate Genomes', 'End Joining Candidate Genomes', logger, 'info')


    # yara indexer of joined candidate genomes
    #~/bin/yara_indexer acc_don.fa -o Acc_Don_index

    yaraIndex = '{}{}_{}'.format(args.outdir, accref, donref)
    if args.b_yaraI and (args.b_new or not checkExistence(logger, 'Yara Indexer', yaraIndex+'.lf.drp', yaraIndex+'.lf.drs', yaraIndex+'.lf.drv', yaraIndex+'.lf.pst', yaraIndex+'.rid.concat',
         yaraIndex+'.rid.limits', yaraIndex+'.sa.ind', yaraIndex+'.sa.len', yaraIndex+'.sa.val', yaraIndex+'.txt.concat', yaraIndex+'.txt.limits', yaraIndex+'.txt.size')):
        cmd = ['{}yara_indexer'.format(rundir), candidateRef, '-o', yaraIndex]
        printAndWrite('Start Yara Indexer', 'Start Yara Indexer', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Yara Indexer. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find yara_indexer in working directory. Trying global next.')
            cmd = ['yara_indexer', candidateRef, '-o', yaraIndex]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)
            except sp.CalledProcessError:
                printAndWrite('Error in Yara Indexer. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Yara Indexer', 'End Yara Indexer', logger, 'info')

    # yara mapper -all -e 10
    # ~/bin/yara_mapper -o *all.sam reads.1.fasta reads.2.fasta -a -e 10

    mapping = '{}{}{}_{}_{}.sam'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_yaraM and (args.b_new or not checkExistence(logger,'Yara Mapping',mapping)): 
        cmd = [x for x in ['{}yara_mapper'.format(rundir), yaraIndex, args.read1fasta, args.read2fasta, '-o', mapping, '-e', args.yara_e, '-ll', args.yara_ll, '-le', args.yara_le, '-t', args.yara_t] if x is not None]
        printAndWrite('Start Yara Mapping', 'Start Yara Mapping', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd , logger)
        except sp.CalledProcessError as err:
            printAndWrite('Error in Yara Mapper. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find yara_mapper in working directory. Trying global next.')
            cmd = [x for x in ['yara_mapper', yaraIndex, args.read1fasta, args.read2fasta, '-o', mapping, '-e', args.yara_e, '-ll', args.yara_ll, '-le', args.yara_le, '-t', args.yara_t] if x is not None]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)
            except sp.CalledProcessError:
                printAndWrite('Error in Yara Mapper. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Yara Mapper', 'End Yara Mapper', logger, 'info')


    # samtools --> sort unmapped bam
    # samtools view -b -F 2 -S *all.sam > *unmapped.bam
    # samtools sort -n *unmapped.bam *unmapped.sort

    unmappedBam = '{}{}{}_{}_{}_unmapped.bam'.format(args.outdir, args.task, readname, accref, donref)
    unmappedSortedBam = '{}{}{}_{}_{}_unmapped.sort'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_sortMapping and (args.b_new or not checkExistence(logger, 'Sorting Mapping', unmappedBam, unmappedSortedBam+'.bam')):
        cmd = ['{}samtools'.format(rundir), 'view', '-b', '-F', args.sam_F, '-S', mapping, '-o', unmappedBam]
        printAndWrite('Start Samtools Extract Unmapped to BAM', 'Start Samtools Extract Unmapped to BAM', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError: 
            printAndWrite('Error in Samtools Extract Unmapped to BAM. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find samtools in working directory. Trying global next.')
            cmd = ['samtools', 'view', '-b', '-F', args.sam_F, '-S', mapping, '-o', unmappedBam]
            logger.debug(' '.join(cmd))
            try:
                sp.call(cmd)           
            except sp.CalledProcessError:
                printAndWrite('Error in Samtools Extract Unmapped to BAM. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Samtools Extract Unmapped to BAM', 'End Samtools Extract Unmapped to BAM', logger, 'info')

        cmd = ['{}samtools'.format(rundir), 'sort', '-n', unmappedBam, unmappedSortedBam]
        printAndWrite('Start Samtools Sort Unmapped BAM', 'Start Samtools Sort Unmapped BAM', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Samtools Sort Unmapped BAM . Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError: 
            logger.warning('Could not find samtools in working directory. Trying global next.')
            cmd = ['samtools', 'sort', '-n', unmappedBam, unmappedSortedBam]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Samtools Sort Unmapped BAM . Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Samtools Sort Unmapped BAM', 'End Samtools Sort Unmapped BAM', logger, 'info')
    unmappedSortedBam += '.bam'

    # bedtools --> unmapped fastq
    # bedtools bamtofastq -i *unmapped.sort.bam -fq *unmapped.sort.1.fastq -fq2 *unmapped.sort.2.fastq

    unmappedSortedFastq1 = '{}{}{}_unmapped.sort.1.fastq'.format(args.outdir, args.task, readname)
    unmappedSortedFastq2 = None
    if (args.read2fasta is not None):
        unmappedSortedFastq2 = '{}{}{}_unmapped.sort.2.fastq'.format(args.outdir, args.task, readname)
    # Only checking existence of unmappedSortedFastq1 to avoid crashing for single-end (and for paired-end both have to be present anyways)
    if args.b_bed and (args.b_new or not checkExistence(logger, 'Converting to Unmapped FASTQ', unmappedSortedFastq1)):
        cmd = [x for x in ['{}bedtools'.format(rundir), 'bamtofastq', '-i', unmappedSortedBam, '-fq', unmappedSortedFastq1, '-fq2' if args.read2fasta is not None else None, unmappedSortedFastq2] if x is not None]
        printAndWrite('Start Converting to Unmapped FASTQ', 'Start Converting to Unmapped FASTQ', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Converting to Unmapped FASTQ. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find bedtools in working directory. Trying global next.')
            cmd = [x for x in ['bedtools', 'bamtofastq', '-i', unmappedSortedBam, '-fq', unmappedSortedFastq1, '-fq2' if args.read2fasta is not None else None, unmappedSortedFastq2] if x is not None]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Converting to Unmapped FASTQ. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Converting to Unmapped FASTQ', 'End Converting to Unmapped FASTQ', logger, 'info')


    # samtools view -bS *all.sam -o *all.bam

    mappedBam = '{}{}{}_{}_{}.bam'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_sortBamQ and checkExistence(logger, 'Generating', mapping) and (args.b_new or not checkExistence(logger, 'Generating Mapped BAM', mappedBam)):
        cmd = ['{}samtools'.format(rundir), 'view', '-b', '-S', mapping, '-o', mappedBam]
        printAndWrite('Start Converting Mapped Sam to BAM', 'Start Converting Mapped Sam to BAM', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Converting Mapped Sam to BAM. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find samtools in working directory. Trying global next.')
            cmd = ['samtools', 'view', '-b', '-S', mapping, '-o', mappedBam]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Converting Mapped Sam to BAM. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Converting Mapped Sam to BAM', 'End Converting Mapped Sam to BAM', logger, 'info')


    # samtools --> sort mapped bam (qname)
    # samtools sort -n *all.bam *all.sort-qname

    mappedSortedQBam = '{}{}{}_{}_{}.sort-qname'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_sortBamQ and (args.b_new or not checkExistence(logger, 'Sorting BAM', mappedSortedQBam+'.bam')):
        cmd = ['{}samtools'.format(rundir), 'sort', '-n', mappedBam, mappedSortedQBam]
        printAndWrite('Start Sorting BAM', 'Start Sorting BAM', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Sorting BAM (qname). Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find samtools in working directory. Trying global next.')
            cmd = ['samtools', 'sort', '-n', mappedBam, mappedSortedQBam]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Sorting BAM (qname). Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Sorting BAM (qname)', 'End Sorting BAM (qname)', logger, 'info')
    mappedSortedQBam += '.bam'


    # sak (fastq to fasta)
    # ~/bin/sak reads.1.fastq -o reads.1.fa
    # ~/bin/sak reads.2.fastq -o reads.2.fa 

    unmappedSortedFasta1 = '{}{}{}_{}_{}_unmapped.sort.1.fa'.format(args.outdir, args.task, readname, accref, donref)
    unmappedSortedFasta2 = None
    if (args.read2fasta is not None):
        unmappedSortedFasta2 = '{}{}{}_{}_{}_unmapped.sort.2.fa'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_sak and (args.b_new or not checkExistence(logger,'Generating FASTA',unmappedSortedFasta1)):
        cmd = ['{}sak'.format(rundir), unmappedSortedFastq1, '-o', unmappedSortedFasta1]
        printAndWrite('Start Swiss-Army-Knife (Sak) Mate 1', 'Start Swiss-Army-Knife (Sak) Mate 1', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Swiss-Army-Knife (Sak) Mate 1. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find Swiss-Army-Knife (Sak) in working directory. Trying global next.')
            cmd = ['sak', unmappedSortedFastq1, '-o', unmappedSortedFasta1]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Swiss-Army-Knife (Sak) Mate 1. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Swiss-Army-Knife (Sak) Mate 1', 'End Swiss-Army-Knife (Sak) Mate 1', logger, 'info')
    if (unmappedSortedFasta2 is not None) and args.b_sak and (args.b_new or not checkExistence(logger,'Generating FASTA',unmappedSortedFasta2)):
        cmd = ['{}sak'.format(rundir), unmappedSortedFastq2, '-o', unmappedSortedFasta2]
        printAndWrite('Start Swiss-Army-Knife (Sak) Mate 2', 'Start Swiss-Army-Knife (Sak) Mate 2', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Swiss-Army-Knife (Sak) Mate 2. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find Swiss-Army-Knife (Sak) in working directory. Trying global next.')
            cmd = ['sak', unmappedSortedFastq2, '-o', unmappedSortedFasta2]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Swiss-Army-Knife (Sak) Mate 2. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Swiss-Army-Knife (Sak) Mate 2', 'End Swiss-Army-Knife (Sak) Mate 2', logger, 'info')


    # gustaf_mate_joining
    # ~/bin/gustaf_mate_joining reads.1.fa reads.2.fa -o reads.joined.fa

    gmj = unmappedSortedFasta1
    if (args.read2fasta is not None):
        gmj = '{}{}{}_{}_{}_unmapped.sort.joined.fa'.format(args.outdir, args.task, readname, accref, donref)
        if args.b_gmj and (args.b_new or not checkExistence(logger,'Gustaf Mate Joining',gmj)):
            cmd = ['{}gustaf_mate_joining'.format(rundir), unmappedSortedFasta1, unmappedSortedFasta2, '-o', gmj]
            printAndWrite('Start Gustaf Mate Joining', 'Start Gustaf Mate Joining', logger, 'info')
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)
            except sp.CalledProcessError:
                printAndWrite('Error in Gustaf Mate Joining. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
            except OSError:
                logger.warning('Could not find Gustaf Mate Joining in working directory. Trying global next.')
                cmd = ['gustaf_mate_joining', unmappedSortedFasta1, unmappedSortedFasta2, '-o', gmj]
                logger.debug(' '.join(cmd))
                try:
                    call(cmd, logger)           
                except sp.CalledProcessError:
                    printAndWrite('Error in Gustaf Mate Joining. Exiting.', 'Exiting.', logger, 'exception')
                    sys.exit(1)
            printAndWrite('End Gustaf Mate Joining', 'End Gustaf Mate Joining', logger, 'info')


    # stellar
    # ~/bin/stellar acc_don.fa reads.joined.fa -o stellar_*.joined.gff -l 30

    stellar = '{}stellar_{}{}_{}_{}.gff'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_stellar and (args.b_new or not checkExistence(logger, 'Stellar', stellar)):
        cmd = ['{}stellar'.format(rundir), candidateRef, gmj, '-o', stellar, '-l', args.stellar_l]
        printAndWrite('Start STELLAR', 'Start STELLAR', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in STELLAR. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find Stellar in working directory. Trying global next.')
            cmd = ['stellar', candidateRef, gmj, '-o', stellar, '-l', args.stellar_l]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in STELLAR. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End STELLAR', 'End STELLAR', logger, 'info')


    # gustaf
    # ~/bin/gustaf acc_don.fa reads.1.fa read.2.fa -m stellar_*.joined.gff -gff gustaf_*.gff -vcf gustaf_*.vcf -st 3 -bth 70 -gth 100 -ith 50 -ll 500 -le 50 -nth 5

    gustaf_vcf = '{}gustaf_{}{}_{}_{}.vcf'.format(args.outdir, args.task, readname, accref, donref)
    gustaf_gff = '{}gustaf_{}{}_{}_{}.gff'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_gustaf and (args.b_new or not checkExistence(logger,'Gustaf',gustaf_vcf,gustaf_gff)):
        cmd = [x for x in ['{}gustaf'.format(rundir), candidateRef, unmappedSortedFasta1, unmappedSortedFasta2, '-m', stellar, '-gff', gustaf_gff, '-vcf', gustaf_vcf, '-st', args.g_st, '-bth', args.g_bth,
            '-gth', args.g_gth, '-ith', args.g_ith, '-ll', args.g_ll, '-le', args.g_le, '-cbp' if args.g_cbp is True else None, '-nth', args.g_nth] if x is not None]
        printAndWrite('Start Gustaf', 'Start Gustaf', logger, 'info')
        logger.debug(' '.join(cmd))
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
            printAndWrite('Error in Gustaf. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        except OSError:
            logger.warning('Could not find Gustaf in working directory. Trying global next.')
            cmd = [x for x in ['gustaf', candidateRef, unmappedSortedFasta1, unmappedSortedFasta2, '-m', stellar, '-gff', gustaf_gff, '-vcf', gustaf_vcf, '-st', args.g_st, '-bth', args.g_bth,
                '-gth', args.g_gth, '-ith', args.g_ith, '-ll', args.g_ll, '-le', args.g_le, '-cbp' if args.g_cbp is True else None, '-nth', args.g_nth] if x is not None]
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)           
            except sp.CalledProcessError:
                printAndWrite('Error in Gustaf. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
        printAndWrite('End Gustaf', 'End Gustaf', logger, 'info')

    # extract single HGT boundaries
    # grep "Trans" gustaf_*.vcf | grep "_2\>"| cut -f 1,2,4,5| sed 's/[][]/\t/g;s/\:/\t/g;s/\t\+/\t/g' > hgt_cand_*.txt

    hgt_cand = '{}hgt_cand_{}{}_{}_{}.txt'.format(args.outdir, args.task, readname, accref, donref)
    if args.b_cand and (args.b_new or not checkExistence(logger, 'Extracting Single Boundaries', hgt_cand)):
        cmd = 'grep \"Trans\" {} | grep \"_2\>\"| sed \'s/;/\\t/g;s/=/\\t/g\'| cut -f 1,2,4,5,15|  sed \'s/[][]/\\t/g;s/\:/\\t/g;s/\\t\+/\\t/g\'| awk \'!seen[$0]++\' >  {}'.format(gustaf_vcf, hgt_cand)
        printAndWrite('Start Extracting Single Boundaries', 'Start Extracting Single Boundaries', logger, 'info')
        logger.debug(cmd)
        try:
            call(cmd, logger, shell=True)
        except sp.CalledProcessError:
            printAndWrite('Error in Extracting Single Boundaries. Exiting.', 'Exiting.', logger, 'exception')
            sys.exit(1)
        printAndWrite('End Extracting Single Boundaries', 'End Extracting Single Boundaries', logger, 'info')

    # Phage screening
    phagesam = None
    if (args.phage_ref is not None):
        # yara indexer of phage genomes
        # ~/bin/yara_indexer Phage_genomes.fa -o Phage_index
        yaraIndex = '{}Phage_index'.format(args.phage_dir)
        if args.b_phage and (args.b_new or not checkExistence(logger, 'Yara Indexer', yaraIndex+'.lf.drp', yaraIndex+'.lf.drs', yaraIndex+'.lf.drv', yaraIndex+'.lf.pst', yaraIndex+'.rid.concat',
             yaraIndex+'.rid.limits', yaraIndex+'.sa.ind', yaraIndex+'.sa.len', yaraIndex+'.sa.val', yaraIndex+'.txt.concat', yaraIndex+'.txt.limits', yaraIndex+'.txt.size')):
            cmd = ['{}yara_indexer'.format(rundir), args.phage_ref, '-o', yaraIndex]
            printAndWrite('Start Yara Indexer Phage database', 'Start Yara Indexer Phage database', logger, 'info')
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)
            except sp.CalledProcessError:
                printAndWrite('Error in Yara Indexer. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
            except OSError:
                logger.warning('Could not find Yara Indexer in working directory. Trying global next.')
                cmd = ['yara_indexer', args.phage_ref, '-o', yaraIndex]
                logger.debug(' '.join(cmd))
                try:
                    call(cmd, logger)           
                except sp.CalledProcessError:
                    printAndWrite('Error in Yara Indexer. Exiting.', 'Exiting.', logger, 'exception')
                    sys.exit(1)
            printAndWrite('End Yara Indexer', 'End Yara Indexer', logger, 'info')

        # yara mapper of phage genomes
        # ~/bin/yara_mapper -o *phage_all.sam reads.1.fasta reads.2.fasta -a -e 10

        phagemapping = '{}Phage_{}.sam'.format(args.outdir, readname)
        if args.b_phage and (args.b_new or not checkExistence(logger, 'Yara Mapping', phagemapping)): 
            cmd = [x for x in ['{}yara_mapper'.format(rundir), yaraIndex, args.read1fasta, args.read2fasta, '-o', phagemapping, '-e', args.yara_e, '-ll', args.yara_ll, '-le', args.yara_le, '-t', args.yara_t] if x is not None]
            printAndWrite('Start Yara Mapping Phage database', 'Start Yara Mapping Phage database', logger, 'info')
            logger.debug(' '.join(cmd))
            try:
                call(cmd , logger)
            except sp.CalledProcessError as err:
                printAndWrite('Error in Yara Mapper. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
            except OSError:
                logger.warning('Could not find Yara Mapper in working directory. Trying global next.')
                cmd = [x for x in ['yara_mapper', yaraIndex, args.read1fasta, args.read2fasta, '-o', phagemapping, '-e', args.yara_e, '-ll', args.yara_ll, '-le', args.yara_le, '-t', args.yara_t] if x is not None]
                logger.debug(' '.join(cmd))
                try:
                    call(cmd, logger)           
                except sp.CalledProcessError:
                    printAndWrite('Error in Yara Mapper. Exiting.', 'Exiting.', logger, 'exception')
                    sys.exit(1)
            printAndWrite('End Yara Mapper', 'End Yara Mapper', logger, 'info')

        # extracting only mapped phage reads
        # ~/bin/samtools view -h -S -F 4 Phage_all.sam > Phage_all_mapped.sam
        phagesam = '{}Phage_{}_mapped.sam'.format(args.outdir, readname)
        if args.b_phage and (args.b_new or not checkExistence(logger,'Sorting Mapping', phagesam)):
            cmd = ['{}samtools'.format(rundir), 'view', '-h', '-F', '4', '-S', phagemapping, '-o', phagesam]
            printAndWrite('Start Samtools extraction of mapped reads', 'Start Samtools extraction of mapped reads', logger, 'info')
            logger.debug(' '.join(cmd))
            try:
                sp.call(cmd)
            except sp.CalledProcessError:
                printAndWrite('Error in Samtools extraction of mapped reads. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
            except OSError:
                logger.warning('Could not find samtools in working directory. Trying global next.')
                cmd = ['samtools', 'view', '-h', '-F', '4', '-S', phagemapping, '-o', phagesam]
                logger.debug(' '.join(cmd))
                try:
                    call(cmd, logger)           
                except sp.CalledProcessError:
                    printAndWrite('Error in Samtools extraction of mapped reads. Exiting.', 'Exiting.', logger, 'exception')
                    sys.exit(1)
            printAndWrite('End Samtools extraction of mapped reads', 'End Samtools extraction of mapped reads', logger, 'info')
 
    # HGT candidate evaluation
    # python ~/bin/hgt_eval.py -o hgt_eval*.vcf  --min_size 500 --max_size 55000 Ecoli_K12_mod_HPylori_1322000-1350000_mod.sort-qname.bam hgt_cand*.txt  "gi|170079663|ref|NC_010473.1|" "gi|766541424|dbj|AP014710.1|"
    #  #!usr/bin/env python as first line in hgt_eval.py => ./hgt_eval.py will run python hgt_eval.py. (Or import it as module.)

    if (args.donor is not None):
        hgteval_vcf = '{}hgt_eval_{}{}_{}_{}.vcf'.format(args.outdir, args.task, readname, accref, donref)
        if args.b_eval and (args.b_new or not checkExistence(logger,'Sorting Mapping', hgteval_vcf)):
            cmd = [x for x in ['python', '{}hgt_eval.py'.format(rundir), mappedSortedQBam, hgt_cand, args.acceptor, args.donor, '--phagefile' if args.phage_ref is not None else None, phagesam,
             '-o', hgteval_vcf, '--min-size', args.h_min, '--max-size', args.h_max, '--tolerance', args.h_tol, '--pair-support' if args.h_pairs is False else None, '--num-boot-regions', args.h_bootnum, '--boot-sens', args.h_bootsens] if x is not None]
            # hgt_eval.func(parameter)
            printAndWrite('Start HGT evaluation', 'Start HGT evaluation', logger, 'info')
            logger.debug(' '.join(cmd))
            try:
                call(cmd, logger)
            except sp.CalledProcessError:
                printAndWrite('Error in HGT evaluation. Exiting.', 'Exiting.', logger, 'exception')
                sys.exit(1)
            printAndWrite('End HGT evaluation', 'End HGT evaluation', logger, 'info')
    else:
        for don_gi, don_id in get_GIs(args.donor2):
            hgteval_vcf = '{}hgt_eval_{}{}_{}_{}.vcf'.format(args.outdir, args.task, readname, accref, don_id)
            if args.b_eval and (args.b_new or not checkExistence(logger,'Sorting Mapping', hgteval_vcf)):
                cmd = [x for x in ['python', '{}hgt_eval.py'.format(rundir), mappedSortedQBam, hgt_cand, args.acceptor, don_gi, '--phagefile' if args.phage_ref is not None else None, phagesam,
                 '-o', hgteval_vcf, '--min-size', args.h_min, '--max-size', args.h_max, '--tolerance', args.h_tol, '--pair-support' if args.h_pairs is False else None, '--num-boot-regions', args.h_bootnum, '--boot-sens', args.h_bootsens] if x is not None]
                # hgt_eval.func(parameter)
                printAndWrite('Start HGT evaluation', 'Start HGT evaluation', logger, 'info')
                logger.debug(' '.join(cmd))
                try:
                    call(cmd, logger)
                except sp.CalledProcessError:
                    printAndWrite('Error in HGT evaluation. Exiting.', 'Exiting.', logger, 'exception')
                    sys.exit(1)
                printAndWrite('End HGT evaluation', 'End HGT evaluation', logger, 'info')

    print ("Total time: ", time.time() - tstart)

if __name__ == '__main__':
    args = parser().parse_args()
    pipeline(args)
