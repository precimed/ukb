#!/usr/bin/env python
'''
(c) 2016-2018 Oleksandr Frei
Various utilities for UK Biobank data
'''

import pandas as pd
import numpy as np
import os
import time, sys, traceback
import argparse
import glob
import socket
import getpass
import six
import itertools
import csv

__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* ukb_helper.py: utilities for UK Biobank data\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2020 Oleksandr Frei\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

def check_input_file(file):
    if not os.path.isfile(file):
        raise ValueError("Input file does not exist: {f}".format(f=file))

def check_output_file(file):
    # Create target folder if it doesn't exist
    output_dir = os.path.dirname(file)
    if output_dir and not os.path.isdir(output_dir): os.makedirs(output_dir)  # ensure that output folder exists

def parse_args(args):
    parser = argparse.ArgumentParser(description="A collection of various utilities for GWAS summary statistics.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("--out", type=str, default=None, help="prefix for the resulting files (<out>.csv, <out>.counts1.txt (counting number non-missing values), <out>.counts2.txt (counting number of non-zero values for all pairs of fields), etc)")
    parent_parser.add_argument("--log", type=str, default=None, help="filename for the log file. Default is <out>.log")
    parent_parser.add_argument("--log-append", action="store_true", default=False, help="append to existing log file. Default is to erase previous log file if it exists.")
    
    subparsers = parser.add_subparsers()

    # 'csv' utility : load raw summary statistics file and convert it into a standardized format
    parser_pheno = subparsers.add_parser("pheno", parents=[parent_parser],
        help='Aggregate a list of data files across multiple baskets.')

    parser_pheno.add_argument("--input", type=str, nargs='+', default=[],
        help="A list of input files, or a mask containing a * wildcard for searching input .csv files.")
    parser_pheno.add_argument("--input-list", type=str,
        help="A text file containing the list of of input files")
    parser_pheno.add_argument("--fields", type=str, nargs='+', default=[],
        help="A text file containing the list of fields to look for. Does not accept wild cards, but for a data field '12224-2.0' it is OK to specify 12224 or 12224-2, to search for all variants of the field. 'eid' column is always added automatically, but it's also acceptable to include 'eid' in the --fields list.")
    parser_pheno.add_argument("--keep", type=str, default=None,
        help="accepts space/tab-delimited text file, without header, with individual IDs in the first column, and removes all unlisted samples from the analysis")
    parser_pheno.add_argument("--remove", type=str, default=None,
        help="accepts the same sort of file as --keep, and removes all listed subjects")
    parser_pheno.add_argument("--allow-copies", action="store_true", default=False, help="When data field is present in multiple input files, "
        "by default the field from the file with largest ID is used. When --allow-copies is specified, "
        "all data field copies are retained. To avoid ambiguity, we prefix all data field names with the ID of the file that it comes from.")
    parser_pheno.add_argument("--dry-run", action="store_true", default=False, help="Just produce the .log file, skipping all actions.")
    parser_pheno.add_argument("--quote-none", action="store_true", default=False, help="Sets pandas.to_csv(quoting=QUOTE_NONE). See pandas documentation for more details.")
    parser_pheno.add_argument("--skip-counts2", action="store_true", default=False, help="Do not produce <out>.counts2.txt file")

    parser_pheno.set_defaults(func=make_pheno)

    return parser.parse_args(args)

def make_pheno(args, log):
    check_output_file(args.out)
    if args.keep is not None: check_input_file(args.keep)
    keep = None if (args.keep is None) else pd.read_csv(args.keep, header=None, delim_whitespace=True, usecols=[0], names=['eid'], dtype=str)

    if args.remove is not None: check_input_file(args.remove)
    remove = None if (args.remove is None) else set(pd.read_csv(args.remove, header=None, delim_whitespace=True, usecols=[0], names=['eid'], dtype=str)['eid'].values)

    if (args.input_list is not None) and (len(args.input) > 0):
        raise ValueError("Can not use --input-list together with --input")
    if args.input_list is not None:
        check_input_file(args.input_list)
        args.input = [x.strip() for x in open(args.input_list, 'r').readlines() if (len(x.strip()) > 0)]
        log.log('--input-list file contains {} files'.format(len(args.input)))
    else:
        if len(args.input) == 1:
            args.input = glob.glob(args.input[0])
            log.log('--input mask matches {} files'.format(len(args.input)))
        else:
            log.log('--input contains {} files'.format(len(args.input)))
    args_input_csv = [x for x in args.input if x.endswith(".csv")]
    if (len(args_input_csv) == 0) or (len(args_input_csv) != len(args.input)):
        raise(ValueError('--input is either empty, or have files with extention other than .csv'))
    
    empty_list = []
    for f in args.input:
        check_input_file(f)
        if os.stat(f).st_size == 0:
            log.log("WARNING: {} file apears to be empty and will be ignored".format(f))
            empty_list.append(f)
    args.input = [f for f in args.input if (f not in empty_list)]

    if len(args.fields) == 1:
        if os.path.isfile(args.fields[0]):
            args.fields = [x.strip() for x in open(args.fields[0], 'r').readlines() if (len(x.strip()) > 0)]
    args.fields = [x for x in args.fields if (x != 'eid')]
    log.log('--fields file contain {} fields'.format(len(args.fields)))

    input_to_id_list = []
    for x in args.input:
        try:
            input_to_id_list.append((x, int(os.path.splitext(os.path.basename(x))[0].replace('ukb', ''))))
        except:
            raise(ValueError("Unable to extract an integer ID from {} (e.g., ukb28289.csv converts to 28289)".format(x)))
    input_to_id = dict(input_to_id_list)
    if len(set(input_to_id.values())) != len(input_to_id_list):
        raise ValueError('Duplicated IDs in --input: {}'.format(', '.join([str(y) for x, y in input_to_id_list])))

    # challendge:
    # columns are named as follows: 23170-0.0
    # user-provided fields can look as 23170 or 23170-0
    # how do we figure this out?

    expands_df = None
    cols_df = None
    for f in args.input:
        cols = list(pd.read_csv(f, sep=',', nrows=0).columns)
        log.log('Found {} fields in {} (file ID {})'.format(len(cols), f, input_to_id[f]))
        expands = [[x, x] for x in cols] + [[x.split('.', 1)[0], x] for x in cols] + [[x.split('-', 1)[0], x] for x in cols]
        df = pd.DataFrame(expands, columns=['key', 'val'])
        df['val'] = ['{} '.format(x) for x in df['val'].values]
        expands_df = df if (expands_df is None) else pd.concat([expands_df, df])
        cols=pd.DataFrame({'cols':[x for x in cols if (x != 'eid')]});        cols['file'] = f;        cols['id'] = input_to_id[f]
        cols_df = cols if (cols_df is None) else pd.concat([cols_df, cols])
    expands_df.drop_duplicates(inplace=True)
    df = expands_df.groupby('key').agg({'val':'sum'}).reset_index()
    expands = dict(zip(df['key'].values, [x.split() for x in df['val'].values]))
    cols_df.sort_values(['cols', 'id'], inplace=True)
    cols_df_dedup = cols_df.drop_duplicates(['cols'], keep='last')
    
    for field in args.fields:
        if field not in expands: log.log('WARNING: Field {} not found'.format(field)); continue
        #if (len(expands[field]) == 1) and (expands[field][0] == field): continue
        log.log('Field {} expands to: {}'.format(field, ', '.join(expands[field])))
    args.fields = list(itertools.chain(*[expands[x] for x in args.fields if (x in expands)]))
    log.log("Final list of fields: {}".format(', '.join(args.fields)))
    if len(args.fields) == 0: raise ValueError("no fields found!")
    for field in args.fields:
        if args.allow_copies: break
        if (cols_df['cols'] == field).sum() > 1:
            log.log('WARNING: field {} is present in multiple files ({}), and will be taken from {}'.format(field,
                ', '.join([str(x) for x in cols_df.loc[cols_df['cols'] == field, 'id'].values]),
                ', '.join([str(x) for x in cols_df_dedup.loc[cols_df_dedup['cols'] == field, 'id'].values])))

    # make final list
    cols_df_dedup = pd.merge(cols_df_dedup, pd.DataFrame({'cols':args.fields}), how='inner', on='cols')
    cols_df = pd.merge(cols_df, pd.DataFrame({'cols':args.fields}), how='inner', on='cols')

    df_merged = None
    for f in args.input:
        usecols = list(cols_df.loc[cols_df['file'] == f, 'cols'].values) if args.allow_copies else list(cols_df_dedup.loc[cols_df_dedup['file'] == f, 'cols'].values)
        if len(usecols) == 0: continue

        log.log('from {} reading {}...'.format(f, ', '.join(usecols)))
        if args.dry_run: continue
        df = pd.read_csv(f, sep=',', usecols=['eid'] + usecols, dtype=str, encoding= 'unicode_escape')
        if args.allow_copies: df.columns = [('{}_ukb{}'.format(x, input_to_id[f]) if (x != 'eid') else 'eid') for x in df.columns]
        log.log('done, {} subjects, {} fields found'.format(df.shape[0], df.shape[1]))
        
        if keep is not None:
            df = pd.merge(df, keep, how='inner', on='eid')
            log.log('keep {} subjects due to --keep'.format(df.shape[0], df.shape[1]))

        if remove is not None:
            df['remove'] = [(x in remove) for x in df['eid'].values]
            df = df[~df['remove']].copy()
            log.log('keep {} subjects due to --remove'.format(df.shape[0], df.shape[1]))

        if len(df) == 0: continue

        df_merged = df if (df_merged is None) else pd.merge(df_merged, df, how='outer', on='eid')
        log.log('after merging, combined data so far has {} subjects and {} fields'.format(df_merged.shape[0], df_merged.shape[1]))

    if args.dry_run:
        log.log('Done (--dry-run)')
        return

    log.log('Saving combined data frame to {}.csv ...'.format(args.out))
    df_merged.to_csv('{}.csv'.format(args.out), sep=',', index=False, quoting=(csv.QUOTE_NONE if args.quote_none else csv.QUOTE_MINIMAL))

    df_merged.notnull().sum().reset_index().to_csv('{}.counts1.txt'.format(args.out), sep='\t', header=None, index=False)

    if not args.skip_counts2:
        v=df_merged.notnull().astype(float).values; m = np.dot(np.transpose(v), v); 
        pd.DataFrame(m, columns=df_merged.columns, index=df_merged.columns).astype(int).to_csv('{}.counts2.txt'.format(args.out), sep='\t')

    log.log('Done.')

### =================================================================================
###                                Misc stuff and helpers
### =================================================================================
def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = six.moves.reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh, mode):
        self.fh = fh
        self.log_fh = open(fh, mode) if (fh is not None) else None

        # remove error file from previous run if it exists
        try:
            os.remove(fh + '.error')
        except OSError:
            pass

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            self.log_fh.flush()

    def error(self, msg):
        '''
        Print to log file, error file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            with open(self.fh + '.error', 'w') as error_fh:
                error_fh.write(str(msg).rstrip() + '\n')

### =================================================================================
###                                Main section
### ================================================================================= 
if __name__ == "__main__":
    args = parse_args(sys.argv[1:])

    if args.out is None:
        raise ValueError('--out is required.')

    log = Logger(args.log if args.log else (args.out + '.log' if (args.out != '-') else None), 'a' if args.log_append else 'w')
    start_time = time.time()

    try:
        defaults = vars(parse_args([sys.argv[1]]))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
        header = MASTHEAD
        header += "Call: \n"
        header += './ukb_helper.py {} \\\n'.format(sys.argv[1])
        options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
        header += '\n'.join(options).replace('True','').replace('False','')
        header = header[0:-1]+'\n'
        log.log(header)
        log.log('Beginning analysis at {T} by {U}, host {H}'.format(T=time.ctime(), U=getpass.getuser(), H=socket.gethostname()))

        # run the analysis
        args.func(args, log)

    except Exception:
        log.error( traceback.format_exc() )
        raise

    finally:
        log.log('Analysis finished at {T}'.format(T=time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
