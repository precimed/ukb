# ukb
Helper tools to pre-process UK Biobank data.

``ukb_helper.py pheno`` functionality allows to look for a specific set of data fields across .csv files from ``<UKBDATA>/phenotypes/Baskets`` folder, and output the result as a .csv file. It also outputs the number of non-null values per field, and counts pairwise overlap. It accepts flexible format, so instead of `` 50-0.0`` it's OK to specify ``50-0`` or ``50`` (in the later case the tool will read ``50-0.0``, ``50-1.0``, and ``50-2.0`` fields). In the case when a field is present in multiple csv files, the file with largest ID will be used. ``--allow-copies`` allows to change this behavior, and output copies of a fields.

You feedback is very welcome! Please submit your feedback via github tickets - I might not have time to address it straightaway, but eventually I'll get to it. Happy UKBanking!

P.S. if you happen to have ukb-related scripts that are handy for data logistics, feel free to upload them "as it is" to this repo. It doesn't matter if they are user-friendly or not, either way it's nice to keep it in one place.

# demo
```
>python /home/oleksandr/precimed/ukb/ukb_helper.py pheno --input "/space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb2*csv" --fields 50 20216 20016 20510 21022 73
8 --out ~/ukbtest
***********************************************************************
* ukb_helper.py: utilities for UK Biobank data
* Version 1.0.0
* (C) 2020 Oleksandr Frei
* Norwegian Centre for Mental Disorders Research / University of Oslo
* GNU General Public License v3  
***********************************************************************
Call:
./ukb_helper.py pheno \
        --out /home/oleksandr/ukbtest \
        --input ['/space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb*csv'] \
        --fields ['50', '20216', '20016', '20510', '21022', '738'] 

Beginning analysis at Wed Feb 26 12:57:36 2020 by oleksandr, host ip24.ucsd.edu
--input mask matches 19 files
WARNING: /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb20712.csv file apears to be empty and will be ignored
WARNING: /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb20711.csv file apears to be empty and will be ignored
WARNING: /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb20713.csv file apears to be empty and will be ignored
WARNING: /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb21700.csv file apears to be empty and will be ignored
--fields file contain 6 fields
Found 149 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb27125.csv (file ID 27125)
Found 212 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26122.csv (file ID 26122)
Found 1399 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26121.csv (file ID 26121)
Found 1525 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb23401.csv (file ID 23401)
Found 5139 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb24843.csv (file ID 24843)
Found 33 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26346.csv (file ID 26346)
Found 49 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26124.csv (file ID 26124)
Found 4470 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb29266.csv (file ID 29266)
Found 5136 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb22124.csv (file ID 22124)
Found 212 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb23402.csv (file ID 23402)
Found 49 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb29060.csv (file ID 29060)
Found 4467 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb27107.csv (file ID 27107)
Found 1585 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26123.csv (file ID 26123)
Found 11 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb28289.csv (file ID 28289)
Found 149 fields in /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb24434.csv (file ID 24434)
Field 50 expands to: 50-0.0, 50-1.0, 50-2.0
Field 20216 expands to: 20216-2.0
Field 20016 expands to: 20016-0.0, 20016-1.0, 20016-2.0
Field 20510 expands to: 20510-0.0
Field 21022 expands to: 21022-0.0
Field 738 expands to: 738-0.0, 738-1.0, 738-2.0
Final list of fields: 50-0.0, 50-1.0, 50-2.0, 20216-2.0, 20016-0.0, 20016-1.0, 20016-2.0, 20510-0.0, 21022-0.0, 738-0.0, 738-1.0, 738-2.0
WARNING: field 50-0.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 50-1.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 50-2.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 20216-2.0 is present in multiple files (22124, 24843, 29266), and will be taken from 29266
WARNING: field 20016-0.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 20016-1.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 20016-2.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 20510-0.0 is present in multiple files (23401, 26123), and will be taken from 26123
WARNING: field 21022-0.0 is present in multiple files (22124, 24843, 27107, 29266), and will be taken from 29266
WARNING: field 738-0.0 is present in multiple files (23401, 26123), and will be taken from 26123
WARNING: field 738-1.0 is present in multiple files (23401, 26123), and will be taken from 26123
WARNING: field 738-2.0 is present in multiple files (23401, 26123), and will be taken from 26123
from /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb29266.csv reading 20016-0.0, 20016-1.0, 20016-2.0, 20216-2.0, 21022-0.0, 50-0.0, 50-1.0, 50-2.0...
done, 502536 subjects, 9 fields found
after merging, combined data so far has 502536 subjects and 9 fields
from /space/gwas-syn1/1/data/GWAS/UKBioBank/phenotypes/Baskets/ukb26123.csv reading 20510-0.0, 738-0.0, 738-1.0, 738-2.0...
done, 502539 subjects, 5 fields found
after merging, combined data so far has 502539 subjects and 13 fields
Saving combined data frame to /home/oleksandr/ukbtest.csv ...
Done.
Analysis finished at Wed Feb 26 12:59:45 2020
Total time elapsed: 2.0m:9.050000000000011s

```

# usage
```
>python ukb_helper.py pheno --help
usage: ukb_helper.py pheno [-h] [--out OUT] [--log LOG] [--log-append]
                           [--input INPUT [INPUT ...]]
                           [--input-list INPUT_LIST]
                           [--fields FIELDS [FIELDS ...]] [--keep KEEP]
                           [--remove REMOVE] [--allow-copies] [--dry-run]
                           [--quote-none] [--skip-counts2]

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             prefix for the resulting files (<out>.csv,
                        <out>.counts1.txt (counting number non-missing
                        values), <out>.counts2.txt (counting number of non-
                        zero values for all pairs of fields), etc)
  --log LOG             filename for the log file. Default is <out>.log
  --log-append          append to existing log file. Default is to erase
                        previous log file if it exists.
  --input INPUT [INPUT ...]
                        A list of input files, or a mask containing a *
                        wildcard for searching input .csv files.
  --input-list INPUT_LIST
                        A text file containing the list of of input files
  --fields FIELDS [FIELDS ...]
                        A text file containing the list of fields to look for.
                        Does not accept wild cards, but for a data field
                        '12224-2.0' it is OK to specify 12224 or 12224-2, to
                        search for all variants of the field. 'eid' column is
                        always added automatically, but it's also acceptable
                        to include 'eid' in the --fields list.
  --keep KEEP           accepts space/tab-delimited text file, without header,
                        with individual IDs in the first column, and removes
                        all unlisted samples from the analysis
  --remove REMOVE       accepts the same sort of file as --keep, and removes
                        all listed subjects
  --allow-copies        When data field is present in multiple input files, by
                        default the field from the file with largest ID is
                        used. When --allow-copies is specified, all data field
                        copies are retained. To avoid ambiguity, we prefix all
                        data field names with the ID of the file that it comes
                        from.
  --dry-run             Just produce the .log file, skipping all actions.
  --quote-none          Sets pandas.to_csv(quoting=QUOTE_NONE). See pandas
                        documentation for more details.
  --skip-counts2        Do not produce <out>.counts2.txt file

```
