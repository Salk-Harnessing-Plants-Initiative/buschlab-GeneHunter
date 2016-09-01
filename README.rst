genehunter
==========

This is a small program to create and query a SQLITE database containing data from 
the gff files.



Usage
-----

usage: genehunter [-h] {createdb,extractloc,extractagi,hunt} ...

tair database suite

positional arguments:
  {createdb,extractloc,extractagi,hunt}
                        subcommand help
    createdb            create new database from gff files
    extractloc          extract elements by locus
    extractagi          extract elements by AGI
    hunt                run gene hunter on pvals and extract gene information

optional arguments:
  -h, --help            show this help message and exit


Creating a database:
usage: genehunter createdb [-h] --gff GFF [--desc DESC] [--aliases ALIASES]
                            [--sorf SORF] -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF             path to gff file
  --desc DESC           path to functional descriptions file
  --aliases ALIASES     path to gene aliases file
  --sorf SORF           path to sORF file
  -o OUTPUT, --output OUTPUT
                        path to output file

Querying database by locus:
usage: genehunter extractloc [-h] --db DB [--loc1 LOC1] [--loc2 LOC2]
                              [--chr CHR] [-c] [-i] [--file FILE] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --db DB               path to database
  --loc1 LOC1           locus 1 (has different meanings depending on other
                        options)
  --loc2 LOC2           locus 2 (different meanings depending on other
                        options)
  --chr CHR             chromosome
  -c                    centered around loc1, interval = +/-loc2
  -i                    from locus interval. loc1=start, loc2=end
  --file FILE           read positions from file. The file should contain
                        lines with 4 columns each. (chromosome, position,
                        upstream interval,downstream interval). If only 2
                        columns are present, the interval will be pos-loc1 to
                        pos+loc2
  -o OUTPUT, --output OUTPUT
                        path to output file. Will print to stdout if omitted


Querying database by AGI:
usage: genehunter extractagi [-h] --db DB [--agi AGI] [--file FILE]
                              [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --db DB               path to database
  --agi AGI             single AGI
  --file FILE           read AGIs from file. The file must contain one AGI per
                        line.
  -o OUTPUT, --output OUTPUT
                        path to output file. Will print to stdout if omitted


Using gene hunter function:
usage: genehunter hunt [-h] --db DB --dir DIR [--name NAME] [-u UDISTANCE]
                          [-d DDISTANCE] [-P PVALUE_THRESHOLD]
                          [-M MINOR_ALLELE_COUNT] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  --db DB               path to tair10 database
  --dir DIR             directory where to look for pval files
  --name NAME           name identifier for the pval files. Unix like globs
                        are allow (e.g. "name*")
  -u UDISTANCE, --udistance UDISTANCE
                        maximal upstream distance from TSS (default=4000)
  -d DDISTANCE, --ddistance DDISTANCE
                        maximal downstream distance from TSS (default=4000)
  -P PVALUE_THRESHOLD, --pvalue_threshold PVALUE_THRESHOLD
                        SNP p-value threshold (default=1.0e-6)
  -M MINOR_ALLELE_COUNT, --minor_allele_count MINOR_ALLELE_COUNT
                        minor allele count threshold (default=10)
  -o OUTPUT, --output OUTPUT
                        path to output file. Will print to stdout if omitted.

