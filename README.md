## MaSC: Mappability-sensitive cross-correlation
MaSC.pl Version 1.2.1: A Perl-based reference implementation of Mappability-Sensitive Cross-Correlation (MaSC) for estimating mean fragment length of single-end short-read sequencing data

### Citation
This repo contains software accompanying the paper ["MaSC: Mappability-Sensitive Cross-Correlation to Estimate Fragment Length for Short Read Sequencing Data"](https://academic.oup.com/bioinformatics/article/29/4/444/200320), _Bioinformatics,_ 2012. This software is written in Perl and is available free for non-commerical use.

Please cite the above publication if you use the described algorithm or this software.

### Description
Given a set of mapped reads and corresponding mappability files, this software provides fragment-length estimation using both naive and mappability-adjusted (MaSC) cross-correlation. Please see the above paper for further details about the methods.

### Synopsis

```console
    MaSC.pl --help --verbose --mappability_path=/mappability/dir/
    --chrom_length_file=/dir/chrom_lens.txt --input_bed=/dir/input.bed
    --prefix=myprefix --smooth_win_size=15 --min_shift=0 --max_shift=400
```

OPTIONS
    Optional:

    --verbose,-h

    This switch turns on verbose output which outputs program status while
    it is running.

    --help,-h

    Prints detailed help and terminates

    Required:
    --mappability_path,-ma
    The directory containing the organism and read-length appropriate
    mappability wiggle files. These files should indicate only positions
    where the genome is uniquely mappable, see the example files provided
    with the test data available at
    http://www.perkinslab.ca/pubs/RPPP2012.html . This data can be generated
    using the UCSC Table Browser, extracting the Mapability track for the
    genome of interest and filtering for Mapability=1. These file names
    should be in the form <prefix>_<chromosome_identifier>.map e.g.
    hg19_36mer_chr10.map where the identifiers correspond to those contained
    in both the CHROMOSOME_LENGTHS_FILE and the READ_BED_FILE. Only the
    chromosomes referenced in MAPPABILITY_DIRECTORY files will be used to
    calculate the genomic length from the CHROMOSOME_LENGTHS_FILE (option
    --chrom_length_file, see below).

    --chrom_length_file,-ch
    A tab-delimited file containing the lengths of each of the chromosomes.
    The first column should contain the chromosome identifiers and the
    second should contain the chromosome lengths.

    --input_bed,-i
    A bed file containing the reads. The chromosomes referenced should
    correspond to those in the MAPPABILITY_DIRECTORY files.

    --prefix,-p
    This determines the names of the files that will be output in the
    directory in which MaSC.pl is run. The software outputs a tab delimited
    table (<PREFIX>_MaSC.txt) and PNG figure (<PREFIX>_MaSC.png) of the
    naive and MaSC correlation values as a function of shift (d).

    Optional:
    --smooth_win_size,-s
    One half of the size of the window to use in the sliding average
    calculations; For a value of n specified here, the total window size
    will be 2*n+1. Default is 15.

    --min_shift,-min
    Minimum shift (d) to calculate. Default is 0.

    --min_shift,-min
    Maximum shift (d) to calculate. Default is 400.
    
