#!/usr/bin/perl 

use strict;
use GD::Graph::lines;
use Bit::Vector;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use constant version => "1.2.1";
#reminder; must change version in POD manually
use constant default_sliding_win_size => 15;
$|=1;

# Release notes
# 1.0 -Original release
# 1.1 -Changed graph size, line colour and fonts for improved readability
#     -Omit chromosome from analysis if no reads for that chromosome
#     -Output fragment length estimate on command line
# 1.2 -Nov 29, 2012
#     -Corrections to Pod manual
#     -Validation of mappability filename before parsing
# 1.2.1 -Feb 13, 2013  
#     -updating docs/readme clarifying mapability
#     -Fix docs on --chrom-length-file
# Code Written by: G. Palidwor

sub sliding_average($$){
	my $array_ref = shift;
	my $window = shift;
	my $winsize = 2*$window+1;
	my @avg;
	for(my $i=0;$i<scalar(@$array_ref);$i++){
		my $i1;
		my $i2;
		$i1=$i-$window;
		if($i1<0){
			$i1 = 0;
		}
		$i2=$i+$window;
		if($i2>=(scalar(@$array_ref))){
			$i2 = scalar(@$array_ref)-1;
		}
		my $val = sum(@$array_ref[$i1..$i2])/($i2-$i1+1);
		push(@avg,$val);
	}
	return(@avg);
}


sub max_index($){
	my $array_ref = shift;
	my $max_index = 0;
	$$array_ref[$max_index] > $$array_ref[$_] or $max_index = $_ for 1 .. scalar(@$array_ref); 
	return $max_index;
}
my $read_len;

my $verbose = 0;
my $mappability_path;
my $chrom_length_file="";
my $input_bed;
my $output_file_prefix;
my $smooth_win_size=15;
my $min_shift=0;
my $max_shift=400;
my $help = 0;
GetOptions(	"verbose!"=>\$verbose,
				"mappability_path=s"=>\$mappability_path,
				"chrom_length_file:s"=>\$chrom_length_file,
				"input_bed=s"=>\$input_bed,
				"prefix=s"=>\$output_file_prefix,
				"smooth_win_size:i"=>\$smooth_win_size,
				"min_shift:i"=>\$min_shift,
				"max_shift:i"=>\$max_shift,
				"help|?"=>\$help
				) or pod2usage(-verbose=>0,-exitvalue=>1);

pod2usage(-verbose=>2,-exitval=>2) if $help;
if(!defined($mappability_path) || !defined($input_bed) || !defined($output_file_prefix)){
	print STDERR "Error: --mappability_path not specified\n" unless defined($mappability_path);
	print STDERR "Error: --input_bed not specified\n" unless defined($input_bed);
	print STDERR "Error: --prefix not specified\n"  unless defined($output_file_prefix);
   pod2usage(-verbose=>0,-exitval=>2);
}elsif($min_shift >= $max_shift){
 	print STDERR "Error: --min_shift must be less than --max_shift.\n";
   print STDERR "Error: Currently --min_shift=$min_shift and --max_shift=$max_shift\n";
   pod2usage(-verbose=>0,-exitval=>2);
}elsif($smooth_win_size<1){
	print STDERR "Error: --smooth_win_size must be greater than zero.\n";
   pod2usage(-verbose=>0,-exitval=>2);
}elsif($chrom_length_file eq ""){
	print STDERR "Error: --chrom_length_file must be specified\n";
	pod2usage(-verbose=>0,-exitval=>2);
}elsif(!(-d $mappability_path)){
	print STDERR "Error: Provided mappability path is not a directory $mappability_path\n";
   pod2usage(-verbose=>0,-exitval=>2);
}elsif(!(-e $mappability_path)){
	print STDERR "Error: Provided mappability path is not readable $mappability_path\n";
   pod2usage(-verbose=>0,-exitval=>2);
}elsif(!opendir(MAPPABILITY_DIR, $mappability_path)){
	print STDERR "Error: can't open directory --mappability_path $mappability_path\n";
   pod2usage(-verbose=>0,-exitval=>2);
}


my @map_files = readdir(MAPPABILITY_DIR); 

#get the chromosome lengths
my %chrom_length_hash;
open CHROM_LENGTHS, "<$chrom_length_file" or die "can't open chromosome length file: '$chrom_length_file': $!\n";
while(<CHROM_LENGTHS>){
   my ($chrom,$chrlen) = split(/\s+/);
	$chrom_length_hash{$chrom} = $chrlen;
}

#get the mapability
my %chr_map_hash;
my %pos_strand_hash;
my %neg_strand_hash;
my $total_genome_len = 0;
my $shift_range=$max_shift - $min_shift;

foreach my $map_file (@map_files){
	if($map_file eq "." || $map_file eq ".." || $map_file !~ /([^\._]*)\.map$/){
		next;
	}
	#parse chromosome identifier from mapability file name	
   $map_file =~ /([^\._]*)\.map$/;
	my $key = $1;
	if($verbose){
		print "Parsing mapability file: $map_file\t$key\n";
	}
	if($key eq ""){
		die "Fatal Error: unable to determine chromosome identifier for $map_file.\n" . 
          "Please review MaSC documentation for correct filename format.\n";
	}
   if(!exists($chrom_length_hash{$key})){
		die "Fatal Error: chromsome identifier $key from $map_file does not match any " .
          "chromosome length records. \n";	
	}

	$total_genome_len += $chrom_length_hash{$key};

   my $pos_vector = Bit::Vector->new($chrom_length_hash{$key}+$shift_range);
	$pos_strand_hash{$key} = \$pos_vector;

   my $neg_vector = Bit::Vector->new($chrom_length_hash{$key}+$shift_range);
	$neg_strand_hash{$key} = \$neg_vector;

   my $vector = Bit::Vector->new($chrom_length_hash{$key}+$shift_range);
  
	open CHR_MAP, "<$mappability_path/$map_file" or die "Can't open $mappability_path/$map_file: $!\n";

	#skip the first line which is a wiggle track description
	readline(*CHR_MAP);

	while(<CHR_MAP>){
		my ($chr,$start,$end,$score) = split(/\s+/);
		$vector->Interval_Fill($start,$end);
	}
	$chr_map_hash{$key}=\$vector;
}

open INPUT_BED, "<$input_bed" or die "can't open $input_bed: $!";
my $bed_line_num = 0;
if($verbose){
	print "Parsing bed file: $input_bed\n";
}

my %unknown_chrom;
while(<INPUT_BED>){
	$bed_line_num++;
	chomp;
	my ($chr,$start,$end,$id,$val,$strand) = split(/\t/);
	if($read_len == 0){
		$read_len = $end - $start;
	}
	if(!defined($neg_strand_hash{$chr})){
		#warn "Warning: unknown chromosome identifier '$chr' in BED file $input_bed line number $bed_line_num\n";
		if(!exists($unknown_chrom{$chr})){
			$unknown_chrom{$chr}=0;
		}
		$unknown_chrom{$chr}++;
		next;
	}
	if($strand eq "-"){
		${$neg_strand_hash{$chr}}->Bit_On($end);
	}
	if($strand eq "+"){
		${$pos_strand_hash{$chr}}->Bit_On($start);
	}
}


my $vlen = $max_shift - $min_shift + 1;
my @r1;$r1[$vlen]=0;
my $total_pos_reads=0;
my $total_neg_reads=0;

my @total_mdr_count=0;$total_mdr_count[$vlen]=0;
my @total_fgmdr_count=0;$total_fgmdr_count[$vlen]=0;
my @total_pos_mdr=0;$total_pos_mdr[$vlen]=0;
my @total_neg_mdr=0;$total_neg_mdr[$vlen]=0;

foreach my $chr (keys %neg_strand_hash){
	if($verbose){ 
  		print "Processing reads for: $chr\n";
	}
	#check for chromsomes with no reads, remove them from
	#the calculation
	if(${$pos_strand_hash{$chr}}->Norm()==0 && ${$pos_strand_hash{$chr}}->Norm()==0){
		if($verbose){
			print "Skipping calculations for $chr: no mapped reads found\n";
		}
		next;
	}
	my $Mdr = ${$pos_strand_hash{$chr}}->Shadow();
	my $working_vector2 = $Mdr->Shadow();
	my $working_vector3 = $Mdr->Shadow();
	my $working_vector4 = $Mdr->Shadow();
	my $neg_map_vector = ${$chr_map_hash{$chr}};
	my $pos_map_vector = $neg_map_vector->Clone();

	#shift the end of the read to the start for MaSC
	for(my $i=0;$i<$read_len;$i++){
		$neg_map_vector->rotate_left();
	}

	my $neg_vector = $neg_strand_hash{$chr};
 	my $pos_vector = $pos_strand_hash{$chr};

	$total_neg_reads+=${$neg_vector}->Norm();
	$total_pos_reads+=${$pos_vector}->Norm();

	for(my $i=0;$i<$min_shift;$i++){
		#rotate vector
		${$neg_vector}->rotate_right();
		$neg_map_vector->rotate_right();
	}

	for(my $d=0;$d<=($max_shift-$min_shift);$d++){
		#neg_pos map overlap
		$working_vector2->And($neg_map_vector,$pos_map_vector);
		$total_mdr_count[$d] += $working_vector2->Norm();
		
		#neg pos strand overlap
		$working_vector3->And(${$neg_vector},${$pos_vector});
		$r1[$d] += $working_vector3->Norm();

		#strand and map overlap
		$working_vector4->And($working_vector3,$working_vector2);	
		$total_fgmdr_count[$d] += $working_vector4->Norm();

		#neg strand mdr
		$working_vector4->And(${$neg_vector},$working_vector2);
		$total_neg_mdr[$d] += $working_vector4->Norm();
			
		#pos strand mdr
		$working_vector4->And(${$pos_vector},$working_vector2);
		$total_pos_mdr[$d] += $working_vector4->Norm();
	
		#rotate vector
		${$neg_vector}->rotate_right();
		$neg_map_vector->rotate_right();
	}
}

my $neg_mean=$total_neg_reads/$total_genome_len;
my $pos_mean=$total_pos_reads/$total_genome_len;
my $neg_var=$neg_mean-$neg_mean**2;
my $pos_var=$pos_mean-$pos_mean**2;

my $csv_file = $output_file_prefix . "_MaSC.txt";
my $png_file = $output_file_prefix . "_MaSC.png";
open TXT, ">$csv_file" or die "can't open $csv_file: $!\n";
print TXT "d\tCorrelation\tMean Correlation\tMaSC\tMean MaSC\n";
my @corr_val;
my @masc_val;
my @d_val;
for(my $d=0;$d<$vlen;$d++){
	my $pos_mean_mdr=$total_pos_mdr[$d]/$total_mdr_count[$d];
	my $neg_mean_mdr=$total_neg_mdr[$d]/$total_mdr_count[$d];
	my $pos_var_mdr=$pos_mean_mdr-$pos_mean_mdr**2;
	my $neg_var_mdr=$neg_mean_mdr-$neg_mean_mdr**2;
	my $corr=(($r1[$d]/$total_genome_len - $neg_mean*$pos_mean)/sqrt($neg_var*$pos_var));
	my $masc=($total_fgmdr_count[$d]/$total_mdr_count[$d]  - $neg_mean_mdr*$pos_mean_mdr)/sqrt($neg_var_mdr*$pos_var_mdr);
	push(@d_val,$min_shift+$d);
	push(@corr_val,$corr);
	push(@masc_val,$masc);
}

my @avg_corr = sliding_average(\@corr_val,$smooth_win_size);
my $max_corr_index = $d_val[max_index(\@avg_corr)];
my @avg_masc = sliding_average(\@masc_val,$smooth_win_size);
my $max_masc_index = $d_val[max_index(\@avg_masc)];

for(my $i=0;$i<$vlen;$i++){
	print TXT "$d_val[$i]\t$corr_val[$i]\t$avg_corr[$i]\t$masc_val[$i]\t$avg_masc[$i]\n";
}
close TXT;
my @data = (
	\@d_val,
	\@corr_val,
	\@avg_corr,
	\@masc_val,
	\@avg_masc
);

my $my_graph = GD::Graph::lines->new(800,600);
$my_graph->set(
	x_label => 'Shift',
	y_label => 'Cross Correlation',
	x_label_skip => 50,
	title => $output_file_prefix,
	legend_placement => 'RM',
	bgclr => 'white',
	transparent => '0',
	r_margin=>30,
	l_margin=>10,
	t_margin=>10,
	b_margin=>10,
	x_label_position => 0.5,
	dclrs=>[qw(pink red #00FFFF #0000FF)]
);
$my_graph->set_x_axis_font(GD::Font->MediumBold);
$my_graph->set_y_axis_font(GD::Font->MediumBold);
$my_graph->set_x_label_font(GD::Font->Large);
$my_graph->set_y_label_font(GD::Font->Large);
$my_graph->set_title_font(GD::Font->Giant);
$my_graph->set_legend_font(GD::Font->Small);
	
$my_graph->set_legend(("Naive FL:" . $max_corr_index,"Avg Naive","MaSC FL:" . $max_masc_index,"Avg MaSC"));

my $gd = $my_graph->plot(\@data) or die $my_graph->error;
open(IMG, ">$png_file") or die $!;
binmode IMG;
print IMG $gd->png;

if(scalar(keys(%unknown_chrom))>0){
	print("\n");
	print("Warning! The following unknown chromosome identifiers were found and ignored:\n");
	foreach my $uchrom (sort(keys(%unknown_chrom))){
		print($uchrom . ": " . $unknown_chrom{$uchrom} . " instances\n");
	}
	print("\n")
}

print("\n");
print("Naive Fragment Length:" . $max_corr_index . "\n");
print("MaSC Fragment Length:" . $max_masc_index . "\n")
__END__

=head1 DESCRIPTION

This document describes the usage of the Perl-based Mappability-Sensitive Cross-Correlation (MaSC) software that accompanies the paper:

   Ramachandran P., Palidwor G., Porter C.J., Perkins T.J., "MaSC: 
   Mappability-Sensitive Cross-Correlation for Estimating Mean Fragment Length 
   of Single-End Short Read Sequencing Data", Bioinformatics, Volume 29, 
   Issue 4, pp. 444-450, 2013.

Please cite this publication if you use the described algorithm or this software.

=head1 NAME 

MaSC.pl Version 1.2.1: A Perl-based reference implementation of Mappability-Sensitive Cross-Correlation (MaSC).

Given a set of mapped reads and corresponding mappability files, this software provides fragment-length estimation using both naive and mappability-adjusted (MaSC) cross-correlation. Please see the above paper for further details about the methods.

=head1 SYNOPSIS

MaSC.pl --help --verbose --mappability_path=/mappability/dir/ --chrom_length_file=/dir/chrom_lens.txt --input_bed=/dir/input.bed --prefix=myprefix --smooth_win_size=15 --min_shift=0 --max_shift=400

=head1 OPTIONS

Optional:

--verbose,-h

This switch turns on verbose output which outputs program status while it is running.

--help,-h

Prints detailed help and terminates

=item B<Required:>

=item B<--mappability_path,-ma>

The directory containing the organism and read-length appropriate mappability wiggle files.  These files should indicate only positions where the genome is uniquely mappable, see the example files provided with the test data available at http://www.perkinslab.ca/pubs/RPPP2012.html . This data can be generated using the UCSC Table Browser, extracting the Mapability track for the genome of interest and filtering for Mapability=1. These file names should be in the form <prefix>_<chromosome_identifier>.map e.g. hg19_36mer_chr10.map where the identifiers correspond to those contained in both the CHROMOSOME_LENGTHS_FILE and the READ_BED_FILE. Only the chromosomes referenced in MAPPABILITY_DIRECTORY files will be used to calculate the genomic length from the CHROMOSOME_LENGTHS_FILE (option --chrom_length_file, see below).

=item B<--chrom_length_file,-ch>

A tab-delimited file containing the lengths of each of the chromosomes. The first column should contain the chromosome identifiers and the second should contain the chromosome lengths.

=item B<--input_bed,-i>

A bed file containing the reads. The chromosomes referenced should correspond to those in the MAPPABILITY_DIRECTORY files.

=item B<--prefix,-p>

This determines the names of the files that will be output in the directory in which MaSC.pl is run.  The software outputs a tab delimited table (<PREFIX>_MaSC.txt) and PNG figure (<PREFIX>_MaSC.png) of the naive and MaSC correlation values as a function of shift (d). 

=item B<Optional:>

=item B<--smooth_win_size,-s>

One half of the size of the window to use in the sliding average calculations; For a value of n specified here, the total window size will be 2*n+1. Default is 15.

=item B<--min_shift,-min>

Minimum shift (d) to calculate. Default is 0.

=item B<--min_shift,-min>

Maximum shift (d) to calculate. Default is 400.

=head1 CONTACT INFORMATION

Address bug reports and comments to:
theodore.j.perkins@gmail.com

=cut
