#!/usr/bin/perl -w

# Description
# Extract scars from paired-end fastqs.
# The protocol generates reads where R1 has a 6-base barcode and 6-base UMI, R2 the scar. 
# Before read extraction, we take the UMI and barcode-information from the R1 and append
# it to the R2, and map the R2. The R2 is forward mapped.

# Input and usage:
# ./scar_CIGAR_sc_v2.pl -R1 [input R1]
#						-R2 [input R2]
#						-op [output prefix]
#						-t  [threads for mapping]
#						-r  [reference file to map to]
#						-bc [barcode list]
#                 		-g  [requested gene]
#				  		-k  [amount of bases in R2 to consider (starting from the primer)]

# Output
# A fastq-file named [prefix].fastq that combines the UMI and barcode-information from 
# R1 (in the readname field) together with the R2 sequence and quality.
# A samfile named [prefix].sam that results from mapping this fastq to the reference.
# A file named [prefix]_scars.txt that contains all scars found in the samfile, one line
# per scar.
# A count table named [prefix]_cscars_coutc.txt in which all scar reads are assigned
# to their barcodes - this is before doing a UMI-correction.
# A count table named [prefix]_cscars_couta.txt in which all scar abundances are assigned
# to their barcodes - this is after doing a UMI-correction and this file should be used
# for all downstream analysis.
# A statistics file named [prefix]_stats.txt with the total amount of reads, the number
# of reads mapped, and the number of reads over applied thresholds (thresholds are
# currently hardcoded but can be changed in the code below).

# Created by: B. Spanjaard
# Created on: 23/4/2015.
# Last update: 4/10/2016

# Dependencies
use strict;
use warnings;
use Getopt::Long;
use POSIX;

# Globals - could change this to input if necessary
my $barcodelength = 16;
my $UMIlength = 10;
my $barcodestart = 0;
my $UMI_number = 4 ** $UMIlength;
my $min_CIGAR_reads = 20;
my $min_CIGAR_fraction = 0.05;
my $protocol = 4;

my $include_barcode_rewrite = 1;
my $include_map = 1;
my $include_barcode_read = 1;
my $include_write_scars = 1;


# Subs
sub FindEndLocation
{
    # Given the position of the leftmost matching base and the cigar score,
    # returns the genomic location of the rightmost base of the read.
    
    my $position = $_[0];
    my $cigar = $_[1];

    my @cigar_numbers = split(/[A-Z]/, $cigar);
    # Calculate the total number of bases in the cigar. This includes deleted
    # bases. Even though these are not in the read, they are in the genome
    # and should therefore be incorporated in the genomic position.
    my $frag_length_w_D = 0;
    foreach my $cigar_number (@cigar_numbers)
    {
	$frag_length_w_D += $cigar_number;
    }
    
    my @cigar_chars = split(/\d+/, $cigar);
    # The number of mismatches before the leftmost matching base should not be counted
    # in the length of the whole read. The number of mismatches after the rightmost
    # matching base should also be subtracted from the read.
    my $s_count = 0;
#	if cigar starts with S, add number belonging to first S; then, if cigar ends with S,
#	add number belonging to last S.
#	Pseudo: while $cigar =~ /([0-9]*S)/g{@x = split("S", $1); $s_count = $s_count + $x[0];}

    if($cigar_chars[1] eq "M")
    {
	$s_count = 0;
    } elsif ($cigar_chars[1] eq "S")
    {
	my @cigar_S_split = split("S", $cigar);
	$s_count = $cigar_S_split[0];
    } else
    {
	print("Unexpected cigar! First character $cigar_chars[1].\n");
    }
    my $start_location = $position - $s_count;

    my @i_split = split("I", $cigar);
    pop(@i_split);
    # The number of insertions does not have influence on the genomic end
    # location, but they do take up bases in the read.
    my $total_i_count = 0;
    foreach my $part (@i_split)
    {
	my ($i_count) = $part =~ /(\d+)$/;
    # print(" $i_count");
	$total_i_count += $i_count;
    }
    
    my $end_location = $start_location + $frag_length_w_D - 1 - 
	$total_i_count;
    
    return $end_location;
}

sub GetBarcode
{
	my $R1_bc_UMI = $_[0];
	
	my $barcode = "";
	$barcode = substr($R1_bc_UMI, $UMIlength + 5, $barcodelength);
	
	return $barcode;
}

sub GetUMI
{
	my $R1_bc_UMI = $_[0];
	
	my $UMI = "";
	$UMI = substr($R1_bc_UMI, 5, $UMIlength);
	
	return $UMI;
}


# Get input parameters and open files
my $in_R1_name = "";
my $in_R2_name = "";
my $out_prefix = "";
my $threads = "";
my $reference = "";
my $barcode_name = "";
my $gene = "";
my $R2_keep = "";
my $min_UMI_count = 1;
my $length = "";
my $primer_sequence = "";
GetOptions
(
	'R1=s' => \$in_R1_name,
	'R2=s' => \$in_R2_name,
	'op=s' => \$out_prefix,
	't=s' => \$threads,
	'r=s' => \$reference,
    'bc=s' => \$barcode_name,
    'g=s' => \$gene,
    'k=i' => \$R2_keep,
    'mU=i' => \$min_UMI_count,
    'l=i' => \$length,
    'ps=s' => \$primer_sequence
);
open(my $R1_file, $in_R1_name)
	or die "Cannot find $in_R1_name\n";
open(my $R2_file, $in_R2_name)
	or die "Cannot find $in_R2_name\n";
open(my $fastq_out_file, ">".$out_prefix.".fastq");
open(my $barcodefile, $barcode_name)
    or die "Cannot find $barcode_name\n";

# Write R2 information and R1 barcode and UMI to new fastq
if ($include_barcode_rewrite == 1)
{
while(my $R1_header = <$R1_file>)
{
	chomp($R1_header);
	my $R1_sequence = <$R1_file>;
	chomp($R1_sequence);
	my $R1_plus = <$R1_file>;
	my $R1_qual = <$R1_file>;
	
	my $R1_UMI_bc = substr($R1_sequence, $barcodestart, ($barcodelength + $UMIlength));
	
	my $R2_header = <$R2_file>;
	chomp($R2_header);
	my @R2_header_fields = split(" ", $R2_header);
	my $R2_out_header = $R2_header_fields[0]." BC:Z:".$R1_UMI_bc;
	print $fastq_out_file "$R2_out_header\n";
	
	my $R2_sequence = <$R2_file>;
	chomp($R2_sequence);
	my $R2_out = substr($R2_sequence, 0, $R2_keep);
	print $fastq_out_file "$R2_out\n";
	my $R2_plus = <$R2_file>;
	chomp($R2_plus);
	print $fastq_out_file "$R2_plus\n";
	my $R2_qual = <$R2_file>;
	chomp($R2_qual);
	my $R2_qual_out = substr($R2_qual, 0, $R2_keep);
	print $fastq_out_file "$R2_qual_out\n";
}
}
close $fastq_out_file;

# Map new R2
my $sam_name = $out_prefix.".sam";
if ($include_map == 1)
{
my $map_command = "bwa mem -t ".$threads." -v 0 -C ".$reference." ".$out_prefix.".fastq > ".$sam_name;
system($map_command) == 0 or die "Could not execute ".$map_command."\n";
}

# Set allowed barcodes
my %bc_hash = ("Unknown", 0);
my $max_barcode = 0;
if ($include_barcode_read == 1)
{
while(my $barcodeline = <$barcodefile>)
{
    chomp($barcodeline);
    my @barcodefields = split(" ", $barcodeline);
    my $barcode = $barcodefields[1];
    $bc_hash{$barcode} = $barcodefields[0];
    if($bc_hash{$barcode} > $max_barcode)
    {
	$max_barcode = $bc_hash{$barcode};
    }
    # print "Read in barcode $bc_hash{$barcode}: $barcode\n";
}
}
close($barcodefile);

# Read samfile and write CIGAR-scars to output-file
my %CIGAR_hash = (); # Counts all reads in a CIGAR including duplicates, keys are CIGARs. 
my %scar_UMI_hash = (); # Counts all reads including duplicates, keys are 
# CIGAR_sequence_location, barcode and UMI.
my %scar_hash = (); # Counts all reads including duplicates, keys are CIGAR_sequence_location and
# barcode.
my %scar_UMI_unique_hash = (); # Counts all unique reads, keys are
# CIGAR_sequence_location and barcode.
my $reads = 0;
my $mapped_reads = 0;
my $mapped_reads_barcode = 0;
my $mapped_reads_length = 0;
my %primer_hash = ($primer_sequence => 1,);
my $primer_length = length($primer_sequence);
my $mapped_reads_length_with_PCRp =0;
if ($include_write_scars == 1)
{
my $scarfile_out_name = $out_prefix."_scars.txt";
open (my $scarfile_out, ">".$scarfile_out_name);
open(my $samfile_in, $sam_name)
    or die "Cannot find $sam_name\n";
print $scarfile_out "CIGAR\tBarcode\tUMI\tLeft location\tFlag\tSNP-code\tSequence\n";
while (my $samline = <$samfile_in>)
{
    # Find read
    chomp($samline);
    if(substr($samline, 0 , 1) eq "\@")
    {
	next;
    }
    my @R1_fields = split(" ", $samline);
    if($R1_fields[1] > 256)
    {
	next;
    }
    $reads++;

    # Search for the barcode-UMI comment.
    my $bc_UMI_index = 13;
    while(substr($R1_fields[$bc_UMI_index], 0, 5) ne "BC:Z:")
    {
    	$bc_UMI_index++;
    }
    
    my $barcode = GetBarcode($R1_fields[$bc_UMI_index]);
    my $UMI = GetUMI($R1_fields[$bc_UMI_index]);
    
    # Filter mapped reads
    if(!($R1_fields[2] eq $gene))
    {
    	next;
    }
    $mapped_reads++;
    
    # Filter reads with correct barcode and no "N" in UMI.
    if(!(@bc_hash{$barcode}) | $UMI =~ "N")
    {
		next;
    }
    $mapped_reads_barcode++;
    
    # Filter reads with correct length.
    my $read_length = "";
   	$read_length = length($R1_fields[9]);
    if($read_length != $length)
    {
    	next;
    }
    $mapped_reads_length++;
    
    # Filter reads with correct PCR-primer.
    # Old: assume primer is at beginning of mapped read.
    #my $read_primer = substr($R1_fields[9], 0, $primer_length);
    # New: determine whether read is mapped forward or reverse.
    my $map_direction = int($R1_fields[1]/16) % 2; # 0 if forward (does not contain a 16), 1 if reverse (contains a 16).
    my $read_primer = substr($R1_fields[9], 0, $primer_length);
    # If forward: PCR-primer is at beginning of mapped read.
    # If reverse: PCR-primer is at beginning of reverse-complemented mapped read.
    if ($map_direction == 1)
    {
	my $full_read = reverse $R1_fields[9];
	$full_read =~ tr/ATGCatgc/TACGtacg/;
	$read_primer = substr($full_read, 0, $primer_length);
    }
    if(!exists($primer_hash{$read_primer}))
    {
    	next;
    }
    $mapped_reads_length_with_PCRp++;
    
    # Add scar to tables.
    my $bc_number = $bc_hash{$barcode};
    my $cigar = $R1_fields[5];
    my $location = $R1_fields[3];
    my $flag = $R1_fields[1];
    my $SNP = $R1_fields[12];
    my $sequence = $R1_fields[9];
    	
    print $scarfile_out "$cigar\t$barcode\t$UMI\t$location\t$flag\t$SNP\t$sequence\n";
    if(exists $CIGAR_hash{$cigar})
    {
    	# Count the number of times we have seen this CIGAR.
		$CIGAR_hash{$cigar}++;
    }else
    {
    	$CIGAR_hash{$cigar} = 1;
    }
    
    if(exists $scar_UMI_hash{$cigar."_".$sequence."_".$location}{$bc_number}{$UMI})
    {
		# Count the number of times we've seen this scar in this barcode with this UMI.
		$scar_UMI_hash{$cigar."_".$sequence."_".$location}{$bc_number}{$UMI}++;
    }else
    {
		$scar_UMI_hash{$cigar."_".$sequence."_".$location}{$bc_number}{$UMI} = 1;
		if(exists $scar_UMI_unique_hash{$cigar."_".$sequence."_".$location}{$bc_number})
		{
			# Count the number of times we have seen a unique UMI for this scar in this
			# barcode.
			$scar_UMI_unique_hash{$cigar."_".$sequence."_".$location}{$bc_number}++;
		}else
		{
			$scar_UMI_unique_hash{$cigar."_".$sequence."_".$location}{$bc_number} = 1;
		}
	}
	
	if(exists $scar_hash{$cigar."_".$sequence."_".$location}{$bc_number})
	{
		# Count the number of times we have seen this scar in this barcode.
		$scar_hash{$cigar."_".$sequence."_".$location}{$bc_number}++;
	}else
	{
		$scar_hash{$cigar."_".$sequence."_".$location}{$bc_number} = 1;
	}
}
close($samfile_in);
}

# Discard UMIs that do not have at least $min_UMI_count reads
my $scarfile_out_filtered_name = $out_prefix."_scars_filter.txt";
open (my $scarfile_out_filtered, ">".$scarfile_out_filtered_name);
print $scarfile_out_filtered "Barcode\tUMI\tCIGAR\tSequence\tLeft location\tReads\n";
my $reads_UMI_discarded = 0;
my $reads_UMI_kept = 0;
foreach my $scaratloc (keys %scar_UMI_hash)
{
	foreach my $barcode (keys %{$scar_UMI_hash{$scaratloc}})
	{
		foreach my $UMI (keys %{$scar_UMI_hash{$scaratloc}{$barcode}})
		{
			if($scar_UMI_hash{$scaratloc}{$barcode}{$UMI} < $min_UMI_count)
			{
				# Remove one UMI from the total UMI count of this scar.
				$scar_UMI_unique_hash{$scaratloc}{$barcode}--;
		
				# Change the number of reads for a given scar in a given barcode if we
				# discard UMIs.
				$scar_hash{$scaratloc}{$barcode} -= $scar_UMI_hash{$scaratloc}{$barcode}{$UMI};
		
				# Add to the total number of reads discarded.
				my $discarded_reads = $scar_UMI_hash{$scaratloc}{$barcode}{$UMI};
				$reads_UMI_discarded += $discarded_reads;
				# print "Discarded UMI $UMI for scar $scaratloc in cell $barcode: $discarded_reads discarded.\n";
				# print "$reads_UMI_discarded reads discarded in total.\n";
			}else
			{
				print $scarfile_out_filtered "$barcode\t$UMI\t";
				my @scaratloc_fields = split("_", $scaratloc);
				print $scarfile_out_filtered "$scaratloc_fields[0]\t$scaratloc_fields[1]\t$scaratloc_fields[2]\t";
				print $scarfile_out_filtered "$scar_UMI_hash{$scaratloc}{$barcode}{$UMI}\n";

				my $kept_reads = $scar_UMI_hash{$scaratloc}{$barcode}{$UMI};
				$reads_UMI_kept += $kept_reads;
				print $scarfile_out_filtered "Logging additional $kept_reads as kept, total reads kept $reads_UMI_kept.\n";
			}
		}
	}
}


# Print out CIGAR counts
# Calculation of reads kept/discarded is already done above
# my $reads_kept = -$reads_UMI_discarded;
my $readfile_out_name = $out_prefix."_cscars_coutc.txt";
open (my $readfile_out, ">".$readfile_out_name);
print $readfile_out "CIGAR\tSequence\tLocation\t";
my $scarfile_out_name = $out_prefix."_cscars_couta.txt";
open (my $scarfile_out, ">".$scarfile_out_name);
print $scarfile_out "CIGAR\tSequence\tLocation\t";
for (my $bc_number = 1; $bc_number <= $max_barcode; $bc_number++)
{
    print $readfile_out "$bc_number\t";
    print $scarfile_out "$bc_number\t";
}
print $readfile_out "\n";
print $scarfile_out "\n";
for my $scaratloc (keys %scar_hash)
{
    my @scaratloc_fields = split("_", $scaratloc);
    # REMOVE THE CONSTRAINTS ON CIGARs AND PERCENTAGE OF CIGARS.
	# if ($CIGAR_hash{$scaratloc_fields[0]} < $min_CIGAR_reads)
    # {
     	# print "CIGAR $scaratloc_fields[0] has $CIGAR_hash{$scaratloc_fields[0]} reads which is smaller than $min_CIGAR_reads.\n";
    # 	next;
    #}
    # print "CIGAR $scaratloc_fields[0] has $CIGAR_hash{$scaratloc_fields[0]} reads which is larger than $min_CIGAR_reads.\n";
        
    # my $threshold = ceil($min_CIGAR_fraction * $CIGAR_hash{$scaratloc_fields[0]});
    # my $seq_specific_sum = 0;
    # foreach my $barcode (keys %{$scar_hash{$scaratloc}})
    # {
    #	$seq_specific_sum += $scar_hash{$scaratloc}{$barcode};
    # }
	# if($seq_specific_sum < $threshold)
	# {
		# print "Scar $scaratloc has $seq_specific_sum which is smaller than $threshold.\n";
	#	next;
	# }
	# print "Scar $scaratloc has $seq_specific_sum which is larger than $threshold.\n";

	# Calculation of reads kept/discarded is already done above
	# $reads_kept += $seq_specific_sum;
    
    my $scar_reads_flag = 0;
    # print $readfile_out "Set scar reads flag to $scar_reads_flag\n";
    my $readfile_line = "$scaratloc_fields[0]\t$scaratloc_fields[1]\t$scaratloc_fields[2]\t";
    my $scarfile_line = "$scaratloc_fields[0]\t$scaratloc_fields[1]\t$scaratloc_fields[2]\t";
    # print $readfile_out "$scaratloc_fields[0]\t$scaratloc_fields[1]\t$scaratloc_fields[2]\t";
    # print $scarfile_out "$scaratloc_fields[0]\t$scaratloc_fields[1]\t$scaratloc_fields[2]\t";
    for (my $bc_number = 1; $bc_number <= $max_barcode; $bc_number++)
    {
		{
			no warnings 'uninitialized';
			# if(exists $scar_hash{$scaratloc}{$bc_number})
			if($scar_hash{$scaratloc}{$bc_number} > 0)
			{
	    		$scar_reads_flag++;
	    		# print $readfile_out "On barcode $bc_number, increased scar reads flag to $scar_reads_flag.\n";
	    		my $reads = $scar_hash{$scaratloc}{$bc_number};
	    		my $unique_reads = $scar_UMI_unique_hash{$scaratloc}{$bc_number};
	    		# print "Scar $scaratloc for barcode $bc_number has $reads reads, of which $unique_reads unique.\n";
	    		my $corrected_reads;
	    		if($unique_reads == $UMI_number)
	    		{
	    			$corrected_reads = "Inf";
	    		}else
	    		{
	    			$corrected_reads = log(1 - $unique_reads/$UMI_number)/log(1 - 1/$UMI_number);
	    		}
	    		$readfile_line = $readfile_line."$reads\t";
	    		# print $readfile_out "Found $reads reads, set print line to $readfile_line.\n";
	    		$scarfile_line = $scarfile_line."$corrected_reads\t";
	    		# print $readfile_out "$reads\t";
	    		# print $scarfile_out "$corrected_reads\t";
			}else
			{
	    		$readfile_line = $readfile_line."0\t";
	    		$scarfile_line = $scarfile_line."0\t";
	   			# print $readfile_out "0\t";
	   			# print $scarfile_out "0\t";
			}
    	}
    }
    if($scar_reads_flag > 0)
    {
    	# Add scar to the output files if we have found at least one read in one barcode.
    	# print $readfile_out "Scar reads flag is $scar_reads_flag, which is greater than 0.\n";
    	print $readfile_out $readfile_line."\n";
    	print $scarfile_out $scarfile_line."\n";
    }
}
close($readfile_out);
close($scarfile_out);

# Write statistics file
my $statfile_out_name = $out_prefix."_stats.txt";
open (my $stat_outfile, ">".$statfile_out_name);
print $stat_outfile "Reads: $reads\n";
my $reads_map_perc = 100 * $mapped_reads/$reads;
print $stat_outfile "Mapped to $gene: $mapped_reads ($reads_map_perc %)\n";
my $reads_map_barcode_perc = 100 * $mapped_reads_barcode/$reads;
print $stat_outfile "Mapped to $gene with correct barcode and no N in UMI: $mapped_reads_barcode ($reads_map_barcode_perc %)\n";
my $reads_map_length_perc = 100 * $mapped_reads_length/$reads;
print $stat_outfile "Mapped to $gene with right length: $mapped_reads_length ($reads_map_length_perc %)\n";
my $reads_map_length_PCRp_perc = 100 * $mapped_reads_length_with_PCRp/$reads;
print $stat_outfile "Mapped to $gene with right length and primer: $mapped_reads_length_with_PCRp ($reads_map_length_PCRp_perc %)\n";
my $reads_UMI_discarded_perc = 100 * $reads_UMI_discarded/$reads;
my $reads_UMI_kept_perc = 100 * $reads_UMI_kept/$reads;
# print $stat_outfile "Reads over thresholds ($min_CIGAR_reads per CIGAR/loc, $min_CIGAR_fraction of CIGAR/loc reads, at least $min_UMI_count reads per UMI): $reads_kept ($reads_kept_perc %)\n";
print $stat_outfile "Reads over thresholds (at least $min_UMI_count reads per UMI): $reads_UMI_kept ($reads_UMI_kept_perc %)\n";
print $stat_outfile "Reads discarded (less than $min_UMI_count reads per UMI): $reads_UMI_discarded ($reads_UMI_discarded_perc %)\n";
close($stat_outfile);
