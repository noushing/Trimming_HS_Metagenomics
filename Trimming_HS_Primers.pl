#Created by Noushin Ghaffari(nghaffari@tamu.edu) on 11/19/2018.
#Copyright 2018 __Texas A&M University__AgriLife Genomics & Bioinformatics__ All rights reserved.
#This code searches for a provided primer string, and trims anything before it. The goal is to trim the heterogeneity spacers (HS) bases before the primer in a 16s rRNA metagenomics amplicon libarary.


#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;

use List::Util;
use POSIX;

my ($read_identifier, $sequence, $qual_scor_str);

my @Temp_arr = undef;

my @arr_FASTQ_R1 = undef;
my @arr_log = undef;


my $output_R1_file_name = ($ARGV[0] . "_TrimmedHS_R1_FASTQ.fastq");

unlink $output_R1_file_name;


my $total_reads_to_be_written_in_File = 4;

my $loc_primer=0;
my $read_number=1;

print strftime "%Y-%m-%d %H:%M:%S\n", localtime();

while (($read_identifier, $sequence, $qual_scor_str) = Read_Sep_FASTQ(\*STDIN, \@Temp_arr))
{
    my $primer = $ARGV[1];
    
    #print("Read Number: ",$read_number."\n");
    
  
    #Finding the postion of first base of the primer
    $loc_primer=index $sequence, $primer;
    
    
    #This if checks to ensure reads that have primers attached to them to be saved in the output, if not, we only keep info of those reads in the log file, but not in the output fastq file.
    if ($loc_primer >= 0)
    {
        if (defined($arr_FASTQ_R1[0]))
        {
            #Making read1
            push(@arr_FASTQ_R1,"@".$read_identifier."\n");
            push(@arr_FASTQ_R1,substr($sequence,$loc_primer,length($sequence))."\n");
            push(@arr_FASTQ_R1,"+\n");
            push(@arr_FASTQ_R1,substr($qual_scor_str,$loc_primer,length($sequence))."\n");
            
            
        }
        else
        {
            #Making read1
            $arr_FASTQ_R1[0] = "@".$read_identifier."\n";
            push(@arr_FASTQ_R1,substr($sequence,$loc_primer,length($sequence))."\n");
            push(@arr_FASTQ_R1,"+\n");
            push(@arr_FASTQ_R1,substr($qual_scor_str,$loc_primer,length($sequence))."\n");
            
        }
        
        
        if (defined($arr_log[0]))
        {
            #Saving trimmed bases in a log file
            push(@arr_log,$loc_primer." , ".substr($sequence,0,$loc_primer)."\n");
            
        }
        else
        {
            #Saving trimmed bases in a log file
            $arr_log[0] = $loc_primer." , ".substr($sequence,0,$loc_primer)."\n";
        }
        
        
        if(scalar(@arr_FASTQ_R1) == $total_reads_to_be_written_in_File)
        {
            #writing the fastq file
            Writing_FASTQ_Files();
            @arr_FASTQ_R1 = undef;
            #print("\nWriting to the output files ...\n");
            
        }
        #writing the log file
        Writing_log_File();
        @arr_log = undef();
    }
    #This else will save the info of reads with no primer in the log file.
    else
    {
        if (defined($arr_log[0]))
        {
            #Saving trimmed bases in a log file
            push(@arr_log,$loc_primer." , ".substr($sequence,0,$loc_primer)."\n");
            
        }
        else
        {
            #Saving trimmed bases in a log file
            $arr_log[0] = $loc_primer." , ".substr($sequence,0,$loc_primer)."\n";
        }
        
    }
    
    
    
    $read_number++;
}


sub Writing_FASTQ_Files
{
    
    #Saving Read 1
    open FILER1, ">>", $output_R1_file_name or die $!;
    print FILER1 @arr_FASTQ_R1;
    close FILER1;
    
}

sub Writing_log_File
{
    
    open FILERlog, ">>", $ARGV[0]."_log.txt" or die $!;
    print FILERlog @arr_log;
    close FILERlog;
    
}


#Modified version of: https://github.com/lh3/readfq/blob/master/readfq.pl
sub Read_Sep_FASTQ {
    
    my ($input_file, $Temp_arr) = @_;
    
    
    #NG: To fix the error of newer Perl versions (May 2012) #@$Temp_arr = [undef, 0] if (!defined(@$Temp_arr));
    return if ($Temp_arr->[1]);
    if (!defined($Temp_arr->[0])) {
        while (<$input_file>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $Temp_arr->[0] = $_;
                last;
            }
        }
        if (!defined($Temp_arr->[0])) {
            $Temp_arr->[1] = 1;
            return;
        }
    }
    #NG: Changed to the following to include everything in the identifier, including the PE 1 or 2 info after the whitespaces.
    my $read_identifier = substr($_,1,length($_));
    my $sequence = '';
    my $c;
    $Temp_arr->[0] = undef;
    while (<$input_file>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $sequence .= $_;
    }
    $Temp_arr->[0] = $_;
    $Temp_arr->[1] = 1 if (!defined($Temp_arr->[0]));
    return ($read_identifier, $sequence) if ($c ne '+');
    my $qual_scor_str = '';
    while (<$input_file>) {
        chomp;
        $qual_scor_str .= $_;
        if (length($qual_scor_str) >= length($sequence)) {
            $Temp_arr->[0] = undef;
            return ($read_identifier, $sequence, $qual_scor_str);
        }
    }
    $Temp_arr->[1] = 1;
    return ($read_identifier, $sequence);
}



