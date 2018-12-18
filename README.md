# Trimming_HS_Metagenomics

This repository is made for hosting codes that can be used for trimming the heterogeneity spacers (HS) bases before the primer in a 16s rRNA metagenomics amplicon library. User should modify the Running Script for his/her own project, and that script with call the Perl code. The code searches for a provided primer string, and trims anything before it. The goal is to keep reads that have the premier attached to them, without any HS bases before the primer.

Within the Running Script, user should modify the following three parameters:

input_data_loc="location to FastQ.gz folder"
perl_code="location to Trimming_HS_Primers.pl scripts"
primer="primer sequence"

input_data_loc should have all the input *.fastq.gz files for trimming.

Then, by running the script, a perl code is called to search for the primer and everything before the primer is trimmed off. Reads that do not have primer attached also are discarded.


Input read:

HSHSHSHSHSHSPPPPPPPPPPPPPPPPPATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
-HS bases -- primer bases -- Read

In the example above, program will search for "PPPPPPPPPPPPPPPPP" in every give read, trims off all the "HSHSHSHSHSHS" bases, and the resulting read will be:

PPPPPPPPPPPPPPPPPATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
primer bases  -- Read

The length of the HS or the primer sequence will not affect the results, since the program will look for the exact match for input primer sequence.

Example command for only 1 input fasta file:

perl Trimming_HS_Primers.pl < input.fastq sample_ID "PPPPPPPPPPPPPPPPP"

Example command for running the Running Scrip, which will point to the location of multiple input fasta files:

sh Running_Script.sh
