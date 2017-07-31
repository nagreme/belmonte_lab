#!/usr/bin/perl

#-----------------------------Description------------------------------
# AUTHOR: Nadège Pulgar-Vidal
# DATE: July 21th 2016
#
# PROGRAM: SeqAnalysisPipe (SAP bash shorcut)
#
# DESCRIPTION: Gets reads files (.fastq), a genome annotation file (.gff/.gff3), 
#	       and a fasta file (.fa) from the user and then builds
#	       and runs the commands for trimmomatic, bowtie2-build,
#	       tophat2, cuffquant, cuffnorm, and cuffdiff. Prompting
#	       system is interactive and fairly robust. Running cuffdiff
#	       is optional and there is an option to only build commands
#	       and not run them. Either way a copy of the built commands
#	       is saved to a text file: "commands.txt".
#----------------------------------------------------------------------


#might require installation of additional perl modules (the includes ones) on other computers


#-------------------------------Includes-------------------------------

use strict;
use warnings; 

use Readonly;
use String::Util qw(trim);
use Scalar::Util qw(looks_like_number);

#------------------------------Constants-------------------------------

#Default settings
#Trimmomatic
use constant TRIMMOMATIC_JAR => "/home/mark/bioinformaticsv2/Programs/Trimmomatic-0.33/trimmomatic-0.33.jar";
use constant MIN_THREADS => 1;
use constant MAX_THREADS => 20;
use constant SE => "SE";
use constant PE => "PE";
use constant SE_ADAPTERS => "/home/mark/bioinformaticsv2/Programs/Trimmomatic-0.33/adapters/TruSeq3-SE.fa";
use constant PE_ADAPTERS => "/home/mark/bioinformaticsv2/Programs/Trimmomatic-0.33/adapters/TruSeq3-PE-2.fa";
use constant ADAPTER_SETTINGS => ":2:30:10";
use constant HEADCROP => 9;
use constant LEADING => 30;
use constant TRAILING => 30;
use constant SLIDING_WINDOW => "4:30";
use constant MINLEN => 50;
use constant AVGQUAL => 30;
#Tophat2
use constant TOPHAT_CONSTANT_OPTIONS => "--b2-very-sensitive --no-coverage-search"; 

#Output directories
use constant TRIM_OUT_DIR => "./trimmomatic_output";
use constant TOPHAT_OUT_DIR => "./tophat2_output";
use constant CQ_OUT_DIR => "./cuffquant_output";
use constant CN_OUT_DIR => "./cuffnorm_output";
use constant CD_OUT_DIR => "./cuffdiff_output";

#Output extensions
use constant TRIM_OUT_EXT => ".fq";
use constant TRIM_PAIRED_EXT => "_paired";
use constant TRIM_UNPAIRED_EXT => "_unpaired";
use constant TOPHAT_OUT_DIR_EXT => "_aligned";

#Input filenames/extensions
use constant TRIM_INPUT_EXT => ".fastq";
use constant CQ_INPUT_FILE => "accepted_hits.bam"; #input filename for cuffquant
use constant CND_INPUT_FILE => "abundances.cxb"; #cuffnorm and cuffdiff use the same input files

#Other
use constant CMDS_FILENAME => "commands.txt";
use constant THREADS_OVERRIDE => "CTRL";
use constant ALL_ALIGN_SUMS_FILENAME => "all_align_summaries.txt";
use constant ALIGN_SUM_FILENAME => "align_summary.txt";

#------------------------------Variables-------------------------------
#choices
my $use_cuffdiff; #user response to prompt to use cuffdiff
my $run_cmds; #user response to prompt asking if the cmds should be run or just built

#Data
my @reads_files_ext; #list of reads file names
my @reads_files; #list of reads file names without extensions
my @reads_files_pairs; #list of pairs of reads files for PE trimmomatic without extensions
my $genome_file; #name of the genome annotation file (.gff/.gtf/.gff3)
my $fasta_file; #name of the fasta file (.fa/.fna) 
#if you change the trimmomatic jar or adapaters file change the name in the trim_settings function prompt (in the "these are the seetings we will use" display)

#Output Directories
my @combined_dirs; #list of names of the combined alignment dirs (PE, tophat)

#Settings (defaults)
#General
my $threads = 15;

#Trimmomatic
my $end_type = SE;
my $adapters_file = SE_ADAPTERS;
my $adapters_settings = ADAPTER_SETTINGS;
my $headcrop = HEADCROP;
my $leading = LEADING;
my $trailing = TRAILING;
my $sliding_window = SLIDING_WINDOW;
my $minlen = MINLEN;
my $avgqual = AVGQUAL;

#Tophat2
my $genome_index_base;

#CuffNorm
my $cn_labels; # string of comma separated cuffquant labels
my $cn_files; #string of comma and space separated cuffquant input files

#CuffDiff
my $cd_labels; # string of comma separated cuffdiff labels
my $cd_files; #string of comma and space separated cuffdiff input files


#---------------------------------Main---------------------------------

#Clear
system("clear");

#Start
print_banner(); #done

#Get everything we need from the user
get_genome_and_fasta(); #done
get_reads_files(); #done
#get_threads(); #done
#get_trim_setting(); #done
get_cn_labels_and_files(); #done


#make cuffdiff optional
print "\n\n###################################
Do you want to run Cuffdiff after Cuffnorm? (y/n)
(Press enter for no)\n";
trim($use_cuffdiff = <>);

if ($use_cuffdiff =~ m/^[Yy]/)
{
get_cd_labels_and_files(); #done
}#if (use cuffdiff)


#Open the text file that we'll print all the used cmds in
open(CMDS_FILE,">".CMDS_FILENAME) or die "Couldn't open ".CMDS_FILENAME.": $!\n";

#Print a little header to the file with the input files provided by the user
print_cmds_file_header();


#make running the commands optional: offer the option to build only
#prompt
print "\n\nDo you want to run all the commands or just build them? (run/build) (They will be saved to an output text file called 'commands.txt')\n(Press enter for run)\n";

trim($run_cmds = <>);

#if the user only wants to build commands
if ($run_cmds =~ m/^[Bb]/)
{
print "\nYou have chosen to only build the commands to use. 
Please check the 'commands.txt' file when the pipeline finishes running.\n";

#sleep for 6 seconds so the user can read the message before it is pushed up by other prints
sleep(6);
}#if

#Run everything!
run_trim(); #initial test done (SE,PE)
run_bowtie(); #initial test done (SE,PE)
run_tophat(); #initial test done (SE,PE)

#concatenate the alignment summary files if we ran tophat
unless ($run_cmds =~ m/^[Bb]/)
{
concatenate_align_summaries();
}#unless (only build)

run_cuffquant(); #intial test done (SE,PE)
run_cuffnorm(); #intial test done (SE,PE)


#Note: run_bowtie() checks if the index is already built to avoid 
#a lot of potentially extra work. It does this by checking the 
#existence of two of the files the index should contain (&&) in the 
#current directory. 


if ($use_cuffdiff =~ m/^[Yy]/)
{
run_cuffdiff(); #intial test done (SE,PE)
}#if (use cuffdiff)


#RScripts? Non catégorique.


#Print the end time of the cmds file
print_cmds_files_end();

#close the cmds_file because we're good programmers
close(CMDS_FILE) or warn "Couldn't close file properly\n";


#let the user know everything is over
print_end_message();



#-------------------------Function definitions-------------------------

#----------------------------------------------------------
sub print_banner
{
print "\n*************************************************************
Welcome to the Sequencing Analysis Pipeline!\n";

print "\nThe pipeline puts your reads through trimmomatic, tophat2 (bowtie2-build), cuffquant, cuffnorm, and cuffdiff subsequently. 

The required files are: reads files (.fastq, short intuitive names are recommended), a genome annotation file (.gff), and a fasta file (.fa). 

Please note that the genome annotation file must have the same name as the fasta file before the extension and that all the required files must be in the current directory.\n
*************************************************************\n";
}#print_banner


#----------------------------------------------------------
sub print_cmds_file_header
{
my $date_string = localtime(); #hold the current date and time

#Print the date and time
print CMDS_FILE "Time and date created: $date_string\n\n";

#print a header
print CMDS_FILE "This files contains the list of all commands run by the Sequencing Analysis Pipeline.\n";

#and the user files used to build the commands
print CMDS_FILE "\nThese are the input files that were used:
Genome annotation file: ".trim($genome_file)."
Fasta file: ".trim($fasta_file)."
Reads file(s): @reads_files_ext


************
* Commands *
************\n\n";

}#print_cmds_file_header



#----------------------------------------------------------
sub get_genome_and_fasta
{
my $valid_input = 0; #valid input flag
my $gen_file; #local genome annotation file 
my $fa_file; #local fasta file 
my $base_name; #local genome-index base 
my $user_resp; #user response to final check


print "\n\n###################################
Genome and fasta\n";

#one loop until valid input
do
{
#Prompt for genome (scalar)
print "\nEnter the name of the genome annotation file (.gff):\n";

#read genome annotation file name
trim($gen_file = <>);

#check existance of genome annotation file (continue only if valid)
if (-f trim($gen_file))
{
#Prompt for fasta (scalar)
print "\n\nEnter the name of the fasta file (.fa):\nPlease note that it must have the same name as the genome annotation file before the extension.\n";

#read fasta file name
trim($fa_file = <>);

#check existance of fasta file (continue only if valid)
if (-f trim($fa_file))
{
#check name match before extension
#extract base name from one (without extension)
$base_name = substr($gen_file, 0, rindex($gen_file,".")); #extract chars up to the last .

#look for "base_name." in the other (using index or regex)
if (index($fa_file,"$base_name.") >= 0)
{
#get the user to validate (final check)
#Prompt
print "\n\nThese are the files you have chosen:
Genome annotation file: ".trim($gen_file)."
Fasta file: $fa_file

Is this information correct? (y/n)
(Enter means yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle loop flag
$valid_input = 1;

}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}#if (name match)
else
{
printf "\nThe names of the genome and fasta files have to match before the extension. Try again.\n";
}#else (name match)

}#if (fa_file exists)
else
{
printf "\nThis file does not exist, try again.\n";
}#else (fa_file exists)

}#if (gen_file exists)
else
{
printf "\nThis file does not exist, try again.\n";
}#else (gen_file exists)

}while(not($valid_input));

#save genome and fasta files
$genome_file = $gen_file;
$fasta_file = $fa_file;

#save genome index base
$genome_index_base = $base_name;

}#get_genome_and _fasta


#----------------------------------------------------------
sub get_reads_files
{
my $valid_input = 0; #flag
my @file_list; #hold names of reads files for checking (function)
my $user_input; #hold the string of input from the user
my $check_return; #hold a message return from the check function
my $user_resp; #user response to final check
my $i; #for loop iterator

print "\n\n###################################
Reads files

Note: If you plan on using paired end (PE) settings for trimmomatic make sure to include all the reads files you will need.\n";

#loop until valid input
do
{
#prompt for reads file list, space separated (scalar, split to array)
print "\nEnter the reads file(s) (.fastq) in one line, *space seperated*:
(You must enter at least one file)\n";

#read in as string
trim($user_input = <>);

#split to array
@file_list = split(" ",$user_input);

#if there was input
if (not(scalar @file_list == 0))
{

#check existance of each (function)
$check_return = check_file_list(@file_list);

#if everything okay
if (not($check_return)) #true if check_return is blank
{
#get the user to validate (final check)
#Prompt
print "\n\nThese are the file(s) you entered:
@file_list

Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle loop flag
$valid_input = 1;
}#if (user check)
else
{
print "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}#if (file list check)
else
{
printf "\n$check_return Try again.\n";
}#else (file list check)

}#if (no input)
else
{
print "\nYou need to enter at least one input file.\n";
}#else (no input)

}while(not($valid_input));

#save array of file names with extensions
@reads_files_ext = @file_list;

#save array of file names without extensions
for ($i = 0; $i < scalar @file_list; $i++) 
{
$reads_files[$i] = substr($file_list[$i], 0, rindex($file_list[$i], ".")); #extract chars up to the extension
}#foreach

}#get_reads_files 


#----------------------------------------------------------
sub check_file_list #doesn't currently check for duplicate entries (deal with it)(for now) 
{
my @list = @_; 
my $messg = ""; #empty message means succes (all files exists)
my $item; #loop variable

#check that each file in the list exists
foreach $item (@list)
{
#if a file does not exist
if (not(-f $item))
{
$messg = "$item does not exist."; #save the error
last; #break out of the loop (sorry)
}#if
}#foreach

return $messg; #blank if eveything okay, else conatins error message
}#check_file_list


#----------------------------------------------------------
sub get_threads
{
my $user_check = 0; #flag
my $user_threads; #holds user input for number of threads
my $user_resp; #user response to final check

print "\n\n###################################
Threads\n";

#do while loop for user check
do
{
#prompt for threads (if nothing use default)
print "\nEnter the number of threads you would like to use (integers ",MIN_THREADS,"-",MAX_THREADS," inclusive):\n";

#read number of threads from user
trim($user_threads = <>);

#if it is a number
if (looks_like_number(trim($user_threads)))
{
#error check and adjust if needed to closest reasonable value
if ($user_threads < MIN_THREADS)
{
print "This number is too small so we will use ",MIN_THREADS," thread.\n";
$threads = MIN_THREADS;
}#if (too small)

elsif ($user_threads > MAX_THREADS)
{
print "This number is too big so we will use ",MAX_THREADS," threads.\n";
$threads = MAX_THREADS;
}#elsif (too big)

else
{
$threads = int(trim($user_threads));
}#else (valid number of threads)

}#if (it's a number)

#if it's not a number check for the override code
elsif (trim($user_threads) eq THREADS_OVERRIDE)
{
#If triggered take whatever the user enters next as the threads value
print "\nOverride triggered. Threads will be set to whatever you enter next:\n";

$user_threads = <>; 
$threads = trim($user_threads);
}#elsif (not a number, override)

#if the user wants max threads
elsif (trim($user_threads) =~ m/[Mm][Aa][Xx]/)
{
$threads = MAX_THREADS;
}#elsif (max threads selected)

#if the user wants min threads
elsif (trim($user_threads) =~ m/[Mm][Ii][Nn]/)
{
$threads = MIN_THREADS;
}#elsif (min threads selected)

else #not a number, nor a special code
{
print "\nThat's not a number so the previous value will be used.\n";
}#else (non numeric input)

#get the user to validate (final check)
#Prompt
print "\n\nThis is the number of threads that will be used:
$threads

Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle loop flag
$user_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}while(not($user_check));

}#get_threads


#----------------------------------------------------------
sub get_trim_setting
{
my $user_resp; #hold user response to prompts (y/n)
my $user_check = 0; #flag
my $user_input; #hold whatever the user inputs

print "\n\n###################################
Trimmomatic

(More info: http://www.usadellab.org/cms/?page=trimmomatic)\n";


$adapters_file =~ m/^.*\//; #put the filename without path into $'

#show default trimmomatic settings
print "\nThese are the settings trimmomatic will use:
(non customizable)
trimmomatic-0.33.jar
phred33

(customizable)
$end_type
$' (ILLUMINACLIP:adapaters file)(linked to end type)
ILLUMINACLIP:adapters$adapters_settings (:seed_mismatches:palindrome_clip_threshold:simple_clip_threshold)
HEADCROP:$headcrop 
LEADING:$leading 
TRAILING:$trailing 
SLIDINGWINDOW:$sliding_window (window_size:required_quality)
MINLEN:$minlen
AVGQUAL:$avgqual\n";


#prompt for: use custom trim settings? 
print "\nDo you want to use these settings? (y/n)
(Press enter for yes, if you select no you will be prompted for custom settings)\n";

#read response
trim($user_resp = <>);

#no means changes setting everything else means use default
if ($user_resp =~ m/^[Nn]/)
{
#if yes: loop (do-while) with check at end (display chosen settings, check for correction)
do
{
#prompt for each setting sequentially

print "\nNote: If you do not want to change a setting you can leave it blank and press enter to move on to the next setting.\n";


get_end_type();
get_adapter_settings();
get_headcrop();
get_leading();
get_trailing();
get_sliding_window();
get_minlen();
get_avgqual();


$adapters_file =~ m/^.*\//; #put the filename without path into $'

#get the user to validate (final check)
#Prompt
print "\n\nThese are the trim settings that will be used:
(non customizable)
trimmomatic-0.33.jar
phred33

(customizable)
$end_type
$' (ILLUMINACLIP:adapaters file)(linked to end type)
ILLUMINACLIP:adapters$adapters_settings (:seed_mismathces:palindrome_clip_threshold:simple_clip_threshold)
HEADCROP:$headcrop 
LEADING:$leading 
TRAILING:$trailing 
SLIDINGWINDOW:$sliding_window (window_size:required_quality)
MINLEN:$minlen
AVGQUAL:$avgqual

Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle loop flag
$user_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}while(not($user_check));

#If the user selected PE settings prompt them for file pairs
if ($end_type =~ m/PE/)
{
get_PE_pairs();
}#if (PE)

}#if (change settings)


}#get_trim_setting


#all the gets for the trim settings
#----------------------------------------------------------
sub get_end_type
{
my $user_input;

#Prompt for end type
print "\nEnter the end type you want. Use SE for single end and PE for paired end (currently $end_type)
(Note: If you select PE you will be prompted to enter the pairs of read files)\n";

#read end type from user
trim($user_input = <>);

#check end type validity
check_end_type($user_input);
}#get_end_type


#----------------------------------------------------------
sub get_adapter_settings
{
my $user_input;

#Prompt for adapater settings
print "\nEnter the adapter settings you want (please include the first colon and follow the :#:#:# format) (currently $adapters_settings)\n";

#read adapter setting from the user
trim($user_input = <>);

#check adapater settings validity
check_adapter_settings($user_input);
}#get_adapter_settings


#----------------------------------------------------------
sub get_headcrop
{
my $user_input;

#Prompt for headcrop
print "\nEnter the headcrop value you want (>= 0) (currently $headcrop):\n";

#read headcrop from user
trim($user_input = <>);

#check headcrop settings validity
check_headcrop_settings($user_input);
}#get_headcrop


#----------------------------------------------------------
sub get_leading
{
my $user_input;

#Prompt for leading
print "\nEnter the leading value you want (>= 0) (currently $leading):\n";

#read leading from user
trim($user_input = <>);

#check leading settings validity
check_leading_settings($user_input);
}#get_headcrop


#----------------------------------------------------------
sub get_trailing
{
my $user_input;

#Prompt for trailing
print "\nEnter the trailing value you want (>= 0) (currently $trailing):\n";

#read trailing from user
trim($user_input = <>);

#check trailing settings validity
check_trailing_settings($user_input);
}#get_trailing


#----------------------------------------------------------
sub get_sliding_window
{
my $user_input;

#Prompt for sliding window
print "\nEnter the sliding window value you want (please follow the format #:#) (currently $sliding_window):\n";

#read sliding window from user
trim($user_input = <>);

#check sliding window settings validity
check_sliding_window_settings($user_input);
}#get_sliding_window


#----------------------------------------------------------
sub get_minlen
{
my $user_input;

#Prompt for minlen
print "\nEnter the minimum length value you want (>= 1) (currently $minlen):\n";

#read minlen from user
trim($user_input = <>);

#check minlen settings validity
check_minlen_settings($user_input);
}#get_minlen


#----------------------------------------------------------
sub get_avgqual
{
my $user_input;

#Prompt for avgqual
print "\nEnter the average quality value you want (>= 0) (currently $avgqual):\n";

#read trailing from user
trim($user_input = <>);

#check avgqual settings validity
check_avgqual_settings($user_input);
}#get_avgqual


#----------------------------------------------------------
sub get_PE_pairs
{
my $user_check = 0; #flag
my $valid_pair; #flag
my $user_pair; #holds a pair of filenames inputed by user
my $user_resp; #user response to final check
my @pair_files; #holds a pair of reads split into an array
my @curr_pairs; #holds an array of current reads pairs
my $check_result; #return message or error from the check
my $i; #loop iterator

print "\n\nYou will be prompted for the pairs of reads files to use for paired end trimming.\nThe number of pairs to be entered is the number of reads files divided by 2.\n";

#loop until user check 
do
{

#for each pair (for loop, num_iterations = num_reads / 2)
for ($i = 0; $i < int(scalar @reads_files/2); $i++)
{
$valid_pair = 0; #reset flag

#loop until valid pair input
do
{
#prompt
print "\nEnter the reads files for pair #$i without extensions, on one line, sperated by a space, choosing from these filenames:
@reads_files

Reminder of current pairs:\n";
print_array_by_line(@curr_pairs);
print "\n";

#read input 
trim($user_pair = <>);

#split to arr for check
@pair_files = split(" ",$user_pair);

#If there are two files 
if (scalar @pair_files == 2)
{
#the pair contains 2 different files
if ($pair_files[0] ne $pair_files[1])
{
#check that the files are in @reads_files
$check_result = check_files_in_list(@pair_files);

#and they both exist in our list
if (not($check_result)) #true if check_result is blank (no error message)
{
#get the user to validate (final check)
#Prompt
print "\n\nThese are the files you entered for pair #$i:
@pair_files
Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$valid_pair = 1;
}#if (user check)
else
{
printf "\n\nYou said the information was incorrect so try again.\n";
}#else (user check)

}#if (files exist)
else
{
print "\n$check_result Try again.\n";
}#else (error)

}#if (same file entered twice for one pair)
else
{
print "\n\nA pair must contain 2 different files, try again.\n";
}#else (duplicate in pair)

}#if (not 2 input files for pair)
else
{
print "\n\nA pair must contain 2 files, try again.\n";
}#else (not 2 input files for pair)

}while(not($valid_pair));

#if good add to arr curr_pairs
push (@curr_pairs, trim($user_pair));
}#for


#final user check, display curr data
print "\n\nThese are the ".int(scalar @reads_files/2)." pair(s) you've entered:\n";
print_array_by_line(@curr_pairs);

print "\nIs this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle loop flag
$user_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}while(not($user_check));

#save reads_files_pairs
@reads_files_pairs = @curr_pairs;

}#get_PE_pairs



#----------------------------------------------------------
sub print_array_by_line
{
my $item; #hold an item in @_

foreach $item (@_)
{
print "$item\n";
}#foreach

}#print_array_by_line



#check validity (if empty or invalid use defaults and notify user)
#(valid entries override current default, invalid does nothing)
#(save valid entries as you go)
#----------------------------------------------------------
sub check_end_type
{
#if not blank
if (trim($_[0]))
{
#if it matches PE
if ($_[0] =~ m/[Pp]+[Ee]*/)
{
$end_type = PE;
$adapters_file = PE_ADAPTERS;
}#if (PE)

elsif ($_[0] =~ m/[Ss]+[Ee]*/)
{
$end_type = SE;
$adapters_file = SE_ADAPTERS;
}#elsif (SE)

else
{
print "\nInvalid input. End type not updated.\n";
}#else (invalid format)

}#if (not blank)

}#check_end_type


#----------------------------------------------------------
sub check_adapter_settings
{
#if not blank
if (trim($_[0]))
{
#if it follow the right format update the value
if ($_[0] =~ m/:\d+:\d+:\d+/)
{
$adapters_settings = trim($_[0]);
}#if (format)
else
{
print "\nWrong format. Adapter settings not updated.\n";
}#else (format)
}#if (not blank)

} #check_adapter_settings


#----------------------------------------------------------
sub check_headcrop_settings
{
#if not blank
if (trim($_[0]))
{
#if it is a number
if (looks_like_number($_[0]))
{
#if it is within reasonable values
if ($_[0] >= 0)
{
$headcrop = trim($_[0]);
}#if (value)
else
{
print "\nHeadcrop value must be greater than or equal to 0. Headcrop settings not updated.\n";
}#else (value)
}#if (number)
else
{
print "\nInput must be a number. Headcrop settings were not updated.\n";
}#(number)
}#if (not blank)

} #check_headcrop_settings


#----------------------------------------------------------
sub check_leading_settings
{
#if not blank
if (trim($_[0]))
{
#if it is a number
if (looks_like_number($_[0]))
{
#if it is within reasonable values
if ($_[0] >= 0)
{
$leading = trim($_[0]);
}#if (value)
else
{
print "\nLeading value must be greater than or equal to 0. Leading settings were not updated.\n";
}#else (value)
}#if (number)
else
{
print "\nInput must be a number. Leading settings were not updated.\n";
}#(number)
}#if (not blank)

} #check_leading_settings


#----------------------------------------------------------
sub check_trailing_settings
{
#if not blank
if (trim($_[0]))
{
#if it is a number
if (looks_like_number($_[0]))
{
#if it is within reasonable values
if ($_[0] >= 0)
{
$trailing = trim($_[0]);
}#if (value)
else
{
print "\nTrailing value must be greater than or equal to 0. Trailing settings were not updated.\n";
}#else (value)
}#if (number)
else
{
print "\nInput must be a number. Trailing settings were not updated.\n";
}#(number)
}#if (not blank)

} #check_trailing_settings 


#----------------------------------------------------------
sub check_sliding_window_settings
{
#if not blank
if (trim($_[0]))
{
#if it follow the right format update the value
if ($_[0] =~ m/\d+:\d+/)
{
$sliding_window = trim($_[0]);
}#if (format)
else
{
print "\nWrong format. Sliding window settings were not updated.\n";
}#else (format)
}#if (not blank)
} #check_sliding_window_settings


#----------------------------------------------------------
sub check_minlen_settings
{
#if not blank
if (trim($_[0]))
{
#if it is a number
if (looks_like_number($_[0]))
{
#if it is within reasonable values
if ($_[0] >= 0)
{
$minlen = trim($_[0]);
}#if (value)
else
{
print "\nMinimum length value must be greater than or equal to 1. Minimum length settings were not updated.\n";
}#else (value)
}#if (number)
else
{
print "\nInput must be a number. Minimum length settings were not updated.\n";
}#(number)
}#if (not blank)

} #check_minlen_settings


#----------------------------------------------------------
sub check_avgqual_settings
{
#if not blank
if (trim($_[0]))
{
#if it is a number
if (looks_like_number($_[0]))
{
#if it is within reasonable values
if ($_[0] >= 0)
{
$avgqual = trim($_[0]);
}#if (value)
else
{
print "\nAverage quality value must be greater than or equal to 0. Average quality settings not updated.\n";
}#else (value)
}#if (number)
else
{
print "\nInput must be a number. Average quality settings not updated.\n";
}#(number)
}#if (not blank)

} #check_avgqual_settings



#----------------------------------------------------------
sub get_cn_labels_and_files
{
my $user_resp; #hold user response to prompts (y/n)
my $user_check = 0; #flag for final check
my $label_check = 0; #flag for checking labels entry
my @labels; #holds the labels as an array so we can assign files to them
my $labels_str; #holds current labels in comma separated string
my $files_str = ""; #holds current files in comm/space seperated string
my $curr_label; #foreach loop variable, holds one label


print "\n\n###################################
Cuffnorm

(More info: http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/index.html)\n";

#loop until all information (labesl and files) is correct
do
{

#loop until we get the labels we want
do
{
#Prompt for cn_labels
print "\nEnter the label names to use for Cuffnorm in a comma separated list:
(These will be used to group your bioreps. You must enter at least one label)\n";

#read labels from user as string
trim($labels_str = <>);

#split to array
@labels = split(",",$labels_str);

if (trim($labels_str))
{

#check with user that they want to use these labels
print "\n\nThese are the label names you entered for Cuffnorm:
@labels
Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$label_check = 1;
}#if (user check)
else
{
print "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}#if (at least one label)
else
{
print "\nYou must enter at least one label. Try again.\n";
}#else (at least one label)

}while(not($label_check));


#for each label, loop until valid: prompt for files in label (print list of included files without ext.)
foreach $curr_label (@labels)
{

$files_str .= get_label_files($curr_label)." ";
#once label info is good tack on to cn_files string in proper format (comma separated, space at end)

}#foreach (label)


#ask: corrections? (display labels and files: if yes re-do, else return)
print "\n\nThese are the label names and files you entered for Cuffnorm:
(Note: If you selected PE settings you will see a '/' between the files in a pair)
labels: ".trim($labels_str)."
files: ".trim($files_str)."

Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$user_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
#reset file_str
$files_str ="";
}#else (user check)


}while(not($user_check));

#save cn_labels and cn_files
$cn_labels = $labels_str;
$cn_files = $files_str;

}#get_cn_labels_and_files


#----------------------------------------------------------
sub get_label_files #needs a label as only parameter
{
my $label_check = 0; #flag
my $check_result; #holds return message from check function
my $user_input; #hold whatever the user inputs
my $user_resp; #user resonse to final check
my @files; #holds split files

do
{
#Prompt

#Single End
if ($end_type =~ m/SE/)
{
print "\n\nWhich files belong to the label \"".trim($_[0])."\"?\n";
print "Enter these as a comma separated list without file extensions choosing from these:
@reads_files\n";
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
print "\n\nWhich pairs belong to the label \"".trim($_[0])."\"?\n";
print "Enter these as a comma separated list without file extensions choosing from these pairs
(Note: please write the pair names exactly as printed, including the space between the pair):\n";
print_array_by_line(@reads_files_pairs);
print "\n";
}#elsif (PE)


#read files as string
trim($user_input = <>);

#split into array
@files = split(",",trim($user_input));


#Check input
#Single End
if ($end_type =~ m/SE/)
{
#split (,) to temp array for check and get result
$check_result = check_files_in_list(@files);
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
$check_result = check_pairs_in_list(@files);
}#elsif (PE)


#if everything was okay
if (not($check_result)) #true if check_result is blank (no error message)
{
#get the user to validate (final check)
#Prompt
#Single End
if ($end_type =~ m/SE/)
{
print "\n\nThese are the files you entered for label \"".trim($_[0])."\":\n";
print "@files\n";
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
print "\n\nThese are the pairs you entered for label \"".trim($_[0])."\":\n";
print_array_by_line(@files);
}#if (PE)


print "\nIs this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$label_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}#if (no error)
else
{
print "\n$check_result Try again.\n";
}#else (error)

}while(not($label_check));


#adjustment for PE to be able to split the input_files string more precisely
#add a slash between the files in a pair instead of a space
if ($end_type =~ m/PE/)
{
$user_input =~ s/ /\//g;
}#if (PE, file list adjust)

#return the list of files comma separated
return trim($user_input);
}#get_label_files


#----------------------------------------------------------
sub check_pairs_in_list
{
my $pair; #holds one filename from @_
my $found; #flag
my $i; #loop counter
my $messg = "";

#check that there are filenames to check
if (not(scalar @_ == 0))
{
#check that each file (item in @_) is in the reads_files (linear search)
foreach $pair (@_)
{
$found = 0; #reset flag
$i = 0; #reset loop counter

#while we haven't found it in and we haven't fallen off the array
while (not($found) && $i < scalar @reads_files_pairs)
{
#if it matches the current filename
if (trim($pair) eq trim($reads_files_pairs[$i]))
{
#toggle flag
$found = 1; 
}#if (match)

$i++; #increment loop counter
}#while (search)

#if not found, save error message and break out of the loop
if (not($found))
{
$messg = "The pair \"".trim($pair)."\" was not found in list of reads files pairs.";
last; #(so sorry)
}#if (not found)

}#foreach (file)
}#if (files to check)
else
{
$messg = "You must enter at least one pair.";
}#else (no files)

return $messg;
}#check_pairs_in_list


#----------------------------------------------------------
sub check_files_in_list
{
my $file; #holds one filename from @_
my $found; #flag
my $i; #loop counter
my $messg = "";


#check that there are filenames to check
if (not(scalar @_ == 0))
{
#check that each file (item in @_) is in the reads_files (linear search)
foreach $file (@_)
{
$found = 0; #reset flag
$i = 0; #reset loop counter

#while we haven't found it in and we haven't fallen off the array
while (not($found) && $i < scalar @reads_files)
{
#if it matches the current filename
if (trim($file) eq trim($reads_files[$i]))
{
#toggle flag
$found = 1; 
}#if (match)

$i++; #increment loop counter
}#while (search)

#if not found, save error message and break out of the loop
if (not($found))
{
$messg = "The file \"".trim($file)."\" was not found in list of reads files.";
last; #(so sorry)
}#if (not found)

}#foreach (file)
}#if (files to check)
else
{
$messg = "You must enter at least one file.";
}#else (no files)

return $messg;
}#sub_check_files_in_list


#----------------------------------------------------------
sub get_cd_labels_and_files
{
my $user_resp; #hold user response to prompts (y/n)
my $user_check = 0; #flag for final check
my $label_check = 0; #flag for checking labels entry
my @labels; #holds the labels as an array so we can assign files to them
my $labels_str; #holds current labels in comma separated string
my $files_str = ""; #holds current files in comm/space seperated string
my $curr_label; #foreach loop variable, holds one label


print "\n\n###################################
Cuffdiff

(More info: http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html)\n";

#check is the user just wants the same labels and files as with Cuffnorm: it's probable and it'll save quite a bit of prompting
print "\nDo you want to use the same labels and files as the ones you entered for Cuffnorm? (y/n)\n";
trim($user_resp = <>);

if ($user_resp =~ m/^[Yy]/)
{
$cd_labels = $cn_labels;
$cd_files = $cn_files;
}#if (use same files/labels


#otherwise, prompt away!
else
{
#loop until all information (labesl and files) is correct
do
{

#loop until we get the labels we want
do
{
#Prompt for cn_labels
print "\nEnter the label names to use for Cuffdiff in a comma separated list:
(These will be used to group the treatments/bioreps you want to compare)\n";

#read labels from user as string
trim($labels_str = <>);

#split to array
@labels = split(",",$labels_str);

#check with user that they want to use these labels
print "\n\nThese are the label names you entered for Cuffdiff:
@labels
Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$label_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
}#else (user check)

}while(not($label_check));


#for each label, loop until valid: prompt for files in label (print list of included files without ext.)
foreach $curr_label (@labels)
{

$files_str .= get_label_files($curr_label)." ";
#once label info is good tack on to cn_files string in proper format (comma separated, space at end)

}#foreach (label)


#ask: corrections? (display labels and files: if yes re-do, else return)
print "\n\nThese are the label names and files you entered for Cuffdiff:
(Note: If you selected PE settings you will see a '/' between the files in a pair)
labels: ".trim($labels_str)."
files: ".trim($files_str)."

Is this information correct? (y/n)
(Press enter for yes)\n";

#read response
trim($user_resp = <>);

#(no/n means redo loop, everything else means okay)
if (not($user_resp =~ m/^[Nn]/))
{
#we're good! toggle inner loop flag
$user_check = 1;
}#if (user check)
else
{
printf "\nYou said the information was incorrect so try again.\n";
#reset file_str
$files_str ="";
}#else (user check)


}while(not($user_check));

#save cd_labels and cd_files
$cd_labels = $labels_str;
$cd_files = $files_str;

}#else (dif labels from cuffnorm)

}#get_cd_labels_and_files

 
#----------------------------------------------------------
sub run_trim
{
my $all_trim_cmds; #holds return from build_trim_cmds function

#message: starting trims
print "\n\n\n###################################
Starting trims...\n\n";

#unless we're only building the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#build the trimmomatic output dir otherwise it won't be happy...
mkdir(TRIM_OUT_DIR); #mkdir won't overwrite existing dir and won't crash if it does exist
}#unless (only build)

#get the trim commands built
$all_trim_cmds = build_trim_cmds();

#Write these cmds to the cmds_file
print CMDS_FILE "\nTrimmomatic commands:\n".$all_trim_cmds."\n";

#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#then run the trim commands
system("$all_trim_cmds");
}#if (run_cmds)

#message: trims done
print "\n\nTrims done.\n";
}#run_trim


#----------------------------------------------------------
sub build_trim_cmds
{
my $trim_cmd; #hold the full trim cmd assembled from part 1, 2 and the input output file, for one file
my $io_files; #holds the string of input and output files
my $all_trim_cmds = ""; #hold all the cmds for all files 
my $i; #loop counter

#build the constant part of the trim command (2 parts because the input files are in the middle)
my $trim_cmd_part1 = "java -jar ".TRIMMOMATIC_JAR." $end_type -threads $threads -phred33";
my $trim_cmd_part2 = "ILLUMINACLIP:$adapters_file$adapters_settings HEADCROP:$headcrop LEADING:$leading TRAILING:$trailing SLIDINGWINDOW:$sliding_window MINLEN:$minlen AVGQUAL:$avgqual";


#get the right io files depending on end_type settings

#Single End
if ($end_type =~ m/SE/)
{
#for each input file (for SE only for now)
for ($i = 0; $i < scalar @reads_files_ext; $i++)
{
#get the io files string built
$io_files = build_SE_trim_io_files($i);

#build the full cmd for current file
$trim_cmd = $trim_cmd_part1." ".$io_files." ".$trim_cmd_part2.";";

#and add it to the list of full trim commands
$all_trim_cmds .= $trim_cmd." ";
}#for
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
#for each input file (for SE only for now)
for ($i = 0; $i < scalar @reads_files_pairs; $i++)
{
#get the io files string built
$io_files = build_PE_trim_io_files($i);

#build the full cmd for current file
$trim_cmd = $trim_cmd_part1." ".$io_files." ".$trim_cmd_part2.";";

#and add it to the list of full trim commands
$all_trim_cmds .= $trim_cmd." "." &>> trim_log.txt";
}#for 
}#elsif (PE)


#this should never be triggered but just in case...
#^like it'll never be triggered unless someone changes the setters for end settings in the code
else
{
$all_trim_cmds = "echo '\n\nYou messed with something you shouldn't have didn't you?\n\n'";
}#else

return $all_trim_cmds
}#build_trim_cmds


#----------------------------------------------------------
sub build_SE_trim_io_files #pass it the $i index, returns an io_files string
{
my $io_files; #holds the string of input and output files 
my $input_file; #holds the string of input files
my $output_file; #holds the string of output files


#hold the input filename
$input_file = $reads_files_ext[$_[0]];

#build the output filename (including ./out_dir/filename)
$output_file = TRIM_OUT_DIR."/".$reads_files[$_[0]].TRIM_OUT_EXT; #(.fq is trimmmed extension)

#put the input and output files together
$io_files = $input_file." ".$output_file;

#and return that
return $io_files;
}#build_SE_trim


#----------------------------------------------------------
sub build_PE_trim_io_files #pass it the $i index, returns an io_files string
{
my $io_files; #holds the string of input and output files 
my $input_files; #holds the string of input files
my $output_files; #holds the string of output files
my @pair; #hold the split pair of input filenames without extension

#split the pair
@pair = split(" ",$reads_files_pairs[$_[0]]);

#build the string of input files
$input_files = $pair[0].TRIM_INPUT_EXT." ".$pair[1].TRIM_INPUT_EXT;

#build the string of output files
$output_files = TRIM_OUT_DIR."/".$pair[0].TRIM_PAIRED_EXT.TRIM_OUT_EXT." ".TRIM_OUT_DIR."/".$pair[0].TRIM_UNPAIRED_EXT.TRIM_OUT_EXT." ".TRIM_OUT_DIR."/".$pair[1].TRIM_PAIRED_EXT.TRIM_OUT_EXT." ".TRIM_OUT_DIR."/".$pair[1].TRIM_UNPAIRED_EXT.TRIM_OUT_EXT;

#put them together
$io_files = $input_files." ".$output_files;

#and return that
return $io_files;
}#build_PE_trim


#----------------------------------------------------------
sub run_bowtie
{
my $bowtie_cmd; #hold bowtie2-build cmd to run

#message building index...
print "\n\n\n###################################
Building bowtie index...\n\n";

#build the cmd
$bowtie_cmd = "bowtie2-build ".trim($fasta_file)." ./$genome_index_base";

#Write this cmd to the cmds_file
print CMDS_FILE "\nBowtie2-build command:\n".$bowtie_cmd."\n";


#only check two of the files because yeah... (check the two possible kinds of extensions)
#index not built 
if ((-f $genome_index_base.".1.bt2" && -f $genome_index_base.".rev.2.bt2") || (-f $genome_index_base.".1.bt2l" && -f $genome_index_base.".rev.2.bt2l"))
{
print "\n\nIndex was already built.
(If you want it rebuilt you must delete the existing files)\n";

#If the commands won't be run make a note of it in the commands file
print CMDS_FILE "\nThe Bowtie index is already built so the previous command does not need to be run.\n";
}#if (index already built)


else #index not built 
{
#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#run
system("$bowtie_cmd");
}#unless (run_cmds)

#message: index done
print "\n\nIndex built.\n";

}#else (index not built) 
}#run_bowtie


#----------------------------------------------------------
sub run_tophat
{
my $all_tophat_cmds; #holds return from build_tophat_cmds function

#message: starting tophat
print "\n\n\n###################################
Starting alignments...\n\n";

#get the trim commands built (also runs the cmd to set up the transcriptome index)
$all_tophat_cmds = build_tophat_cmds();

#Write these cmds to the cmds_file
print CMDS_FILE "\nTophat2 commands:\n".$all_tophat_cmds."\n";

#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#then run the trim commands
system("$all_tophat_cmds");
}#unless (run_cmds)

#message: tophat done
print "\n\nAlignments done.\n";
}#run_tophat


#----------------------------------------------------------
sub build_tophat_cmds
{
my $setup_cmd; #hold the setup cmd to build the transcriptome index
my $tophat_cmd_part1; #hold the constant part of the tophat cmd (after setup)
my $tophat_cmd; #hold one complete tophat command
my $input_file; #hold the name of one input file
my $all_tophat_cmds = ""; #hold all of the tophat cmds
my $i; #loop counter
my @pair; #holds a split pair of reads filenames
my $combined_tophat_out_dir; #name of the combined output directory to place the alignment in


#build the constant part of the cmd
$tophat_cmd_part1 = "tophat2 ".TOPHAT_CONSTANT_OPTIONS." -p $threads --transcriptome-index ".TOPHAT_OUT_DIR."/$genome_index_base -o ".TOPHAT_OUT_DIR."/";

#build the setup cmd
$setup_cmd = $tophat_cmd_part1." -G ".trim($genome_file)." $genome_index_base";

print CMDS_FILE "\nTophat2 build transcriptome index command:\n".$setup_cmd."\n";


#get the right io files depending on end_type settings

#Single End
if ($end_type =~ m/SE/)
{
#for each input file (for SE only for now)
for ($i = 0; $i < scalar @reads_files_ext; $i++)
{
#get the input filename (don't forget the dir and the extension)
$input_file = TRIM_OUT_DIR."/".$reads_files[$i].TRIM_OUT_EXT;

#build the cmd
$tophat_cmd = $tophat_cmd_part1.$reads_files[$i].TOPHAT_OUT_DIR_EXT." ".$genome_index_base." ".$input_file.";";

#and tack it onto all_tophat_cmds
$all_tophat_cmds .= $tophat_cmd." ";
}#for
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
#for each input file (for SE only for now)
for ($i = 0; $i < scalar @reads_files_pairs; $i++)
{
#get the io files string built
$input_file = build_PE_tophat_input_files($i);

#split the current pair to make it easier to create the sub output directory name
@pair = split(" ",$reads_files_pairs[$i]);

#build the name of the combined output directory
$combined_tophat_out_dir = $pair[0]."_and_".$pair[1];

#add this directory name to a list for future reference
push (@combined_dirs, $combined_tophat_out_dir);

#build the cmd
$tophat_cmd = $tophat_cmd_part1.$combined_tophat_out_dir.TOPHAT_OUT_DIR_EXT." ".$genome_index_base." ".$input_file.";";

#and add it to the list of full trim commands
$all_tophat_cmds .= $tophat_cmd." ";
}#for 
}#elsif (PE)


#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#run the setup command to build the transcriptome index
system("$setup_cmd");
}#unless (run_cmds)

#return the other ones
return $all_tophat_cmds;
}#build_tophat_cmds


#----------------------------------------------------------
sub build_PE_tophat_input_files 
#pass it the current file pair and it'll return the tophat input files as a formatted string
{
my $input_files; #holds the list of input files
my @pair; #hold the split pair of input filenames without extension

#split the current pair
@pair = split(" ",$reads_files_pairs[$_[0]]);

#build the input file list (r1_paired,r1_unpaired,r2_unpaired r2_paired)
$input_files = TRIM_OUT_DIR."/".$pair[0].TRIM_PAIRED_EXT.TRIM_OUT_EXT.",".TRIM_OUT_DIR."/".$pair[0].TRIM_UNPAIRED_EXT.TRIM_OUT_EXT.",".TRIM_OUT_DIR."/".$pair[1].TRIM_UNPAIRED_EXT.TRIM_OUT_EXT." ".TRIM_OUT_DIR."/".$pair[1].TRIM_PAIRED_EXT.TRIM_OUT_EXT;

#and return it
return $input_files;
}#build_PE_tophat_input_files


#----------------------------------------------------------
sub concatenate_align_summaries
{
my $curr_file_path; #hold current path to the align summary we're getting
my $file; #hold current file we're getting the align summary from
my $pair; #hold a comdbined dir pair
my $line; #hols a line from a summary file

#for f in *.txt; do (echo "$f"; cat $f; echo '') >> ALIGN_SUM_FILENAME; done

#open the output file
open(ALL_SUMS, ">".TOPHAT_OUT_DIR."/".ALL_ALIGN_SUMS_FILENAME) or die "Couldn't open ".ALL_ALIGN_SUMS_FILENAME."\n";

#print a title
print ALL_SUMS "This file contains all the align_summary.txt's concatenated\n\n";


#Single End
if ($end_type =~ m/SE/)
{
foreach $file (@reads_files)
{
#build current path to the align summary file
$curr_file_path = TOPHAT_OUT_DIR."/".$file.TOPHAT_OUT_DIR_EXT."/".ALIGN_SUM_FILENAME;

#print name of the reads file to the output file
print ALL_SUMS $file."\n";

#open the current align_summary file  
open(SINGLE_SUMMARY, "<", $curr_file_path);

#print all of its line to the all_sums file
while ($line = <SINGLE_SUMMARY>)
{
print ALL_SUMS $line;
}#while (more lines in file)

print ALL_SUMS "\n\n";

close(SINGLE_SUMMARY) or die "Couldn't close $curr_file_path\n";
}#foreach (reads_files)
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
foreach $pair (@combined_dirs)
{
#build current path to the align summary file
$curr_file_path = TOPHAT_OUT_DIR."/".$pair.TOPHAT_OUT_DIR_EXT."/".ALIGN_SUM_FILENAME;

#print name of the reads file to the output file
print ALL_SUMS $pair."\n";

#open the current align_summary file  
open(SINGLE_SUMMARY, "<", $curr_file_path);

#print all of its line to the all_sums file
while ($line = <SINGLE_SUMMARY>)
{
print ALL_SUMS $line;
}#while (more lines in file)

print ALL_SUMS "\n\n";

close(SINGLE_SUMMARY) or die "Couldn't close $curr_file_path\n";
}#foreach (pair)
}#elsif (PE)


#close file
close(ALL_SUMS) or die "Couldn't close ".ALL_ALIGN_SUMS_FILENAME;

}#concatenate_align_summaries


#----------------------------------------------------------
sub run_cuffquant
{
my $all_cq_cmds; #holds return from build_cq_cmds function

#message: starting cuffquant
print "\n\n\n###################################
Starting Cuffquant...\n\n";

#get the cuffquant commands built 
$all_cq_cmds = build_cq_cmds();

#Write these cmds to the cmds_file
print CMDS_FILE "\nCuffquant commands:\n".$all_cq_cmds."\n";

#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#then run the cuffquant commands
system("$all_cq_cmds");
}#unless (run_cmds)

#message: cuffquant done
print "\n\nCuffquant done.\n";
}#run_cuffquant


#----------------------------------------------------------
sub build_cq_cmds
{
my $cq_cmd; #hold one cuffquant command
my $all_cq_cmds = ""; #hold all the cuffquant commands
my $cq_cmd_part1; #hold the constant part of the cuffquant cmd
my $input_file; #hold the name of one input file
my $i; #loop counter

#build the constant part
$cq_cmd_part1 = "cuffquant -p $threads -o ".CQ_OUT_DIR."/";


#Single End
if ($end_type =~ m/SE/)
{
#for each file in @reads_files
for ($i = 0; $i < scalar @reads_files; $i++)
{
#build ipnut filename
$input_file = TOPHAT_OUT_DIR."/".$reads_files[$i].TOPHAT_OUT_DIR_EXT."/".CQ_INPUT_FILE;

#build cq_cmd
$cq_cmd = $cq_cmd_part1.$reads_files[$i]." ".trim($genome_file)." ".$input_file.";";

#tack it onto all_cq_cmds
$all_cq_cmds .= $cq_cmd." ";

}#for
}#if (SE)



#Paired End
elsif ($end_type =~ m/PE/)
{
#for each file in combined_dirs
for ($i = 0; $i < scalar @combined_dirs; $i++)
{
#build the input filename
$input_file = TOPHAT_OUT_DIR."/".$combined_dirs[$i].TOPHAT_OUT_DIR_EXT."/".CQ_INPUT_FILE;

#build cq_cmd
$cq_cmd = $cq_cmd_part1.$combined_dirs[$i]." ".trim($genome_file)." ".$input_file.";";

#tack it onto all_cq_cmds
$all_cq_cmds .= $cq_cmd." ";

}#for
}#elsif (PE)


#return the full list of cuffquant cmds
return $all_cq_cmds;
}#build_cq_cmds


#----------------------------------------------------------
sub run_cuffnorm
{
my $cn_cmd; #holds return from build_cn_cmd function

#message: starting cuffnorm
print "\n\n\n###################################
Starting Cuffnorm...\n\n";

#get the cuffnorm commands built 
$cn_cmd = build_cn_cmd();

#Write this cmd to the cmds_file
print CMDS_FILE "\nCuffnorm command:\n".$cn_cmd."\n";

#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#then run the cuffnorm commands
system("$cn_cmd");
}#unless (run_cmds)

#message: cuffnorm done
print "\n\nCuffnorm done.\n";
}#run_cuffnorm


#----------------------------------------------------------
sub build_cn_cmd
{
my $cn_cmd; #hold the cuffnorm command
my $input_file_list; #hold the formatted list of input files


#build the input file list
$input_file_list = build_label_file_list($cn_files);

#build the full cuffnorm cmd
$cn_cmd = "cuffnorm -p $threads --labels ".trim($cn_labels)." -o ".CN_OUT_DIR." ".trim($genome_file)." $input_file_list";

return $cn_cmd;
}#build_cn_cmd


#----------------------------------------------------------
sub build_label_file_list #pass it the unformatted files list (cn or cd)
{
my $input_files = ""; #hold the formatted list of all input files


#Single End
if ($end_type =~ m/SE/)
{
$input_files = build_SE_label_file_list($_[0]); #first and only arg should be cn or cd files
}#if (SE)


#Paired End
elsif ($end_type =~ m/PE/)
{
$input_files = build_PE_label_file_list($_[0]); #first and only arg should be cn or cd files
}#elsif (PE, pairs)


return $input_files;
}#build_label_file_list


#----------------------------------------------------------
sub build_PE_label_file_list
{
my @label_groups; #holds an array of the groups of files belonging to each label
my @label_pairs; #holds an array of filenames (no ext) belonging to one label
my $label_pair; #holds one pair of files from one label group
my $input_file; #hold one input filename
my $input_files = ""; #hold the formatted list of all input files
my $label_group; #hold one group of labels
my $comma_flag; #so we don't add a comma after the last file in a label group
my $i; #loop iterator
my @pair_files; #holds the split of the files within a pair


#split cn_files by space into label groups
@label_groups = split(" ",$_[0]); #this paramter should be $cn_files or $cd_files

#foreach label group 
foreach $label_group (@label_groups)
{
#split by "," into pairs list
@label_pairs = split(",",$label_group);

#reset comma flag
$comma_flag = 0;

foreach $label_pair (@label_pairs)
{
#split by "/" into files of a pair
@pair_files = split("/",$label_pair);

#build the path to the right input file
$input_file = CQ_OUT_DIR."/".$pair_files[0]."_and_".$pair_files[1]."/".CND_INPUT_FILE;

#add a comma to the list of files for this label if it's not the first file in the label
if ($comma_flag)
{
$input_files .= ",";
}#if (add comma)

#otherwise toggle the flag because the next file won't be the first
else
{
$comma_flag = 1;
}#else

#tack this input file onto the input file list
$input_files .= $input_file;

}#foreach (label_pair)

#add a space between label groups
$input_files .= " ";

}#foreach (label_group)

return $input_files;

}#build_PE_label_file_list


#----------------------------------------------------------
sub build_SE_label_file_list
{
my @label_groups; #holds an array of the groups of files belonging to each label
my @label_files; #holds an array of filenames (no ext) belonging to one label
my $input_file; #hold one input filename
my $input_files = ""; #hold the formatted list of all input files
my $label_group; #hold one group of labels
my $label_file; #hold one label file
my $comma_flag; #so we don't add a comma after the last file in a label group


#split cn_files by space into label groups
@label_groups = split(" ",$_[0]); #this paramter should be $cn_files or $cd_files

#foreach label group 
foreach $label_group (@label_groups)
{
#split by "," into files list
@label_files = split(",",$label_group);

#reset comma flag
$comma_flag = 0;

#foreach file in a label
foreach $label_file (@label_files)
{
#build the path to the right input file
$input_file = CQ_OUT_DIR."/".$label_file."/".CND_INPUT_FILE;

#add a comma to the list of files for this label if it's not the first file in the label
if ($comma_flag)
{
$input_files .= ",";
}#if (add comma)

#otherwise toggle the flag because the next file won't be the first
else
{
$comma_flag = 1;
}#else

#tack this input file onto the input file list
$input_files .= $input_file;

}#foreach (file)

#add a space between label groups
$input_files .= " ";

}#foreach (label_group)

return $input_files;
}#buil_SE_label_file_list


#----------------------------------------------------------
sub run_cuffdiff
{
my $cd_cmd; #holds return from build_cd_cmd function

#message: starting cuffdiff
print "\n\n\n###################################
Starting Cuffdiff...\n\n";

#get the cuffdiff commands built 
$cd_cmd = build_cd_cmd();

#Write this cmd to the cmds_file
print CMDS_FILE "\nCuffdiff command:\n".$cd_cmd."\n";

#run the cmds unless the user chose to only build the cmds
unless ($run_cmds =~ m/^[Bb]/)
{
#then run the cuffdiff commands
system("$cd_cmd");
}#unless (run_cmds)

#message: cuffdiff done
print "\n\nCuffdiff done.\n";
}#run_cuffdiff


#----------------------------------------------------------
sub build_cd_cmd
{
my $cd_cmd; #hold the cuffnorm command
my $input_file_list; #hold the formatted list of input files

#build the input file list
$input_file_list = build_label_file_list($cd_files);

#build the full cuffnorm cmd
$cd_cmd = "cuffdiff -p $threads -labels ".trim($cd_labels)." -o ".CD_OUT_DIR." ".trim($genome_file)." $input_file_list";

return $cd_cmd;
}#build_cd_cmd


#----------------------------------------------------------
sub print_end_message
{
print "\n
********************
* End of pipeline. *
********************\n\n";
}#print_end_message


#----------------------------------------------------------
sub print_cmds_files_end
{
#get the end time
my $date_string = localtime(); #hold the current date and time

#print it to the file
print CMDS_FILE "\n\nEnd date and time: $date_string\n\n";
}








