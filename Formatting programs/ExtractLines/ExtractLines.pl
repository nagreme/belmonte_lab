#!/usr/bin/perl

#-----------------------------Description------------------------------
# AUTHOR: Nadège Pulgar-Vidal
# DATE: July 20th 2016
#
# PROGRAM: ExtractLines 
#
# DESCRIPTION: Open file provided by user then write the lines that 
#	       contains target to the output file.
#
# Note: User should change the delimiter and target in the code.
# 	These are called DELIM and TARGET respectively
#----------------------------------------------------------------------

#-------------------------------Includes-------------------------------

use strict;
use warnings; 

use Readonly;
use String::Util qw(trim);
use Scalar::Util qw(looks_like_number);


#------------------------------Constants------------------------------

use constant OUTPUT_FILENAME => "extracted_lines.txt";
use constant DELIM => "\t"; #the delimiter for the input file 
use constant HEADER => 1; #1 is true, 0 is false
#^whether there is a header line in the file that should be conserved regardless of whether or not it contains the target

#\t means tab, \s means space. Please put the delimiter within double quotes.
#If the delimiter is a double quote use this \" inside the double quotes
#Do not add spaces within the double quotes

use constant TARGET => "yes"; #extract lines containing this target

#Change this^ to change what the program looks for in each line
#Keep it in double quotes


#------------------------------Variables------------------------------

my $input_file = $ARGV[0]; #the name of the ipnut file specififed by the user
my $input_line; #holds one line from the input file
my @line_tokens = ""; #holds the line plsit by delimiter
my $token; #hold one toekn from the line
my $num_lines_extracted = 0; #line counter (just to let the user known how many lines were extracted)


#---------------------------------Main--------------------------------

#Open the input and output files
open_files();

print "\nRunning...\n\n";

#If there it a header line at the beginning of the file print it to output regardless
if (HEADER)
{
#read the line
$input_line = <INFILE>;

#print it out
print OUTFILE trim($input_line)."\n";
}#if (header)

#while it has input, read a line 
while ($input_line = <INFILE>)
{
#split by delimiter
@line_tokens = split(DELIM, trim($input_line));

#foreach token in the line
foreach $token (@line_tokens)
{
#check if last token equals the target (case sensitive)
if (trim($token) eq TARGET)
{
#then print it to the output file
print OUTFILE trim($input_line)."\n";

#one more line extracted
$num_lines_extracted++;
}#if (yes)

}#foreach (token in line)

}#while (there's input)

#because we're good programmers
close_files();

print "$num_lines_extracted line(s) extracted\n";

print "\nEnd of processing.\n\n";


#-------------------------Function definitions------------------------

#----------------------------------------------------------
sub open_files
{
#Try opening the provided input file
open(INFILE,"<".$input_file) or die "Couldn't open $input_file: $!\n";

#Open the output file
open(OUTFILE,">".OUTPUT_FILENAME) or die "Couldn't open ".OUTPUT_FILENAME.": $!\n";
}#open_file


#----------------------------------------------------------
sub close_files
{
#Try closing the provided input file
close(INFILE) or die "Couldn't close $input_file: $!\n";

#Close the output file
close(OUTFILE) or die "Couldn't close ".OUTPUT_FILENAME.": $!\n";
}#close_file













