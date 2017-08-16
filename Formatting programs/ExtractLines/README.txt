README for ExtractLines

Usage:
$ExtractLines INPUT_FILENAME

Run from command line providing an input_file name to parse

Extract lines containing target

(change these in the script)
The default target is "yes" (TARGET)
The default delimiter is a tab (DELIM)
The default output filename is "extracted_line.txt" (OUTPUT_FILENAME)
By default, the program assumes your input's first line is a header to be included in the output (HEADER)

All of these can be changed in the source code (ExtractLines.pl) 
(More details are included there)

bash shortcut is ExtractLines 

*** Important Note ***
This is an older script: using grep (terminal) and piping the output to a file is easier and more efficient.
