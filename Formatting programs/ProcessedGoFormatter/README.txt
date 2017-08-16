ProcessedGoFormatter usage

To run: 

$GOFormat inputFileName outputFilename (once you've added go format to a bin somewhere, or put it on your path)

-OR-

$java ProcessedGoFormatter inputFilename outputFilename

inputFilename => name/path of input file
outputFilename => name/path of output file


Notes:
-These file names are expected and necessary. The program will tell you it's missing file names.
-If one of your files cannot be opened (incorrect name/path or does not exist) The program will tell you it cannot find the file.
-GO terms are filtered by p-value: anything greater than 0.001 is excluded from the output file.
-The program will print "End of processing." to the terminal when it's done.

Errors:
-If the program is not compiled (no .class file present, only .java),
before running it compile it with:

$javac ProcessedGoFormatter.java


-If the program prints several lines (stack trace), the exception (error) type will be at the top and will (hopefully) give you an idea of what went wrong. These kinds of error will most likely arise if the format of the input file is wrong, since the program expects a very specific format.
