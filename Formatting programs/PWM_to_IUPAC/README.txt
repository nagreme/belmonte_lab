pwm_to_iupac

@author: Nadège Pulgar-Vidal (pulgarvn@myumanitoba.ca)
@date: May 11, 2016

---------PROGRAM INFORMATION

This program converts text files containing PWMs into their corresponding IUPAC sequence

While it should work with most PWM file formats (provided the correct orientation is selected) I cannot guarantee that it will work with all formats.

A nucleotide’s probability is considered significant if it is greater or equal to 0.1 by default.
This can be changed by specifying a desired threshold value between 0 and 1 when running from terminal.

The general format expected is the motif id or name on the first line then the base probabilities
on following lines, either extending vertically or horizontally. There can be other text lines 
between the motif name and the probabilities but only the first line will be taken as the motif name.
The other lines in between will be skipped.

Additionally:
-Probabilities must be expressed in decimal format and separated by spaces or other whitespace (but not new line characters).
-Individual PWMs must be separated by a blank line. (Extra blank lines between motifs should not be a problem)
-There MUST be a blank line at the end of the file (otherwise it will fall into an infinite loop).

>Vertical format example:   (CONVERSION_TYPE = V or v)

Motif_ID
Pos  A     C      G      T
1   P(A1)  P(C1)  P(G1)  P(T1)
2   P(A2)  P(C2)  P(G2)  P(T2)  
3   P(A3)  P(C3)  P(G3)  P(T3) 
.   .      .      .      .
.   .      .      .      .
.   .      .      .      .
n   P(An)  P(Cn)  P(Gn)  P(Tn)

The prefixed line number is no longer necessary: it's presence is inconsequential. (May 25th)
Exception: If there is a zero it must be expressed with a decimal place (ie 0.000 instead of just 0)
otherwise it will be interpreted as a line number (integer) and will lead to an error.

Vertical example:
M0086_1.02.txt
Pos	A	C	G	T
1	0.287385908013657	0.18498146880198	0.164060209169499	0.363572242713558
2	0.268822931218917	0.245690086061912	0.172342428210344	0.313144588564022
3	0.339996933227662	0.263259777689586	0.151657718112127	0.245085641123334
4	0.153137435214402	0.201021880848225	0.101790628046118	0.5440498088481
5	0.240408827004599	0.109918599076006	0.541436105125234	0.108236708258281
6	0.108236708258281	0.541436105125234	0.109918599076006	0.240408827004599
7	0.5440498088481	0.101790628046118	0.201021880848225	0.153137435214402
8	0.245085641123334	0.151657718112127	0.263259777689586	0.339996933227662
9	0.313144588564022	0.172342428210344	0.245690086061912	0.268822931218917
10	0.363572242713558	0.164060209169499	0.18498146880198	0.287385908013657

M0087_1.02.txt
Pos	A	C	G	T
1	0.210612235083376	0.249639703984584	0.212728923235575	0.327019137696465
2	0.368807863977612	0.131666974575549	0.206494920637766	0.293030240809073
3	0.149879907939015	0.326958862281317	0.179438828162233	0.343722401617435
4	0.59709814845862	0.299490142185392	0.0681752436778613	0.0352364656781267
5	0.0713674230204506	0.0871683673170179	0.0308100675523138	0.810654142110218
6	0.10919780544118	0.00375527292996241	0.882854992808631	0.004191928820227
7	0.0879971867821681	0.740484606287012	0.0159156809400298	0.15560252599079
8	0.743094476064179	0.000719515789834785	0.194576069947403	0.0616099381985832
9	0.192295105365977	0.194227703304044	0.35390594971111	0.259571241618869
10	0.226084737261053	0.210800011799758	0.403796029314724	0.159319221624465

Output:
M0086_1.02.txt   Length: 10
NNNNNNNNNN

M0087_1.02.txt   Length: 10
NNNMTRYRNN

M0464_1.02.txt   Length: 10
NHMAGGGGNN

M0470_1.02.txt   Length: 10
NHRGGGWHNN


>Horizontal format example:   (CONVERSION_TYPE = H or h)

Motif_ID
A   P(A1) P(A2) P(A3) ... P(An)
C   P(C1) P(C2) P(C3) ... P(Cn)
G   P(G1) P(G2) P(G3) ... P(Gn)
T   P(T1) P(T2) P(T3) ... P(Tn)

In the horizontal format the PWM rows must be prefixed with the nucleotide’s letter

Horizontal example:
>MA0271.1	ARG80																					
A		0.89	0	0.161	0	0.44	0															
C		0	0	0	0.301	0	0.256															
G		0	0.301	0	0	0.79	0															
T		0.69	0	0	0	0.69	0.2															
																						
>MA0272.1	ARG81																					
A		0.59	0.6	0.2	0.169	0	0.1	0	0.59													
C		0.12	0	0	0	0.271	0	0.279	0.59													
G		0.83	0.4	0.239	0	0	0	0	0.4													
T		0.44	0.154	0	0	0.1	0.164	0	0.64													
																						
>MA0273.1	ARO80																					
A		0.102	0.33	0.126	0.407	0.39	0.31	0.83	0.4	0.5	0.5	0.1	0.1	0.15	0.67	0.69	0.105	0.269	0.353	0.257	0.236	0.278
C		0.155	0.193	0.232	0.253	0.193	0.198	0.154	0.744	0.231	0.987	0.22	0.5	0.359	0.4	0.68	0.222	0.119	0.21	0.217	0.103	0.252
G		0.531	0.301	0.317	0.203	0.229	0.135	0.336	0.3	0.8	0.1	0.975	0.99	0.295	0.6	0.6	0.372	0.193	0.205	0.178	0.349	0.17
T		0.21	0.174	0.323	0.136	0.186	0.355	0.426	0.185	0.753	0.5	0.1	0.2	0.329	0.318	0.181	0.299	0.418	0.23	0.347	0.31	0.298
		
Output:
>MA0271.1	ARG80   Length: 6
WGACDY

>MA0272.1	ARG81   Length: 8
NDRAYWCN

>MA0273.1	ARO80   Length: 21
NNNNNNNNNNNNNNNNNNNNN


Note: the examples use the default threshold value for significance of 0.1


---------RUNNING FROM TERMINAL

To run the program on this type of file from terminal

$java pwm_to_iupac pwm_filename.txt CONVERSION_TYPE THRESHOLD_VALUE >output_filename.txt

The threshold value is optional. If nothing is specified, default is 0.1
Add the path names if the file is not in your current directory
Sending the output to a file is not necessary but I recommend it
(The default is to print to the terminal window)


---------TROUBLESHOOTING ERRORS

-If something was wrong (invalid file name, unexpected format, etc.) the output file
will contain one line of text with a brief description of the error


-If you forgot to give an input file name and/or a conversion type code as a command line argument you will get an error:
$Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: 0
$	at pwm_to_iupac.main(pwm_to_iupac.java:55)


-If you forgot to give a conversion type code you will get an error:
$Exception in thread "main" java.lang.Exception: Invalid conversion type.
$	at pwm_to_iupac.main(pwm_to_iupac.java:58)


-If you gave an invalid threshold value (not between 0 and 1 inclusive) you will get an error:
$Exception in thread "main" java.lang.Exception: Invalid threshold value.
$	at pwm_to_iupac.main(pwm_to_iupac.java:68)


-If the class file is not present (program is not compiled) you will get an error:
$Error: Could not find or load main class pwm_to_iupac

To fix this, run
$javac pwm_to_iupac.java

Then run the program normally
