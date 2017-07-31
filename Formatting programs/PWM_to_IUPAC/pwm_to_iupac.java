/******************************************************************************
 * PROGRAM: pwm_to_iupac
 * AUTHOR: Nadege Pulgar-Vidal (pulgarvn@myumanitoba.ca)
 * DATE: May 10th, 2016
 * PURPOSE: Take PWM text files as a command line argument and output
 *          the corresponding IUPAC sequences to stdout.
 *          There are some format requirements for two basic text layouts.
 *          Please refer to the usage file for more information
 ******************************************************************************/
import java.util.*; 
import java.io.*;

public class pwm_to_iupac
{
  //---------------------------------CONSTANTS---------------------------------
  
  public static final String A = "A";
  public static final String C = "C";
  public static final String G = "G";
  public static final String T = "T";
  
  //Prime number values for each base for constant time code lookup
  //If you change these please make sure that there are no collisions 
  // among sum combinations (plz don't change it)
  public static final int A_VALUE = 3;
  public static final int C_VALUE = 5;
  public static final int G_VALUE = 7;
  public static final int T_VALUE = 11;
  
  public static final int[] BASE_VALUES = {A_VALUE, C_VALUE, G_VALUE, T_VALUE};
  
  public static final int NUM_BASES = 4;
  
  public static final int CODES_ARR_LEN = A_VALUE + C_VALUE + G_VALUE + T_VALUE;
  
  public static final char[] IUPAC_CODES = {'-','-','-','A','-', //0-4
                                            'C','-','G','M','-', //5-9
                                            'R','T','S','-','W', //10-14
                                            'V','Y','-','K','H', //15-19
                                            '-','D','-','B','-', //20-24
                                            '-','N'};             //25-26
  
  public static final char HORIZONTAL = 'H';
  public static final char VERTICAL = 'V';
  
  private static double sigThreshold = 0.1; //significance threshold

  //------------------------------------MAIN-----------------------------------
  public static void main(String[] args) throws Exception
  {
    //Command line arguments or text-based interaction?
    //Command line arguments for now, becasue error checking 
    
    //Get file name and send it to read method
    String inputFilename = args[0];
    
    //Is the file horizontal or vertical?
    char type = args[1].charAt(0);
    
    //If the user specified a threshold value
    if (args.length > 2)
    {
      sigThreshold = Double.parseDouble(args[2]);
      
      //Check the validity of the threshold value to notify the user of errors
      if (sigThreshold > 1 || sigThreshold < 0)
      {
        throw new Exception("Invalid threshold value.");
      }//if
    }//if
    
    //Do the read an output the sequences
    if (Character.toUpperCase(type) == HORIZONTAL)
    {
      convertHorizontal(inputFilename);
    }//if
    
    else if (Character.toUpperCase(type) == VERTICAL)
    {
      convertVertical(inputFilename);
    }//else if
    
    else
    {
      throw new Exception("Invalid conversion type.");
    }//else
  }//main
  
  
  //----------------------------------METHODS----------------------------------
  
  //-------------------------------------------------------------------
  // convertHorizontal (the probs extend horizontally)
  // PURPOSE: Reads the probabilities from a file in the expected format
  //          and outputs the corresponding iupac sequence.
  //          Expects the probabilities of each base on a seperate line 
  //          ex: A   P(A1) P(A2) P(A3)
  //              C   P(C1) P(C2) P(C3)
  //              G   P(G1) P(G2) P(G3)
  //              T   P(T1) P(T2) P(T3)
  // PARAMETERS: The text file's name
  //-------------------------------------------------------------------
  private static void convertHorizontal(String inputFilename)
  {
    String motifName; //The name of the current motif
    String iupacSequence; //The IUPAC sequenc of the current motif
    double[] probA, probC, probG, probT; //The probabilities of the dif bases
    Scanner inputFile; //The input file
    Scanner lineScanner; //Parse one line of input
    String line; //The current line from the file
    String[] tokens; //Holds a tokenized line
    
    //So the compiler doesn't complain and so that we don't build a sequence 
    //if something is wrong 
    probA = new double[0];
    probC = new double[0];
    probG = new double[0];
    probT = new double[0];
    
    try
    {
      //Open the file to read from
      inputFile = new Scanner(new File(inputFilename));
      
      //Until we reach the end of the file
      while (inputFile.hasNext())
      {
        //The first line should be the motif name
        motifName = inputFile.nextLine().trim();
        
        //Grab the next line
        line = inputFile.nextLine();
        line = line.trim();
        
        //There might be a blank line between the motif and the probabilities  
        //If that is the case skip to the next line
        if (line.length() == 0)
        {
          line = inputFile.nextLine();
          line = line.trim();
        }//if
        
        //Then reset the current sequence string
        iupacSequence = "";
        
        //Keep reading lines until
        while (line.length() > 0)
        {
          //Tokenize line by whitespace
          tokens = line.split("\\s+");
          
          //Assuming the base is the first token
          //Create the probability arrays
          if (tokens[0].equalsIgnoreCase(A))
          {
            probA = makeProbArray(tokens);
          }//if
          
          else if (tokens[0].equalsIgnoreCase(C))
          {
            probC = makeProbArray(tokens);
          }//else if
          
          else if (tokens[0].equalsIgnoreCase(G))
          {
            probG = makeProbArray(tokens);
          }//else if
          
          else if (tokens[0].equalsIgnoreCase(T))
          {
            probT = makeProbArray(tokens);
          }//else if         
          
          //Move on to the next line if there is one
          if (inputFile.hasNextLine())
          {
            line = inputFile.nextLine();
            line = line.trim();
          }//if
          
        }//while
        
        //If we weren't able to get the probability arrays, the file 
        //did not have the right format
        if (probA.length == 0)
          throw new InputMismatchException();
        
        //Once we have the probability arrays build the IUPAC sequence
        for (int i = 0; i < probA.length; i++)
        {
            iupacSequence += computeCode(probA[i], probC[i], probG[i], probT[i]);
        }//for
        
        //Once we have a full sequence print the motif_id, the motif length, and the sequence
        System.out.printf("%s   Length: %d\n%s\n\n", motifName, iupacSequence.length(), iupacSequence);
        
      }//while
      
      //Close the file once done
      inputFile.close();
      
    }//try
    
    //Possible exceptions
    catch (FileNotFoundException e)
    {
      System.out.println("File not found.");
    }//catch
    catch (InputMismatchException e)
    {
      //A tokens mismatch means the format was different...
      System.out.println("Unexpected file format.");
    }//catch
    catch (NoSuchElementException e)
    {
      System.out.println("File ended unexpectedly.");
    }//catch
    catch (Exception e)
    {
      System.out.println(e);
    }//catch
    
  }//convertHorizontal
  
  //-------------------------------------------------------------------
  // makeProbArray
  // PURPOSE: Given an array of tokens and a double array variable 
  //          create an array of probabilities in doubles from the tokens.
  // PARAMETERS: An array of tokens 
  // RETURNS: An array of doubles
  //-------------------------------------------------------------------
  public static double[] makeProbArray(String[] tokens)
  {
    double[] probs = new double[tokens.length-1];
    
    for (int i = 0; i < probs.length; i++)
    {
      //+1 to skip over the base
      probs[i] = Double.parseDouble(tokens[i+1]); //the tokens should all be valid doubles
    }//for
    
    return probs;
  }//makeProbArray
  
  
  //-------------------------------------------------------------------
  // convertVertical (the probs extend vertically)
  // PURPOSE: Reads the probabilities from a file in the expected format
  //          and outputs the corresponding iupac sequence.
  //          Expects the probabilities to be ordered in A C G T order each line
  //          ex:  1   P(A1)  P(C1)  P(G1)  P(T1)
  //               2   P(A2)  P(C2)  P(G2)  P(T2)  
  //               3   P(A3)  P(C3)  P(G3)  P(T3)  etc.
  //          (line number no longer essential)
  // PARAMETERS: The text file's name
  //-------------------------------------------------------------------
  private static void convertVertical(String inputFilename)
  {
    String motifName; //The name of the current motif
    String iupacSequence; //The IUPAC sequenc of the current motif
    double probA, probC, probG, probT; //The probabilities of the dif bases
    Scanner inputFile; //The input file
    Scanner lineScanner; //Parse one line of input
    String line; //The current line from the file
    
    try
    {
      //Open the file to read from
      inputFile = new Scanner(new File(inputFilename));
      
      //Until we reach the end of the file
      while (inputFile.hasNext())
      {
        //The first line should be the motif name
        motifName = inputFile.nextLine();
        
        //skip blank lines and take the first non-blank line as motif name
        while (motifName.trim().length() == 0)
        {
          motifName = inputFile.nextLine();
        }//while
        
        //If the next line doesn't have a line number it isn't the
        //start of the probabilities yet so skip it
        while (!inputFile.hasNextDouble())
        {
          inputFile.nextLine();
        }//while
        
        //Then reset the current sequence string
        iupacSequence = "";
        
        //Then for each next line until we reach a blank line
        line = inputFile.nextLine();
        line = line.trim();
        
        while (line.length() > 0)
        {
          //Set the line Scanner to this line
          lineScanner = new Scanner(line);
          
          if (lineScanner.hasNext())
          { 
            //Then we skip the line number and then get the probabilities
            // line #    P(A)    P(C)    P(G)    P(T)
            //skip the line number if there is one
            if (lineScanner.hasNextInt())
              lineScanner.nextInt();
            probA = lineScanner.nextDouble();
            probC = lineScanner.nextDouble();
            probG = lineScanner.nextDouble();
            probT = lineScanner.nextDouble();
            
            //Get the IUPAC code and add it to the current sequence
            iupacSequence += computeCode(probA, probC, probG, probT);
          }//if
          
          //Move on to the next line if there is one
          if (inputFile.hasNextLine())
          {
            line = inputFile.nextLine();
            line = line.trim();
          }//if
        } //while
        
        //Once we have a full sequence print the motif_id, the motif length, and the sequence
        System.out.printf("%s   Length: %d\n%s\n\n", motifName, iupacSequence.length(), iupacSequence);
        
      }//while
      
      //Close the file once done
      inputFile.close();
    }//try
    
    //Possible exceptions
    catch (FileNotFoundException e)
    {
      System.out.println("File not found.");
    }//catch
    catch (InputMismatchException e)
    {
      //A tokens mismatch means the format was different...
      System.out.println("Unexpected file format.");
    }//catch
    catch (NoSuchElementException e)
    {
      System.out.println("File ended unexpectedly.");
    }//catch
    catch (Exception e)
    {
      System.out.println(e);
    }//catch

  }//convertVertical
  
  
  //-------------------------------------------------------------------
  // computeCode
  // PURPOSE: Given probabilities (doubles) for A,C,G,T, in that order
  //          calculate the corresponding iupac code and return that char.
  //          This method does not check the validity of the probabilities.
  // PARAMATERS: Four doubles representing the probabilitis of A,C,G, and T
  // RETURNS: The iupac code character or - if it was zero or invalid somehow
  //-------------------------------------------------------------------
  private static char computeCode(double probA, double probC, double probG, double probT)
  {
    final double[] probabilities= {probA, probC, probG, probT}; //Grouped for easier setup
    int significantSum = 0; //Sum of values of significant bases
    
    //Compute sum of the values of the significant bases 
    for (int i = 0; i < NUM_BASES; i++)
    {
      //If base is significant
      if (probabilities[i] >= sigThreshold)
      {
        //Add its value to our sum
        significantSum += BASE_VALUES[i];
      }//if
    }//for
    
    //Once we have the sum, look up and return the corresponding code
    return IUPAC_CODES[significantSum];
  }//computeCode
  
  
}//class








