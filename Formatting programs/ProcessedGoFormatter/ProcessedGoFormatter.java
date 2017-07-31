/********************************************************************
 * Class: ProcessedGoFormatter
 * 
 * Author: Nadège Pulgar-Vidal
 * Date: June 3, 2016
 * Purpose: Format go.processed.txt files outputed by SeqEnrich:
 *          Filter by p-value and print groups vertically rather 
 *          than on the same line.
 * 
 ********************************************************************/
import java.io.*;
import java.util.Scanner;

public class ProcessedGoFormatter 
{
  public static final double P_VALUE_THRESHOLD = 0.001; //p-val threshold for filtering 
  
  //position of elements in a tokenized line of input
  public static final int GO_ID = 0;
  public static final int GO_DESCRIPTION = 1;
  public static final int P_VALUE = 6;
  public static final int GENES_START = 7;
  
  public static void main(String[] args)
  {
    Scanner inFile; //GO term file to read from
    PrintWriter outFile; //Output file
    String inputFilename; //To be obtained as a cmd line arg
    String outputFilename; //To be obtained as a cmd line arg
    String line; //Holds a line from the input file
    String[] lineTokens; //Holds tokenized line
    String goID; //Holds the GO term ID
    String goDescription; //Hold the Go term description
    double pValue; //Hold the GO term's p-value
    String geneAccession; //Holds a gene accession
    
    
    if (args.length < 2)
    {
      System.out.println("Missing file names.");
    }//if
    
    else
    {
      inputFilename = args[0];
      outputFilename = args[1];
      
      try
      {
        inFile = new Scanner(new File(inputFilename));
        outFile = new PrintWriter(new File(outputFilename));
        
        //read every line from the input file
        while (inFile.hasNext())
        {
          line = inFile.nextLine();
          
          //tokenize line
          line = line.trim();;
          lineTokens = line.split("\t");
          
          //System.out.println(line);
          
          //extract relevant info
          goID = lineTokens[GO_ID];
          goDescription = lineTokens[GO_DESCRIPTION];
          pValue = Double.parseDouble(lineTokens[P_VALUE]);
          //print this first then extract gene accessions since there's a variable amount
          //System.out.println(pValue);
          
          //filter out less significant p-values
          if (pValue <= P_VALUE_THRESHOLD)
          {
            outFile.println(goID + " : " + goDescription + "\n" + pValue);
            //System.out.println(goID + " : " + goDescription + "\n" + pValue);
            
            for (int i = GENES_START; i < lineTokens.length; i++)
            {
              geneAccession = lineTokens[i];
              outFile.println(geneAccession);
              //System.out.println(geneAccession);
            }//for
            
            outFile.println();
          }//if
                 
        }//while
        
        inFile.close();
        outFile.close();
        
        System.out.println("End of processing.");
      }//try
      
      catch (FileNotFoundException e)
      {
        System.out.println("Error: File not found.");
      }//catch
      
      catch (Exception e)
      {
        e.printStackTrace();
      }//catch
      
    }//else
    
  }//main
}//ProcessedGoFormatter






















