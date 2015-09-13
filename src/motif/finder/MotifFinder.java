//Program: Hafenrichter.java
//Course: COSC 420
//Description: An implementation of the burte force median string algorithm and greedy motif search algorithm
//Author: Brandon Hafenrichter
//Revised: September 15, 2015
//Language: Java
//IDE: Netbeans 8.0.2

//***********************************************************************************************************
//***********************************************************************************************************

//Class Main.java
//Description: Contains the essential methods and calls required to run the two algorithms

//***********************************************************************************************************
//***********************************************************************************************************

package motif.finder;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;


//Method: Main
//Description: Prompts the user for the information required to run the algorithm and runs the algorithm
//Returns: None
public class MotifFinder {

    public static String[] sequence;
    public static ArrayList<CandidateMotif> motifs;

    public static void main(String[] args) {
        getUserAction();
    }

    //***********************************************************************************************************
    //Method: getUserAction
    //Description: The process of getting all of the user input required in order to run each of the 
    //algorithms. Also allows user to exit
    //Parameters: None
    //Calls: retrieveSequence(), searchSequence()
    //Globals: None
    //***********************************************************************************************************
    
    public static void getUserAction() {
        System.out.println("Motif Finder: Brandon Hafenrichter");
        System.out.println("");
        System.out.println("1. Retrieve Sequence from a file");
        System.out.println("2. Search for motif");
        System.out.println("3. Exit the program");
        KeyboardInputClass input = new KeyboardInputClass();
        String decision = input.getKeyboardInput("Type DNA Sequence:");

        switch (decision) {
            case "1": {
                try {
                    retrieveSequence();
                    break;
                } catch (IOException ex) {
                    Logger.getLogger(MotifFinder.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            break;
            case "2":
                searchSequence(sequence);
                break;
            case "3":
                System.exit(-1);
                break;
            default:
                System.out.println("Invalid Input Please try again.");
                getUserAction();
        }
    }
    
    //***********************************************************************************************************
    //Method: searchSequence
    //Description: This is the decider for which algorithm is used to search the input given by the user
    //Parameters: String[] sequence -> The DNA sequence supplied by the user via text file
    //Returns: None
    //Globals: None
    //***********************************************************************************************************
    
    private static void searchSequence(String[] sequence) {
        KeyboardInputClass input = new KeyboardInputClass();
        String length = input.getKeyboardInput("Specify the length of the motif to be found:");
        System.out.println("");
        String algorithm = input.getKeyboardInput("1. Brute Force Median String Search \n2. Greedy Motif Search");

        switch (algorithm) {
            case "1":
                bruteForce(sequence, length, true);
                break;
            case "2":
                greedyMotif(sequence, length);
                break;
        }
    }

    //***********************************************************************************************************
    //Method: retrieveSequence
    //Description: Converts the textfile supplied by the user into a String[] that can be used by the algorithms
    //Parameters: None
    //Returns: None
    //Globals: String[] sequence -> The DNA sequence that will be checked by both algorithms
    //***********************************************************************************************************
    
    private static void retrieveSequence() throws IOException {
        //test
        File file = new File("Haf.txt");
        // creates the file
        file.createNewFile();

        // creates a FileWriter Object
        FileWriter writer = new FileWriter(file);
        // Writes the content to the file
        writer.write(
                "atcacttggaacatcgtctagc\n"
                + "cgttggaaatcgctcgtctagc\n"
                + "atcgatcgtctagcccgtggaa\n"
                + "atcgatcgtccctcggaatagc\n"
                + "atcgccttagaaatcgtctagc\n"
                + "tcgccttgcaaatcgtctagca\n"
                + "atcgatcccttggtagtctagc\n"
                + "tcgatccttggagcgtctagca\n"
                + "ctcgatcgtctagtcttggaag\n"
                + "atctttggaacgatcgtctagc");
        writer.flush();
        writer.close();

        TextFileClass textFile = new TextFileClass();
        textFile.getFileName("Specify the text file to be read:");
        textFile.getFileContents();
        if (textFile.fileName.length() > 0) {
            sequence = textFile.text;
            getUserAction();
        } else {
            System.out.println("Invalid text file. Please try again.");
            retrieveSequence();
        }
    }

    //***********************************************************************************************************
    //Method: bruteForce
    //Description: The first method that searches the DNA sequence via Median String Search
    //Parameters:   String[] sequence -> DNA sequences broken up into strings via array
    //              int length -> The length of the l-mer specified by the user
    //              boolean isPrintResults -> Determines whether we need to print our results, instead of just
    //              returning them for later use (Greedy Algorithm)
    //Returns: bestMotif -> returns the best Motif for the sequence provided
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    
    private static CandidateMotif bruteForce(String[] sequence, String length, boolean isPrintResults) {
        if(isPrintResults){
            System.out.println("Working...(This could take a while)");
        }
        
        if(Integer.parseInt(length) > sequence[0].length()){
            System.out.println("Motif length cannot be larger than the number of elements in the DNA strand.  Please try again.");
            return new CandidateMotif("",0);
        }
        
        //initializes the candidate motifs with scores of 0
        long startTime = new Date().getTime();
        getConsensusMotifs(sequence, Integer.parseInt(length));

        //keeps track of how many operations are done
        int operations = 0;

        //iterate through the sequences and score the motifs
        for (int i = 0; i < motifs.size(); i++) {
            String curMotif = motifs.get(i).sequence;
            //iterate through each string in the string[]
            int[] motifScores = new int[sequence.length];
            for (int j = 0; j < sequence.length; j++) {
                if(sequence[j] == null){
                    break;
                }
                //DNA strand sequence
                String curSequence = sequence[j];
                CandidateMotif bestMotif = new CandidateMotif(curMotif, 99999);
                //System.out.println("Current Sequence: " + curSequence);
                //iterate through and compare the curMotif to elements in curSequence
                for (int k = 0; k < curSequence.length() - Integer.parseInt(length); k++) {
                    String curSnippett = curSequence.substring(k, k + Integer.parseInt(length));
                    //System.out.println("(" + curMotif + ", " + curSnippett+ ")" + "-> " + hammingDistance(curMotif, curSnippett));
                    int score = hammingDistance(curMotif, curSnippett);
                    if (score < bestMotif.score) {
                        bestMotif.score = score;
                    }
                    operations++;
                }
                //set the score of the motif for that specific sequence
                motifScores[j] = bestMotif.score;
            }

            //sum all of the scores for this particular motif based on distances of strands
            int motifScore = 0;
            for (int j = 0; j < motifScores.length; j++) {
                motifScore += motifScores[j];
                operations++;
            }
            motifs.get(i).score = motifScore;
            //System.out.println(motifs.get(i));
        }

        //find the minimum score
        CandidateMotif bestMotif = new CandidateMotif("", 9999);
        for (int i = 0; i < motifs.size(); i++) {
            if (motifs.get(i).score < bestMotif.score) {
                bestMotif = motifs.get(i);
            }
            operations++;
        }
        if (isPrintResults) {
            System.out.println("");
            System.out.println("Best Motif: " + bestMotif.sequence);
            System.out.println("");
            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != null) {
                    printExposedMotifs(sequence[i], bestMotif.sequence, Integer.parseInt(length), bestMotif.score);
                    operations++;
                } else {
                    break;
                }
            }
            System.out.println("");
            System.out.println("Total Hemming Distance from Consensus Motif to Best Motif in each Sequence: " + bestMotif.score);
            long endTime = new Date().getTime();
            System.out.println("");
            System.out.println("Elapsed Time: " + (endTime - startTime) + " ms for " + Math.abs(operations) + " fundamental operations");

        }

        return bestMotif;
    }

    //***********************************************************************************************************
    //Method: greedyMotif
    //Description: The second method that searches the DNA sequence via Greedy Motif Search Algorithm. (Faster).
    //Parameters:   String[] sequence -> DNA sequences broken up into strings via array
    //              int length -> The length of the l-mer specified by the user
    //Returns: None
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    
    private static void greedyMotif(String[] sequence, String length) {
        System.out.println("Working...(This could take a while)");
        generateCandidateMotifs(Integer.parseInt(length), "");

        //scans and finds the best consensus motif in first two strands
        ArrayList<CandidateMotif> bestMotifs = new ArrayList<CandidateMotif>();
        //first DNA strand searched via bruteForce Algorithm
        //bestMotifs.add(bruteForce(Arrays.copyOfRange(sequence, 0, 1), length, false));
        //second DNA strand as well
        //bestMotifs.add(bruteForce(Arrays.copyOfRange(sequence,1,2), length, false));
        for(int i = 0; i < Integer.parseInt(length); i++){
            bestMotifs.add(bruteForce(Arrays.copyOfRange(sequence,i,i+1), length, false));
        }
        
        String[] test = new String[Integer.parseInt(length)];
        for(int i = 0; i < Integer.parseInt(length); i++){
            test[i] = bestMotifs.get(i).sequence;
        }
        bruteForce(test, length, true);
    }

    //***********************************************************************************************************
    //Method: generateCandidateMotifs
    //Description: Generates all of the possible motifs that could be used by the algorithm depending on the score
    //Parameters:   String[] sequence -> DNA sequences broken up into strings via array
    //              int length -> The length of the l-mer specified by the user
    //Returns: ArrayList<CandidateMotif> -> collection of all the possible motifs, populated by the end of the 
    //         algorithm
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    
    public static void generateCandidateMotifs(int maxLength, String curr) {
        motifs = new ArrayList<CandidateMotif>();
        
        char[] alphabet = new char[]{'a', 'c', 'g', 't'};

        // If the current string has reached it's maximum length
        if (curr.length() == maxLength) {
            motifs.add(new CandidateMotif(curr, 0));

            // Else add each letter from the alphabet to new strings and process these new strings again
        } else {
            for (int i = 0; i < alphabet.length; i++) {
                String oldCurr = curr;

                curr += alphabet[i];
                generateCandidateMotifs(maxLength, curr);
                curr = oldCurr;
            }
        }
    }

    //***********************************************************************************************************
    //Method: hammingDistance
    //Description: returns the Hamming Distance of two strings
    //Parameters:   String v -> First string that will be compared to the second
    //              String w -> Second string that will be compared to the first
    //Returns: int score -> the Hamming Distance 
    //Globals: None
    //***********************************************************************************************************
    
    public static int hammingDistance(String v, String w) {
        //handles case
        v = v.toLowerCase();
        w = w.toLowerCase();
        
        int distance = 0;
        if (v.length() == w.length()) {
            for (int i = 0; i < v.length(); i++) {
                char vChar = v.charAt(i);
                char wChar = w.charAt(i);
                if (vChar != wChar) {
                    distance++;
                }
            }
        } else {
            System.out.println("Strings weren't of same length.");
            return 0;
        }
        return distance;
    }

    //***********************************************************************************************************
    //Method: printExposedMotifs
    //Description: Prints out the capitalized sequences that match the best motif
    //Parameters:   String sequence -> DNA strand that is currently being printed out
    //              String motif -> motif that needs to be matched with
    //              int length -> The length of the l-mer specified by the user
    //              int maxDistance -> The score of the best motif
    //Returns: None
    //Globals: None
    //***********************************************************************************************************
    
    private static void printExposedMotifs(String sequence, String motif, int length, int maxDistance) {
        int startingPosition = 0;
        int bestScore = 9999;
        //gets the starting position of the best sequence to capitalize
        for (int i = 0; i < sequence.length() - length; i++) {
            int curDistance = hammingDistance(motif, sequence.substring(i, i + length));
            if (curDistance < bestScore) {
                bestScore = curDistance;
                startingPosition = i;
            }
        }
        String printedString = "";

        printedString += sequence.substring(0, startingPosition);
        printedString += sequence.substring(startingPosition, startingPosition + length).toUpperCase();
        printedString += sequence.substring(startingPosition + length, sequence.length());
        System.out.println(printedString);
    }
}
