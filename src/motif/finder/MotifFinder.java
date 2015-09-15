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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

//Method: Main
//Description: Prompts the user for the information required to run the algorithm and runs the algorithm
//Returns: None
public class MotifFinder {

    public static String[] sequence;    //DNA strands that will be analyzed later in the program
    public static ArrayList<CandidateMotif> motifs; //all of the possible different candidate motifs
    public static int finalScore;   //final Hamming Distance score to be printed out
    public static int operations;   //number of operations each algorithm takes

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
        String decision = input.getKeyboardInput("");

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
        //user hasn't specified text file
        if (sequence == null) {
            System.out.println("ERROR: Please import data to be traversed.");
            getUserAction();
        }

        //gets the users decision
        KeyboardInputClass input = new KeyboardInputClass();
        int length = input.getInteger(true, 0, 0, 9999, "Specify the length of the motif to be found:");
        System.out.println("");
        String algorithm = input.getKeyboardInput("1. Brute Force Median String Search \n2. Greedy Motif Search");

        if (length > sequence[0].length()) {
            System.out.println("The length cannot be greater than the number of elements in each DNA strand.  Please try again.");
            System.out.println("");
            searchSequence(sequence);
        }

        //Creates option tree
        switch (algorithm) {
            case "1":
                bruteForce(sequence, length, true);
                break;
            case "2":
                greedyMotif(sequence, length, true);
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

        TextFileClass textFile = new TextFileClass();
        textFile.getFileName("Specify the text file to be read:");
        textFile.getFileContents();
        if (textFile.fileName.length() > 0 && textFile.lineCount != 0) {
            sequence = textFile.text;

            //Check to see if all strands are the same length
            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != null) {
                    if (sequence[i].length() != sequence[0].length()) {
                        sequence = null;
                        System.out.println("ERROR: DNA strands are not the same length. Please try again.");
                        retrieveSequence();
                    }
                }
            }

            //Checks for invalid characters
            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != null) {
                    for (int j = 0; j < sequence[0].length(); j++) {
                        char cur = sequence[i].charAt(j);
                        if (cur != 'a' && cur != 'g' && cur != 'c' && cur != 't') {
                            sequence = null;
                            System.out.println("ERROR: DNA strands contain invalid characters.  Please try again.");
                            retrieveSequence();
                        }

                    }
                }
            }
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
    private static void bruteForce(String[] sequence, int length, boolean isPrintResults) {
        if (isPrintResults) {
            System.out.println("Working...(This could take a while)");
        }

        //initializes the candidate motifs with scores of 0
        long startTime = new Date().getTime();
        getConsensusMotifs(sequence, length);

        //iterate through the sequences and score the motifs
        for (int i = 0; i < motifs.size(); i++) {
            String curMotif = motifs.get(i).sequence;
            //iterate through each string in the string[]
            int[] motifScores = new int[sequence.length];
            for (int j = 0; j < sequence.length; j++) {
                if (sequence[j] == null) {
                    break;
                }
                //DNA strand sequence
                String curSequence = sequence[j];
                CandidateMotif bestMotif = new CandidateMotif(curMotif, 99999);
                //System.out.println(curSequence);
                //iterate through and compare the curMotif to elements in curSequence
                for (int k = 0; k < curSequence.length() - length; k++) {
                    String curSnippett = curSequence.substring(k, k + length);
                    //System.out.println("(" + curMotif + ", " + curSnippett+ ")" + "-> " + hammingDistance(curMotif, curSnippett));
                    int score = hammingDistance(curMotif, curSnippett);
                    if (score < bestMotif.score) {
                        bestMotif.score = score;
                    }
                    operations++;
                }
                //set the score of the motif for that specific sequence
                motifScores[j] = bestMotif.score;
                operations++;
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

        //Prints the results of the search
        if (isPrintResults) {
            System.out.println("");
            System.out.println("Best Motif: " + bestMotif.sequence);
            System.out.println("");
            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != null) {
                    printExposedMotifs(sequence[i], bestMotif.sequence, length, bestMotif.score);
                    operations++;
                } else {
                    break;
                }
            }
            System.out.println("");
            System.out.println("Total Hemming Distance from Consensus Motif to Best Motif in each Sequence: " + Math.abs(bestMotif.score));
            long endTime = new Date().getTime();
            System.out.println("");
            System.out.println("Elapsed Time: " + (endTime - startTime) + " ms for " + operations + " fundamental operations");

        }

        //run the program again
        getUserAction();
    }

    //***********************************************************************************************************
    //Method: greedyMotif
    //Description: The second method that searches the DNA sequence via Greedy Motif Search Algorithm. (Faster).
    //Parameters:   String[] sequence -> DNA sequences broken up into strings via array
    //              int length -> The length of the l-mer specified by the user
    //Returns: None
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    private static void greedyMotif(String[] sequence, int length, boolean isPrintResults) {
        System.out.println("Working...(This could take a while)");

        long startTime = new Date().getTime();

        //scans and finds the best consensus motif in first two strands
        ArrayList<CandidateMotif> bestMotifs = new ArrayList<CandidateMotif>();
        CandidateMotif s1 = new CandidateMotif("", 0);
        CandidateMotif s2 = new CandidateMotif("", 0);

        //brute forces the first two lines of the dna strand
        for (int i = 0; i < sequence[0].length() - length; i++) {
            String strand1 = sequence[0].substring(i, i + length);
            for (int j = 0; j < sequence[1].length() - length; j++) {
                String strand2 = sequence[1].substring(j, j + length);
                String[] strands = {strand1, strand2};
                CandidateMotif c = score(strands);
                if (c.score > s1.score) {
                    s1 = new CandidateMotif(strand1, c.score);
                    s2 = new CandidateMotif(strand2, c.score);
                }
                operations++;
            }
            operations++;
        }

        //add the best motifs to the list
        bestMotifs.add(s1);
        bestMotifs.add(s2);

        for (int i = 2; i < sequence.length; i++) {
            if (sequence[i] != null) {
                CandidateMotif bestMotif = new CandidateMotif("", 0);
                for (int j = 0; j < sequence[0].length() - length; j++) {
                    String strand = sequence[i].substring(j, j + length);
                    String[] curStrands = new String[bestMotifs.size() + 1];
                    for (int k = 0; k < bestMotifs.size(); k++) {
                        curStrands[k] = bestMotifs.get(k).sequence;
                    }
                    curStrands[bestMotifs.size()] = strand;
                    //System.out.println(Arrays.toString(curStrands));
                    CandidateMotif c = score(curStrands);
                    if (c.score > bestMotif.score) {
                        bestMotif = c;
                    }
                    operations++;
                }
                operations++;
                bestMotifs.add(bestMotif);
            }
        }

        //compiles all the sequences and scores the final one to get the best candidate
        String[] finalMotifs = new String[bestMotifs.size()];
        for (int i = 0; i < finalMotifs.length; i++) {
            String cur = bestMotifs.get(i).sequence;
            finalMotifs[i] = cur;
            operations++;
        }

        //the bestMotif with is score
        CandidateMotif theBestMotif = score(finalMotifs);

        //Prints the results of the search
        if (isPrintResults) {
            System.out.println("");
            System.out.println("Best Motif: " + theBestMotif.sequence);
            System.out.println("");
            int finalScore = 0;
            for (int i = 0; i < sequence.length; i++) {
                if (sequence[i] != null) {
                    finalScore += printExposedMotifs(sequence[i], theBestMotif.sequence, length, theBestMotif.score);
                    operations++;
                } else {
                    break;
                }
            }
            System.out.println("");
            System.out.println("Total Hemming Distance from Consensus Motif to Best Motif in each Sequence: " + Math.abs(finalScore));
            long endTime = new Date().getTime();
            System.out.println("");
            System.out.println("Elapsed Time: " + (endTime - startTime) + " ms for " + operations + " fundamental operations");
        }
        //run the program again
        getUserAction();
    }

    //***********************************************************************************************************
    //Method: score
    //Description: Takes an array of strings and scores them based on the frequencies of the letters
    //Parameters:   String[] strands -> DNA sequences broken up into strings via array
    //Returns: CandidateMotif bestMotif -> returns the best consensus motif along with the score it had
    //Globals: None
    //***********************************************************************************************************
    private static CandidateMotif score(String[] strands) {
        CandidateMotif bestMotif = new CandidateMotif("", 0);
        for (int i = 0; i < strands[0].length(); i++) {
            //generate the column 
            String column = "";
            for (int j = 0; j < strands.length; j++) {
                column += strands[j].charAt(i);
            }

            int[] letterCount = new int[4];

            //count the frequencies of the letters
            for (int j = 0; j < column.length(); j++) {
                char cur = column.charAt(j);
                if (cur == 'a') {
                    letterCount[0]++;
                } else if (cur == 'c') {
                    letterCount[1]++;
                } else if (cur == 'g') {
                    letterCount[2]++;
                } else if (cur == 't') {
                    letterCount[3]++;
                }
                operations++;
            }

            //finds the max score using scoreIndex [0] = index [1] = count
            int[] scoreIndex = new int[]{0, 0};
            for (int j = 0; j < letterCount.length; j++) {
                if (letterCount[j] > scoreIndex[1]) {
                    scoreIndex[1] = letterCount[j];
                    scoreIndex[0] = j;
                }
                operations++;
            }
            operations++;

            //adds to the sequence returned
            if (scoreIndex[0] == 0) {
                bestMotif.sequence += 'a';
            } else if (scoreIndex[0] == 1) {
                bestMotif.sequence += 'c';
            } else if (scoreIndex[0] == 2) {
                bestMotif.sequence += 'g';
            } else if (scoreIndex[0] == 3) {
                bestMotif.sequence += 't';
            }
            bestMotif.score += scoreIndex[1];

        }
        return bestMotif;
    }

    //***********************************************************************************************************
    //Method: getConsensusMotifs
    //Description: Initializes the motifs and calls the method that generates them
    //Parameters:   String[] sequence -> DNA sequences broken up into strings via array
    //              int length -> The length of the l-mer specified by the user
    //Returns: ArrayList<CandidateMotif> motifs-> collection of all the possible motifs, populated by the end of the 
    //         algorithm generated by generateCandidateMotifs()
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    public static ArrayList<CandidateMotif> getConsensusMotifs(String[] sequence, int length) {
        //generates all possible motifs with the given length
        motifs = new ArrayList<CandidateMotif>();
        generateCandidateMotifs(length, "");
        return motifs;
    }

    //***********************************************************************************************************
    //Method: generateCandidateMotifs
    //Description: Generates all of the possible motifs that could be used by the algorithm depending on the score.
    //             Received much help for developing this method on StackOverflow.
    //Parameters:   int length -> The length of the l-mer specified by the user
    //              String current -> Most important piece of recursion call.  It keeps on adding onto it until
    //              it reaches the length specified by the other parameter.  
    //Returns: ArrayList<CandidateMotif> -> collection of all the possible motifs, populated by the end of the 
    //         algorithm
    //Globals: motifs -> Collection of all possible motifs
    //***********************************************************************************************************
    public static void generateCandidateMotifs(int length, String current) {
        //letters we want to populate the string with
        char[] alphabet = new char[]{'a', 'c', 'g', 't'};

        //String has reached the ideal length, add it to our arraylist
        if (current.length() == length) {
            motifs.add(new CandidateMotif(current, 0));

            //string isn't the length we need, continue down the recursive calls and for loops
        } else {
            //we need to do this 4 times (will be called moreso because of the recursive calls)
            for (int i = 0; i < alphabet.length; i++) {
                //hold onto the old string
                String old = current;

                //increment the current string with the alphabet and see if it is the right length
                current += alphabet[i];
                //recursive call down deeper into the algorithm
                generateCandidateMotifs(length, current);

                //current becomes the old because we've come up out of the recursive calls and need to restart the process
                current = old;
            }
        }
    }

    //***********************************************************************************************************
    //Method: hammingDistance
    //Description: returns the Hamming Distance of two strings
    //Parameters:   String v -> First string that will be compared to the second
    //              String w -> Second string that will be compared to the first
    //Returns: int distance -> the Hamming Distance 
    //Globals: None
    //***********************************************************************************************************
    public static int hammingDistance(String v, String w) {
        int distance = 0;
        if (v.length() == w.length()) {
            for (int i = 0; i < v.length(); i++) {
                char vChar = v.charAt(i);
                char wChar = w.charAt(i);
                if (vChar != wChar) {
                    distance++;
                }
                operations++;
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
    //Returns: bestScore -> returns the best score for that specific sequence to be used in the final hamming distance
    //Globals: None
    //***********************************************************************************************************
    private static int printExposedMotifs(String sequence, String motif, int length, int maxDistance) {
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

        //prints the beginning of the string
        if (startingPosition != 0) {
            printedString += sequence.substring(0, startingPosition);
        }
        //prints the capital letters
        printedString += sequence.substring(startingPosition, startingPosition + length).toUpperCase();
        //prints the rest of the string
        printedString += sequence.substring(startingPosition + length, sequence.length());
        System.out.println(printedString);
        operations++;

        return bestScore;
    }
}
