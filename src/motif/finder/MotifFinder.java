package motif.finder;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

public class MotifFinder {

    public static String[] sequence;
    public static ArrayList<CandidateMotif> motifs;
    public static void main(String[] args) {
        getUserAction();
    }

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
    
    private static void searchSequence(String[] sequence) {
        KeyboardInputClass input = new KeyboardInputClass();
        String length = input.getKeyboardInput("Specify the length of the motif to be found:");
        System.out.println("");
        String algorithm = input.getKeyboardInput("1. Brute Force Median String Search \n2. Greedy Motif Search");
        
        switch(algorithm){
            case "1":
                bruteForce(sequence, length);
                break;
            case "2":
                greedyMotif(sequence, length);
                break;
        }
    }

    private static void retrieveSequence() throws IOException {
        //test
        File file = new File("Haf.txt");
        // creates the file
        file.createNewFile();
        // creates a FileWriter Object
        FileWriter writer = new FileWriter(file);
        // Writes the content to the file
        writer.write("atcacttg\n"
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
    
    private static void bruteForce(String[] sequence, String length) {
        System.out.println("Working...(This could take a while)");
        getConsensusMotifs(sequence,Integer.parseInt(length));
    }

    private static void greedyMotif(String[] sequence, String length) {
        System.out.println("Working...(This could take a while)");
        getConsensusMotifs(sequence,Integer.parseInt(length));
        
        throw new UnsupportedOperationException("Greedy Motif."); //To change body of generated methods, choose Tools | Templates.
    }
    
    public static ArrayList<CandidateMotif> getConsensusMotifs(String[] sequence, int length){
        //generates all possible motifs with the given length
        motifs = generateSequences(sequence, length);
        
        return motifs;
    }

    private static ArrayList<CandidateMotif> generateSequences(String[] sequence, int length) {
        motifs = new ArrayList<CandidateMotif>();
        generateCandidateMotifs(length, "");
        return motifs;
    }
    
    public static void generateCandidateMotifs(int maxLength, String curr) {
        char[] alphabet = new char[]{'a','c','g','t'};
        
        // If the current string has reached it's maximum length
        if(curr.length() == maxLength) {
            motifs.add(new CandidateMotif(curr,0));

        // Else add each letter from the alphabet to new strings and process these new strings again
        } else {
            for(int i = 0; i < alphabet.length; i++) {
                String oldCurr = curr;
                
                curr += alphabet[i];
                generateCandidateMotifs(maxLength,curr);
                curr = oldCurr;
            }
        }
    }
    
    public static int HammingDistance(String v, String w){
        int distance = 0;
        if(v.length() == w.length()){
            for(int i = 0; i < v.length(); i++){
                if(v.charAt(i) != w.charAt(i)){
                    distance++;
                }
            }
        }else{
            System.out.println("Strings weren't of same length.");
            return 0;
        }
        return distance;
    }
}
