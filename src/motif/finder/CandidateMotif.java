//***********************************************************************************************************
//***********************************************************************************************************

//Class:    CandidateMotif
//Description: Object that is used to represent a motif and its score.  The score can be either the Hamming Distance or it Alignment Score
//             depending on what is needed.

//***********************************************************************************************************
//***********************************************************************************************************
package motif.finder;

public class CandidateMotif {
    public String sequence; //sequence/pattern for the motif
    public int score; //hamming distance or score this particular motif scored
    
    public CandidateMotif(String sequence, int score){
        this.sequence = sequence;
        this.score = score;
    }
    
    public String toString(){
        return sequence + ", Score:" + score;
    }
}
