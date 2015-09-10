/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package motif.finder;

public class CandidateMotif {
    public String sequence;
    public int score;
    
    public CandidateMotif(String sequence, int score){
        this.sequence = sequence;
        this.score = score;
    }
    
    public String toString(){
        return sequence + ", " + score;
    }
}
