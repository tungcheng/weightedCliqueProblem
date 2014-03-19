/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.model;

import java.util.ArrayList;

/**
 *
 * @author TungNT
 */
public class ProblemSet extends ArrayList<SubProblem> {
    
    public SubProblem getMin() {
        if(!this.isEmpty()) {
            SubProblem minProb = this.get(0);
            for(int i=0; i < this.size(); i++) {
                if(this.get(i).compareTo(minProb) < 0) {
                    minProb = this.get(i);
                }
            }
            return minProb;
        }
        return null;
    }
    
}
