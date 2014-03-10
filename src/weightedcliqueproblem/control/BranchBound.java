/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.control;

import Jama.Matrix;
import weightedcliqueproblem.model.Problem;
import weightedcliqueproblem.model.SubProblem;

/**
 *
 * @author TungNT
 */
public class BranchBound {
    
    public BranchBound(Problem prob) {
        this.bestRecord = Double.MAX_VALUE;
    }
    
    private double bestRecord;
    private Matrix bestX;
    
    private double upperBound;
    
    private void branchBound (SubProblem subProb) {
        
        double x[] = subProb.getInitialGuess();
        subProb.solve(x);
        
        if(isAllInteger(x)) {
            updateBest(x);
        }
        else {
            if(upperBound > bestRecord) {
                // branch: make new sub problem
            }
            else {
                // backtrack
            }
        }
    }
    
    private boolean isAllInteger (double[] x) {
        for(int i=0; i < x.length; i++) {
            if(x[i] == Math.round(x[i])) {
                return false;
            }
        }
        return true;
    }
    
    private void updateBest(double[] x) {
    }
    
}
