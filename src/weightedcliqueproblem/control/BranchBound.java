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
    
    public BranchBound (Problem prob) {
        mainProb = prob;
        bestRecord = (-1)*Double.MAX_VALUE;
        
        int[] startIndex = new int[prob.getN()];
        curX = new double[prob.getN()];
        for(int i=0; i < curX.length; i++) {
            curX[i] = -1;
            startIndex[i] = i;
            // curX[i] = -1 if x[i] is not integer
            // curX[i] = x[i] if x[i] is integer
        }
        
        this.branchBound(new SubProblem(
                prob.getN(),
                1,
                prob.getB(),
                prob.getQ().getArray(),
                prob.get_q().getArray(),
                startIndex
        ));
        
        showResult();
    }
    
    private Problem mainProb;
    private double bestRecord;
    private Matrix bestX;
    
    private double[] curX;
    
    private double upperBound;
    
    private void branchBound (SubProblem subProb) {
        
        double x[] = subProb.getInitialGuess();
        subProb.solve(x);
        for(int i=0; i<x.length; i++) {
            System.out.println(x[i]);
        }
        
        if(isAllInteger(x)) {
            checkBest(x, subProb.getIndexX());
        }
        else {
            upperBound = calUpperBound(subProb.getObjVal());
            if(upperBound > bestRecord) {
                // branch: make new sub problem
                int subIndex = this.findSubIndex(x);
                double[][] subQ = this.createSubProbMatrixQ(subIndex, subProb.getmQ().getArray());
                double[][] sub_q = this.createSubProbMatrix_q(subIndex, subProb.getMq().getArray());
                int[] subIndexX = this.createSubProbIndexX(subIndex, subProb.getIndexX());
                if(subProb.getN() == 1) {
                    curX[subProb.getIndexX()[subIndex]] = 0;
                    checkBest();
                    if(subProb.getB() >= 1) {
                        curX[subProb.getIndexX()[subIndex]] = 1;
                        checkBest();
                    }
                }
                else {
                    curX[subProb.getIndexX()[subIndex]] = 0;
                    branchBound(new SubProblem(
                            subProb.getN()-1,
                            subProb.getM(),
                            subProb.getB(),
                            subQ,
                            sub_q,
                            subIndexX
                    ));
                    curX[subProb.getIndexX()[subIndex]] = 1;
                    branchBound(new SubProblem(
                            subProb.getN()-1,
                            subProb.getM(),
                            subProb.getB()-1,
                            subQ,
                            sub_q,
                            subIndexX
                    ));
                    curX[subProb.getIndexX()[subIndex]] = -1;
                }
            }
            else {
                // backtrack
            }
        }
    }
    
    private int findSubIndex(double[] x) {
        int subIndex = 0;
        double subTemp = x[subIndex];
        double temp;
        for(int i=0; i < x.length; i++) {
            temp = Math.min(x[i]-0, 1-x[i]);
            if(temp > subTemp) {
                subIndex = i;
                subTemp = x[subIndex];
            }
        }
        return subIndex;
    }
    
    private double[][] createSubProbMatrixQ (int subIndex, double[][] preQ) {
        int n = preQ.length;
        double[][] Q = new double[n-1][n-1];
        int i, j;
        i = 0; j = 0;
        for(int pi=0; (pi<n)&&(pi!=subIndex); pi++) {
            for(int pj=0; (pj<n)&&(pj!=subIndex); pj++) {
                Q[i][j] = preQ[pi][pj];
                j++;
            }
            i++;
        }
        return Q;
    }
    
    private double[][] createSubProbMatrix_q (int subIndex, double[][] pre_q) {
        int n = pre_q.length;
        double[][] q = new double[n-1][1];
        int i = 0;
        for(int pi=0; (pi<n)&&(pi!=subIndex); pi++) {
            q[i][0] = pre_q[pi][0];
            i++;
        }
        return q;
    }
    
    private int[] createSubProbIndexX (int subIndex, int[] preX) {
        int n = preX.length;
        int[] indexX = new int[n-1];
        int i = 0;
        for(int pi=0; (pi<n)&&(pi!=subIndex); pi++) {
            indexX[i] = preX[pi];
            i++;
        }
        return indexX;
    }
    
    private boolean isAllInteger (double[] x) {
        for(int i=0; i < x.length; i++) {
            if(x[i] != Math.round(x[i])) {
                return false;
            }
        }
        return true;
    }
    
    private void checkBest (double[] x, int[] indexX) {
        for(int i=0; i < x.length; i++) {
            this.curX[indexX[i]] = x[i];
        }
        checkBest();
    }

    private void checkBest() {
        double[][] temp = new double[curX.length][1];
        for(int i=0; i < curX.length; i++) {
            temp[i][0] = curX[i];
        }
        double curCheck = this.mainProb.getValue(this.curX);
        if(curCheck > this.bestRecord) {
            this.bestRecord = curCheck;
            this.bestX = new Matrix(temp);
        }
    }

    private double calUpperBound(double objVal) {
        double upperBound = 0;
        for(int i=0; i < curX.length; i++) {
            if(curX[i] != -1) {
                upperBound += curX[i];
            }
        }
        upperBound += objVal;
        return upperBound;
    }

    private void showResult() {
        System.out.println("best X[] :");
        for(int i=0; i<curX.length; i++) {
            System.out.print(" " + curX[i]);
        }
        System.out.println("");
        System.out.println("best obj_val = " + this.mainProb.getValue(curX));
    }
    
}
