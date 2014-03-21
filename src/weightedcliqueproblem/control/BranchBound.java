/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.control;

import weightedcliqueproblem.model.Problem;
import weightedcliqueproblem.model.ProblemSet;
import weightedcliqueproblem.model.SubProblem;
import weightedcliqueproblem.tool.Tool;

/**
 *
 * @author TungNT
 */
public class BranchBound {
    
    private Problem mainProb;
    private double upperBound;
    private double[] bestX;
    private ProblemSet probSet;
    
    public BranchBound (Problem prob) {
        mainProb = prob;
        
        bestX = new double[prob.getN()];
        
        this.init();
        
        this.branchBound();
    }
    
    private void init() {
        upperBound = Double.MAX_VALUE;
        probSet = new ProblemSet();
        
        int[] startIndex = new int[mainProb.getN()];
        double[] preX = new double[mainProb.getN()];
        for(int i=0; i < startIndex.length; i++) {
            startIndex[i] = i;
            preX[i] = -1;
        }
        
        SubProblem startProb = new SubProblem(
            mainProb.getN(),
            1,
            mainProb.getB(),
            mainProb.getQ().getArray(),
            mainProb.get_q().getArray(),
            startIndex,
            preX
        );
        
        double[] x = startProb.getInitialGuess();
        startProb.solve(x);
        startProb.setBestX(x);
        Tool.getTool().show("Relax bestX: ", x);
        
        if(Tool.getTool().isAllInteger(x)) {
            // ra ket qua luon
        }
        
        probSet.add(startProb);
    }

    private void branchBound() {
        while(true) {
            SubProblem subProb = probSet.getMin();
            int subIndex = this.findSubIndex(subProb.getBestX());
            
            double[][] subQ = this.createSubProbMatrixQ(subIndex, subProb.getmQ().getArray());
            double[][] sub_q = this.createSubProbMatrix_q(subIndex, subProb.getMq().getArray());
            int[] subIndexX = this.createSubProbIndexX(subIndex, subProb.getIndexX());
            
            if(subProb.getN() > 1) {
                double[] preX_0 = subProb.getPreX().clone();
                preX_0[subProb.getIndexX()[subIndex]] = 0;
                SubProblem subProb_0 = new SubProblem(
                    subProb.getN()-1,
                    subProb.getM(),
                    subProb.getB(),
                    subQ,
                    sub_q,
                    subIndexX,
                    preX_0
                );
                solveSubProblem(subProb_0);

                double[] preX_1 = subProb.getPreX().clone();
                preX_1[subProb.getIndexX()[subIndex]] = 1;
                SubProblem subProb_1 = new SubProblem(
                    subProb.getN()-1,
                    subProb.getM(),
                    subProb.getB()-1,
                    subQ,
                    sub_q,
                    subIndexX,
                    preX_1
                );
                solveSubProblem(subProb_1);
                
                Tool.getTool().show("sub 0 lowerBoud: ", subProb_0.getLowerBound());
                Tool.getTool().show("sub 1 lowerBoud: ", subProb_1.getLowerBound());
                if(subProb_0.getLowerBound() < this.upperBound) {
                    this.probSet.add(subProb_0);
                }
                if(subProb_1.getLowerBound() < this.upperBound) {
                    this.probSet.add(subProb_1);
                }
            }
            
            this.probSet.remove(subProb);
            
            if(this.probSet.isEmpty() 
                    || this.upperBound == this.probSet.getMin().getLowerBound() ) {
                
                this.showResult();
                break;
            }
        }
    }
    
    private void solveSubProblem(SubProblem subProb) {
        double[] x = subProb.getInitialGuess();
        subProb.solve(x);
        subProb.setBestX(x);
        Tool.getTool().show("Relax bestX: ", x);
        calLowerBound(subProb);
        
        if(Tool.getTool().isAllInteger(x)) {
            this.updateBest(subProb);
        }
    }
    
    private void updateBest(SubProblem subProb) {
        double[] preX = subProb.getPreX();
        double[] newX = subProb.getBestX();
        int[] indexX = subProb.getIndexX();
        this.bestX = preX.clone();
        for(int i=0; i<indexX.length; i++) {
            this.bestX[indexX[i]] = newX[i];
        }
        this.upperBound = this.mainProb.getValue(this.bestX);
        Tool.getTool().show("new bestX: ", this.bestX);
        Tool.getTool().show("new upperBound: ", upperBound);
    }
    
    private int findSubIndex(double[] x) {
        int subIndex = 0;
        double subTemp = 0;
        double temp;
        for(int i=0; i < x.length; i++) {
            temp = Math.min(x[i]-0, 1-x[i]);
            if(temp > subTemp) {
                subIndex = i;
                subTemp = temp;
            }
        }
        Tool.getTool().show("SubIndex: ", subIndex);
        return subIndex;
    }
    
    private double[][] createSubProbMatrixQ (int subIndex, double[][] preQ) {
        int n = preQ.length;
        double[][] Q = new double[n-1][n-1];
        int i, j;
        i = 0; j = 0;
        for(int pi=0; pi<n; pi++) {
            if(pi != subIndex) {
                for(int pj=0; pj<n; pj++) {
                    if(pj != subIndex) {
//                        System.out.println("i, j: " + i + ", " + j);
//                        System.out.println("pi, pj: " + pi + ", " + pj);
                        Q[i][j] = preQ[pi][pj];
                        j++;
                        if(j == Q.length) {
                            j = 0;
                        }
                    }
                }
                i++;
                if(i == Q.length) {
                    i = 0;
                }
            }
        }
        return Q;
    }
    
    private double[][] createSubProbMatrix_q (int subIndex, double[][] pre_q) {
        int n = pre_q.length;
        double[][] q = new double[n-1][1];
        int i = 0;
        for(int pi=0; pi<n; pi++) {
            if(pi != subIndex) {
                q[i][0] = pre_q[pi][0];
                i++;
            }
        }
        return q;
    }
    
    private int[] createSubProbIndexX (int subIndex, int[] preX) {
        int n = preX.length;
        int[] indexX = new int[n-1];
        int i = 0;
        for(int pi=0; pi<n; pi++) {
            if(pi != subIndex) {
                indexX[i] = preX[pi];
                i++;
            }
        }
        return indexX;
    }

    private void calLowerBound(SubProblem subProb) {
        double[] preX = subProb.getPreX();
        double[] newX = subProb.getBestX();
        int[] indexX = subProb.getIndexX();
        double[] tempX = preX.clone();
        for(int i=0; i<indexX.length; i++) {
            tempX[indexX[i]] = newX[i];
        }
        Tool.getTool().show("sub prob X: ", tempX);
        double lowerBound = this.mainProb.getValue(tempX);
        subProb.setLowerBound(lowerBound);
    }

    private void showResult() {
        Tool.getTool().show("Last bestX: ", bestX);
        Tool.getTool().show("Best value: ", this.mainProb.getValue(bestX));
    }
    
}
