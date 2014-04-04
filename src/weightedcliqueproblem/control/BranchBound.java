/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.control;

import Jama.Matrix;
import weightedcliqueproblem.model.Problem;
import weightedcliqueproblem.model.ProblemSet;
import weightedcliqueproblem.model.SubDCA;
import weightedcliqueproblem.model.SubProblem;
import weightedcliqueproblem.tool.Tool;

/**
 *
 * @author TungNT
 */
public class BranchBound {
    
    private static final int MAX_ITERVAL = 500;
    private int countIter;
    private long startTime;
    private long endTime;
    
    private final Problem mainProb;
    private double upperBound;
    private double[] bestX;
    private ProblemSet probSet;
    
    public BranchBound (Problem prob) {
        mainProb = prob;
        countIter = 0;
        startTime = Tool.getTool().getTimeNow();
        
        bestX = new double[prob.getN()];
        
        this.init();
        
        this.branchBound();
    }
    
    private void init() {
        Tool.getTool().showLine();
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
            mainProb.getQ(),
            mainProb.get_q(),
            0,
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
        while(countIter < BranchBound.MAX_ITERVAL) {
            countIter++;
            if(countIter % 50 == 0) {
                Runtime.getRuntime().gc();
            }
            if(countIter == BranchBound.MAX_ITERVAL) {
                Tool.getTool().show("Max iterval reach: " + BranchBound.MAX_ITERVAL +" itervals");
                break;
            }
            
            Tool.getTool().showLine();
            Tool.getTool().show("Iterval ID: ", countIter);
            SubProblem subProb = probSet.getMin();
            int subIndex = this.findSubIndex(subProb.getBestX());
            
            Matrix subQ = this.createSubProbMatrixQ(subIndex, subProb.getmQ());
            Matrix sub_q = this.createSubProbMatrix_q(subIndex, subProb.getMq());
            Matrix sub_q1 = this.createSubProbMatrix_q1(sub_q, subProb.getmQ(), subIndex);
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
                    subProb.getD(),
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
                    sub_q1,
                    subProb.getD() + 0.5*subProb.getmQ().get(subIndex, subIndex),
                    subIndexX,
                    preX_1
                );
                solveSubProblem(subProb_1);
                
                Tool.getTool().show("sub 0 lowerBoud: ", subProb_0.getLowerBound());
                Tool.getTool().show("sub 1 lowerBoud: ", subProb_1.getLowerBound());
                if(subProb_0.getLowerBound() < this.upperBound) {
                    Tool.getTool().show("add subProb_0 to problemSet");
                    this.probSet.add(subProb_0);
                }
                if(subProb_1.getLowerBound() < this.upperBound) {
                    Tool.getTool().show("add subProb_1 to problemSet");
                    this.probSet.add(subProb_1);
                }
            }
            
            this.probSet.remove(subProb);
            
            if(this.probSet.isEmpty() 
                    || this.upperBound <= this.probSet.getMin().getLowerBound() ) {
                
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
        calUpperBound(subProb);
        solveDCA(subProb);
        
        if(Tool.getTool().isAllInteger(x) && subProb.getValue(x) <= this.upperBound) {
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
//        this.upperBound = this.mainProb.getValue(this.bestX);
        this.upperBound = subProb.getLowerBound();
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
    
    private Matrix createSubProbMatrixQ (int subIndex, Matrix preQ) {
        int n = preQ.getRowDimension();
        Matrix Q = new Matrix(n-1, n-1);
        int i, j;
        i = 0; j = 0;
        for(int pi=0; pi<n; pi++) {
            if(pi != subIndex) {
                for(int pj=0; pj<n; pj++) {
                    if(pj != subIndex) {
                        Q.set(i, j, preQ.get(pi, pj));
                        j++;
                        if(j == Q.getRowDimension()) {
                            j = 0;
                        }
                    }
                }
                i++;
                if(i == Q.getRowDimension()) {
                    i = 0;
                }
            }
        }
        return Q;
    }
    
    private Matrix createSubProbMatrix_q (int subIndex, Matrix pre_q) {
        int n = pre_q.getRowDimension();
        Matrix q = new Matrix(n-1, 1);
        int i = 0;
        for(int pi=0; pi<n; pi++) {
            if(pi != subIndex) {
                q.set(i, 0, pre_q.get(pi, 0));
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
        double lowerBound = subProb.getObjVal();
        subProb.setLowerBound(lowerBound);
    }

    private void showResult() {
        endTime = Tool.getTool().getTimeNow();
        Tool.getTool().showLine();
        long solveTime = this.endTime - this.startTime;
        Tool.getTool().show("Solve time (ms): ", solveTime);
        Tool.getTool().show("Number of interval: ", this.countIter);
        Tool.getTool().show("Last bestX: ", bestX);
        Tool.getTool().show("Best value: ", this.mainProb.getValue(bestX));
        Tool.getTool().showLine();
    }

    private Matrix createSubProbMatrix_q1(Matrix sub_q, Matrix mQ, int subIndex) {
        int n = mQ.getRowDimension();
        Matrix col = new Matrix(n - 1, 1);
        
        int ti = 0;
        for(int i=0; i<n; i++) {
            if(i != subIndex) {
                col.set(ti, 0, mQ.get(i, subIndex));
                ti++;
            }
        }
        Matrix q1 = sub_q.plus(col);
        return q1;
    }

    private void calUpperBound(SubProblem subProb) {
        double[] x = roundX(subProb.getBestX(), subProb.getB());
        double value = subProb.getValue(x);
        if(value < this.upperBound) {
            this.upperBound = subProb.getValue(x);
            Tool.getTool().show("Upper X: ", x);
            Tool.getTool().show("New upperBound: ", upperBound);
        }
    }

    private double[] roundX(double[] doubleX, int b) {
        double[] intX = new double[doubleX.length];
        int count = 0;
        for(int i=0; i<doubleX.length; i++) {
            if(count < b) {
                if(doubleX[i] < 0.5) {
                    intX[i] = 0;
                }
                else {
                    intX[i] = 1;
                    count++;
                }
            }
            else {
                intX[i] = 0;
            }
        }
        return intX;
    }

    private void solveDCA(SubProblem subProb) {
        double[] x = subProb.getBestX();
        double[] y = new double[x.length];
        double ro = subProb.getMinEig();
        for(int i=0; i<x.length; i++) {
            y[i] = (x[i] >= 0.5) ? (1 - 0.5*ro) : (-1 - 0.5*ro);
        }
        SubDCA subDCA = new SubDCA(
                subProb.getN(),
                subProb.getM(),
                subProb.getB(),
                subProb.getmQ(),
                subProb.getMq(),
                ro,
                y
        );
        double[] x2 = subDCA.getInitialGuess();
        subDCA.solve(x2);
        Tool.getTool().show("DCA x: ", x2);
        Tool.getTool().show("DCA objF: ", subDCA.getObjVal());
    }
    
}
