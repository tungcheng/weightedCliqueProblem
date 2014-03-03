
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.model;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 *
 * @author TungNT
 */
public class Problem {
    int n;
    int b;
    double[][] c;
    Matrix Q;
    Matrix q;
    Matrix I, e;
    
    public Problem(int n) {
        this.n = n;
        c = new double[n][n];
        double[][] temp = new double[n][1];
        for(int i=0; i<n; i++) {
            temp[i][0] = 0;
        }
        q = new Matrix(temp);
        
        temp = new double[n][1];
        for(int i=0; i<n; i++) {
            temp[i][0] = 1;
        }
        e = new Matrix(temp);
        
        temp = new double[n][n];
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                temp[i][j] = (i == j) ? 1 : 0;
            }
        }
        I = new Matrix(temp);
    }
    
    public void setB(int b) {
        this.b = b;
    }
    
    public void setC(double x, int i, int j) {
        this.c[i][j] = x;
    }
    
    public boolean isSymmetric() {
        for(int i=0; i<n; i++) {
            for(int j=i; j<n; j++) {
                if(c[i][j] != c[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }
    
    public void makePositiveSemidefinite() {
        Q = new Matrix(c);
        Q.timesEquals(2); // x2 ma tran Q de co he so 1/2 trong 1/2 * xT * Q * x
        Q.print(4, 1);
        double minEigValue = this.getSmallestEigenvalue();
        if(minEigValue < 0) {
            //double muy = Math.floor(minEigValue);
            double muy = minEigValue - 0.1;
            System.out.println("min eigv: " + minEigValue);
            System.out.println("muy: " + muy);
            Matrix temp = I.times(muy);
            Q.minusEquals(temp);
            
            muy *= 0.5;
            temp = e.times(muy);
            q.plusEquals(temp);
            
            System.out.println("new Q: ");
            Q.print(4, 1);
            System.out.println("new q: ");
            q.print(4, 1);
            minEigValue = this.getSmallestEigenvalue();
            System.out.println("min eigv: " + minEigValue);
        }
    }
    
    private double getSmallestEigenvalue() {
        EigenvalueDecomposition ed = Q.eig();
        Matrix x = ed.getD();
        double min = x.get(0, 0);
        for(int i=0; i<n; i++) {
            if(min > x.get(i, i)) {
                min = x.get(i, i);
            }
        }
        return min;
    }
    
}
