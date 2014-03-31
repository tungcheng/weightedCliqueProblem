
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.model;

//import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 *
 * @author TungNT
 */
public class Problem {
    private int n;
    private int b;
    private double[][] c;
    private Matrix Q;
    private Matrix q;
    private Matrix oldQ;
    
    public Problem(int n) {
        this.n = n;
        c = new double[n][n];
        double[][] temp = new double[n][1];
        for(int i=0; i<n; i++) {
            temp[i][0] = 0;
        }
        q = new Matrix(temp);
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
    
    public void handleInputMatrix() {
        Q = new Matrix(c);
        oldQ = Q.copy();
        Q.timesEquals(2); // x2 ma tran Q de co he so 1/2 trong 1/2 * xT * Q * x
        Q.print(4, 1);
    }

    public int getN() {
        return n;
    }

    public int getB() {
        return b;
    }

    public Matrix getQ() {
        return Q;
    }

    public Matrix get_q() {
        return q;
    }
    
    public double getValue(double[] x) {
        double[][] temp = new double[x.length][1];
        for(int i=0; i<this.n; i++) {
            temp[i][0] = x[i];
        }
        Matrix mX = new Matrix(temp);

        return getValue(mX);
    }
    
    public double getValue(Matrix mX) {
        Matrix t = mX.transpose().times(this.Q);
        Matrix t2 = t.times(mX);
//        Matrix t3 = getInverse(q).times(mX);
        
        double obj_value;
        obj_value = 0.5*t2.get(0, 0);
        return obj_value;
    }
}
