/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.model;

import Jama.Matrix;
import org.coinor.Ipopt;

/**
 *
 * @author TungNT
 */
public class SubProblem extends Ipopt implements Comparable<SubProblem> {
    
    private int n, m, b, nele_jac, nele_hess;;
    private Matrix mQ, mq;
    
    double[][] temp;
    private Matrix mX;
    
    private int[] indexX;
    private double lowerBound;
    private double[] bestX;
    private double[] preX;
    
    public SubProblem(int n, int m, int b, double[][] Q, double[][] q, int[] indexX, double[] preX) {
        System.out.println("n = " + n);
        System.out.println("m = " + m);
        System.out.println("b = " + b);
        this.n = n;
        this.m = m;
        this.b = b;
        this.mQ = new Matrix(Q);
        this.mq = new Matrix(q);
        this.indexX = indexX;
        this.preX = preX;
        this.lowerBound = Double.MIN_VALUE;
        
        this.mQ.print(4, 1);
        this.mq.print(4, 1);
        
        temp = new double[this.n][1];
        mX = new Matrix(temp);
        
        nele_jac = n;
        nele_hess = n*(n+1)/2;
        
        System.out.println("nele_jac = " + nele_jac);
        System.out.println("nele_hess = " + nele_hess);
        
        double x_L[] = new double[n];
        double x_U[] = new double[n];
        for(int i=0; i < x_L.length; i++){
                x_L[i] = 0.0;
                x_U[i] = 1.0;
        }
        
        double g_L[] = new double[m];
        double g_U[] = new double[m];
        g_U[0] = b;
        g_L[0] = (-1)*2e19;

        int index_style = Ipopt.C_STYLE;
        create(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style);
    }

    public double[] getInitialGuess(){
        /* allocate space for the initial point and set the values */
        double x[] = new double[n];
        System.out.println("n = " + n);
        System.out.println("b = " + b);
//        for(int i=0; i<b; i++) {
//            x[i] = 1.0;
//        }
        int min = (b < n) ? b : n;
        for(int i=0; i<min; i++) {
            x[i] = 1.0;
        }
        for(int i=min; i<n; i++) {
            x[i] = 0.0;
        }

        return x;
    }

    protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value) {
        assert n == this.n;
        
        for(int i=0; i<this.n; i++) {
            temp[i][0] = x[i];
        }

        Matrix t = getInverse(mX).times(mQ);
        Matrix t2 = t.times(mX);
        Matrix t3 = getInverse(mq).times(mX);
        
        obj_value[0] = 0.5 * t2.get(0, 0) + t3.get(0, 0);

        return true;
    }

    protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f) {
        assert n == this.n;

        for(int i=0; i<this.n; i++) {
            temp[i][0] = x[i];
        }
        
        Matrix t = getInverse(mQ).plus(mQ);
        Matrix t2 = t.times(mX);
        Matrix t3 = t2.times(0.5);
        Matrix t4 = mq.plus(t3);
        
        for(int i=0; i<this.n; i++) {
            grad_f[i] = t4.get(i, 0);
        }

        return true;
    }

    protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g) {
            assert n == this.n;
            assert m == this.m;
            
            g[0] = 0;
            for(int i=0; i<n; i++) {
                g[0] += x[i];
            }

            return true;
    }

    protected boolean eval_jac_g(int n, double[] x, boolean new_x,
                    int m, int nele_jac, int[] iRow, int[] jCol, double[] values) {
            assert n == this.n;
            assert m == this.m;

            if (values == null) {
                    for(int i=0; i<n; i++) {
                        iRow[i] = 0;
                        jCol[i] = i;
                    }
            }
            else {
                    for(int i=0; i<n; i++) {
                        values[i] = 1;
                    }
            }

            return true;
    }

    protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor, int m, double[] lambda, boolean new_lambda, int nele_hess, int[] iRow, int[] jCol, double[] values) {
            int idx = 0; /* nonzero element counter */
            int row = 0; /* row counter for loop */
            int col = 0; /* col counter for loop */
            if (values == null) {
                    idx=0;
                    for (row = 0; row < n; row++) {
                            for (col = 0; col <= row; col++) {
                                    iRow[idx] = row;
                                    jCol[idx] = col;
                                    idx++;
                            }
                    }

                    assert idx == nele_hess;
                    assert nele_hess == this.nele_hess;
            }
            else {
                Matrix t = getInverse(mQ).plus(mQ);
                Matrix t2 = t.times(0.5);
                row = 0; col = 0;
                for(int i=0; i<this.nele_hess; i++) {
                    values[i] = obj_factor * t2.get(row, col);
                    col++;
                    if(col > row) {
                        row++;
                        col = 0;
                    }
                }
            }
            return true;
    }    

    public int getN() {
        return n;
    }

    public int getM() {
        return m;
    }

    public int getB() {
        return b;
    }

    public Matrix getmQ() {
        return mQ;
    }

    public Matrix getMq() {
        return mq;
    }

    public int[] getIndexX() {
        return indexX;
    }
    
    private Matrix getInverse(Matrix mX) {
        Matrix mXt = new Matrix( mX.getColumnDimension(), mX.getRowDimension());
        for(int i=0; i<mXt.getRowDimension(); i++) {
            for(int j=0; j<mXt.getColumnDimension(); j++) {
                mXt.set(i, j, mX.get(j, i));
            }
        }
        return mXt;
    }

    public double getLowerBound() {
        return lowerBound;
    }

    public void setLowerBound(double lowerBound) {
        this.lowerBound = lowerBound;
    }

    @Override
    public int compareTo(SubProblem o) {
        if(this.lowerBound < o.getLowerBound()) {
            return -1;
        }
        else if(this.lowerBound > o.getLowerBound()) {
            return 1;
        }
        else {
            return 0;
        }
    }

    public double[] getBestX() {
        return bestX;
    }

    public void setBestX(double[] bestX) {
        this.bestX = bestX;
        if(this.b == 0) {
            for(int i=0; i<this.bestX.length; i++) {
                this.bestX[i] = 0;
            }
        }
    }

    public double[] getPreX() {
        return preX;
    }

}