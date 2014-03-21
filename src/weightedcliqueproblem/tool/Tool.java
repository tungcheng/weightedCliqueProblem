/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.tool;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import javax.swing.JTextArea;

/**
 *
 * @author TungNT
 */
public class Tool {
    private final String newline = "\n";
    
    private static Tool instance;
    
    private JTextArea log;
    NumberFormat fm;
    
    public static Tool getTool() {
        if(instance == null) {
            instance = new Tool();
        }
        return instance;
    }
    
    public void setWriter(JTextArea log) {
        this.log = log;
        fm = new DecimalFormat("#0.0000");
    }
    
    public boolean isAllInteger (double[] x) {
        for(int i=0; i < x.length; i++) {
            if(x[i] != Math.round(x[i])) {
                return false;
            }
        }
        return true;
    }
    
    public long getTimeNow() {
        return System.currentTimeMillis();
    }
    
    public void show(String info, double[] x) {
        this.println(info);
        for(int i=0; i<x.length-1; i++) {
            this.print(" " + fm.format(x[i]) + ";");
        }
        this.print(" " + fm.format(x[x.length-1]) + "");
        this.println();
    }
    
    public void show(String info) {
        this.println(info);
    }
    
    public void show(String info, double x) {
        this.print(info);
        this.println(fm.format(x));
    }
    
    public void show(String info, int x) {
        this.print(info);
        this.println(x+"");
    }
    
    public void show(String info, long x) {
        this.print(info);
        this.println(x+"");
    }
    
    public void showLine() {
        this.println("--------------");
    }
    
    public void show(String info, int[] x) {
        this.println(info);
        for(int i=0; i<x.length-1; i++) {
            this.print(" " + x[i] + ";");
        }
        this.print(" " + x[x.length-1] + "");
        this.println();
    }
    
    private void print(String s) {
        this.log.append(s);
    }
    
    private void println(String s) {
        this.log.append(s);
        this.log.append(newline);
    }
    
    private void println() {
        this.log.append(newline);
    }
}
