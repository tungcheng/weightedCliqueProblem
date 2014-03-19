/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.tool;

/**
 *
 * @author TungNT
 */
public class Tool {
    private static Tool instance;
    
    public static Tool getTool() {
        if(instance == null) {
            instance = new Tool();
        }
        return instance;
    }
    
    public boolean isAllInteger (double[] x) {
        for(int i=0; i < x.length; i++) {
            if(x[i] != Math.round(x[i])) {
                return false;
            }
        }
        return true;
    }
    
    public void show(String info, double[] x) {
        System.out.println(info);
        for(int i=0; i<x.length; i++) {
            System.out.print(" " + x[i]);
        }
        System.out.println("");
    }
    
    public void show(String info, double x) {
        System.out.println(info);
        System.out.println(x);
    }
    
    public void show(String info, int[] x) {
        System.out.println(info);
        for(int i=0; i<x.length; i++) {
            System.out.print(" " + x[i]);
        }
        System.out.println("");
    }
}
