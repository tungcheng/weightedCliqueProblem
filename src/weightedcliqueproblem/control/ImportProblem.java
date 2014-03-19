/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package weightedcliqueproblem.control;

import java.io.File;
import java.util.Scanner;
import weightedcliqueproblem.model.Problem;

/**
 *
 * @author TungNT
 */
public class ImportProblem {
    
    public ImportProblem() {
        
    }
    
    public Problem importFromFile(File file) throws Exception {
        Scanner in = new Scanner(file);
        int b = in.nextInt();
        int n = in.nextInt();
        Problem pb = new Problem(n);
        pb.setB(b);
        double t;
        for(int i=0; i<n; i++) {
            for(int j=0; j<n; j++) {
                t = in.nextDouble();
                pb.setC((-1)*t, i, j);
            }
        }
        if(!pb.isSymmetric()) {
            throw new Exception();
        }
        pb.makePositiveSemidefinite();
        return pb;
    }
}
