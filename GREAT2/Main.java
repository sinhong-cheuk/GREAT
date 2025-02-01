import javax.annotation.processing.SupportedSourceVersion;
import java.io.*;

import java.util.Random;


public class Main {

    public static void main(String[] args) throws IOException {
        String inputGraphFile = "/data1/graphFile";        // graph file path
        
        int reservoir_size = Integer.parseInt(args[1]);    // budget
        double alpha = Double.parseDouble(args[0]);        // alpha  
        
        
        Estimator estimator = new Estimator(reservoir_size, alpha);    // triangle counter
        System.out.println("alpha = " + alpha);
        System.out.println("budget = " + reservoir_size);
        double time0 = System.currentTimeMillis();
        run(estimator, inputGraphFile, "\t");
        double time1 = System.currentTimeMillis();
        double elpased_time = (time1 - time0) / 1000.0;
            
        estimator.output();                                // output local triangle file and calculate LAPE
        estimator.computeLAPE();
        
        System.out.println("elpased_time:" + elpased_time + "s");
            
        System.out.println("triangle detection:" + estimator.getDiscoverd_triangles());    
        System.out.println("global triangle estimation:" + String.format("%4f", estimator.getGlobalTriangle()));    
          

    }

    private static void run(Estimator estimator, String inputGraphFile, String delim) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inputGraphFile));

        long lineNum = 0;
        while(true) {
            String line = br.readLine();
            if(line == null) {
                break;
            }

            int[] edge = parseEdge(line, delim);

            estimator.processEdge(edge[0], edge[1]);       // GREAT2 processing streaming edge

            if ((++lineNum) % 100000000 == 0) {
                System.out.println("Number of edges processed: " + lineNum +", estimated number of global triangles: " + String.format("%4f", estimator.getGlobalTriangle()));
            }
        }
        System.out.println("GREAT2 terminated ...");
        System.out.println("Estimated number of global triangles: " + String.format("%4f", estimator.getGlobalTriangle()));

        br.close();
    }
    

    private static int[] parseEdge(String line, String delim) {
        String[] tokens = line.split(delim);
        int src = Integer.valueOf(tokens[0]);
        int dst = Integer.valueOf(tokens[1]);

        return new int[]{src, dst};
    }
}
