import it.unimi.dsi.fastutil.doubles.DoubleOpenHashSet;
import it.unimi.dsi.fastutil.ints.*;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.Int2BooleanOpenHashMap;

import java.io.*;
import java.util.*;


public class Estimator {

    private Int2ObjectOpenHashMap<Int2IntOpenHashMap> neighbors = new Int2ObjectOpenHashMap<>();  // graph composed of the sampled edges
    
    private Int2DoubleOpenHashMap nodeToCount = new Int2DoubleOpenHashMap(); // local triangle counts

    private int maxID = -1;            

    private double globalTriangle = 0;                      // global triangles
    private double pre_global = 0;                          // for calculate triangle estimatation increment per round


    private int k;                                          // size of the reservoir
    private int[][] reservoir;
    private double t = 0;                                   // number of streaming edges processed so far
    private double interval = 0;                            // total interval of new triangles
    private double aver_interval = 0;                       // average interval of new triangles in a round
    private double z;                                       // z = (1 - x)^y
    private int round_bound;                                // when computation round less than round bound, we use the init_alpha, not adaptive
    private double init_alpha;                              // initial alpha for the top-(round_bound) round
    private double p = 1;                                   // sampling probability     
  
    private int next_slot_index = 1;                        // for top-k edges in reservoir sampling
    private int store_index = 0;
    private int empty_slot = 0;                             // current empty slot of the reservoir

    private int[] delete_index;                             // current empty slot of the reservoir
    private int[] remain_index;                             // position of the remaining edges

    private List<Integer> numbers = new ArrayList<>();      // for generating the random index 

    private double cur_round = 1;                           // current computation round of sampling scheme
    private double alpha;                                   // discard probability of a reservoir edge
    private double survive_rate;                            // survive probability of an edge, 1-alpha
        
    // cache of the (1-alpha)^r, accelerate the calculation of a triangle's probability    
    private final int max_round = 10000;
    private double[] survive_rate_array = new double[max_round];     
    private double[][] survive_rate_cache = new double[max_round][max_round];
             
    private double[][] p_and_round;                         // for p_uv and r_uv
    private double[] sample_time;

    private final Random random = new Random();

    private int N;                                          // number of deleted edges

    private long discoverd_triangles = 0;                   
    private int discoverd_triangles_per_round = 0;          // for caculating average interval

    public Estimator(int sizeOfReservoir, double z, int round_bound, double init_alpha) {
        this.reservoir = new int[2][sizeOfReservoir + 1];   // we have special use of index '0'     
        this.k = sizeOfReservoir;
        this.p_and_round = new double[2][sizeOfReservoir + 1];
        
        this.z = z;
        this.round_bound = round_bound;
        this.init_alpha = init_alpha;
        

        survive_rate_array[0] = 1;

        survive_rate_cache[0][0] = 1;
        survive_rate_cache[1][1] = 1;


        for (int i = 1; i < k + 1; i++) {
            numbers.add(i);
        }

        delete_index = new int[N];
        remain_index = new int[k - N];
    }

    public void processEdge(int src, int dst) {
        if(src == dst) { //ignore self loop
            return;
        }

        t++;
        
        // get maxID for output local triangle file
        if (src > maxID) {
            maxID = src;
        }
        
        if (dst > maxID) {
            maxID = dst;
        }
        

        if (t == k) {
            
            alpha = init_alpha;
            this.N = (int) (k * alpha);

            survive_rate = 1 - alpha;
            survive_rate_array[(int)cur_round] = survive_rate;
        }

        
        count(src, dst); //count triangles involved

        if (t < k + 1) {
            // top-k edges in reservoir sampling
            reservoir[0][next_slot_index] = src;
            reservoir[1][next_slot_index] = dst;

            p_and_round[0][next_slot_index] = 1.0;
            p_and_round[1][next_slot_index] = 1.0;


            sample(src, dst, next_slot_index);
            next_slot_index++;

        } else {
            // now the sampling probability turns to be p * (1 - alpha)
            if (empty_slot == 0) {
                // reservoir is full now, we need to randomly discard N edges
                
                // how many triangles we estimated in last round
                double increment = globalTriangle - pre_global;
                pre_global = globalTriangle;

                if (cur_round > round_bound) {
                    System.out.println("round "+ cur_round + " global triangle estimation: " + String.format("%4f", + globalTriangle));
                    System.out.println("round "+ cur_round + " triangle interval: " + String.format("%4f", + interval));
                    System.out.println("round "+ cur_round + " increment: " + String.format("%4f", + increment));
                    System.out.println("round "+ cur_round + " triangle detection: " + String.format("%4d", + discoverd_triangles_per_round));

                    if (discoverd_triangles_per_round != 0) {
                        aver_interval =  interval  / discoverd_triangles_per_round;
                        System.out.println("round "+ cur_round + " average triangle interval: " + aver_interval);
                        System.out.println("round "+ cur_round + " z value: " + String.format("%4f", + z));

                        alpha = generateAlphaByInterval(aver_interval);
                        alpha = Math.max(init_alpha, alpha);
                        
                    } else {
                        System.out.println("round "+ cur_round + " average triangle interval: Infinity");
                        //alpha = Math.max(0.1, alpha);
                        alpha = init_alpha;    
                    }
                    
                    System.out.println("round "+ cur_round + " z value: " + String.format("%4f", + z));

                    System.out.println("round "+ cur_round + " alpha value: " + String.format("%4f", + alpha));
                    System.out.println();

                    cur_round++;
                    interval = 0;
                    discoverd_triangles_per_round = 0;

                    this.N = (int) (k * alpha);

                    survive_rate = 1 - alpha;

                    // update the survive_rate_cache, then we can calculate a triangle's probability faster
                    survive_rate_array[(int)cur_round] = survive_rate;
                    survive_rate_cache[(int)cur_round][(int)cur_round] = survive_rate;

                    for (int i = 0; i < cur_round; i++) {
                        survive_rate_cache[i][(int)cur_round] = survive_rate * survive_rate_cache[i][(int)cur_round - 1];
                    }
                    survive_rate_cache[(int)cur_round + 1][(int)cur_round] = 1;
                    
                    // p = p * (1 - alpha);      // another sampling probablity
                    
                    
                    System.out.println();
                } else {
                    System.out.println("round "+ cur_round + " global triangle estimation: " + String.format("%4f", + globalTriangle));
                    System.out.println("round "+ cur_round + " triangle interval: " + String.format("%4f", + interval));
                    System.out.println("round "+ cur_round + " increment: " + String.format("%4f", + increment));
                    System.out.println("round "+ cur_round + " triangle detection: " + String.format("%4d", + discoverd_triangles_per_round));
                    aver_interval =  interval  / discoverd_triangles_per_round;
                    System.out.println("round "+ cur_round + " aver_interval: " + aver_interval);

                    System.out.println("round "+ cur_round + " average triangle interval: " + aver_interval);
                    System.out.println("round "+ cur_round + " z value: " + String.format("%4f", + z));
                    System.out.println("round "+ cur_round + " alpha value: " + String.format("%4f", + alpha));
                    System.out.println();


                    cur_round++;
                    interval = 0;
                    discoverd_triangles_per_round = 0;
                    alpha = init_alpha;

                    this.N = (int) (k * alpha);

                    survive_rate = 1 - alpha;

                    // update the survive_rate_cache, then we can calculate a triangle's probability faster
                    survive_rate_array[(int)cur_round] = survive_rate;
                    survive_rate_cache[(int)cur_round][(int)cur_round] = survive_rate;

                    for (int i = 0; i < cur_round; i++) {
                        survive_rate_cache[i][(int)cur_round] = survive_rate * survive_rate_cache[i][(int)cur_round - 1];
                    }
                    survive_rate_cache[(int)cur_round + 1][(int)cur_round] = 1;
                    
                    // p = p * (1 - alpha);      // another sampling probablity
                    System.out.println();
                }


                
                randomIndex();

                if (alpha <= 0.5) {
                    // remove edges
                    for (int i = 0; i < N; i++) {
                        int index_tobeRemove = delete_index[i];

                        int src_tobeRemove = reservoir[0][index_tobeRemove];
                        int dst_tobeRemove = reservoir[1][index_tobeRemove];

                        deleteEdge(src_tobeRemove, dst_tobeRemove);
                    }

                } else {
                    //when alpha > 0.5, save the remaining would be faster
                    Int2ObjectOpenHashMap<Int2IntOpenHashMap> temp_neighbors = new Int2ObjectOpenHashMap<>();

                    int remaining = k - N;

                    for (int i = 0; i < remaining; i++) {
                        int index_tobeRemain = remain_index[i];

                        int src_tobeRemain = reservoir[0][index_tobeRemain];
                        int dst_tobeRemain = reservoir[1][index_tobeRemain];

                        temp_neighbors.computeIfAbsent(src_tobeRemain, k -> new Int2IntOpenHashMap())
                                .put(dst_tobeRemain, index_tobeRemain);

                        temp_neighbors.computeIfAbsent(dst_tobeRemain, k -> new Int2IntOpenHashMap())
                                .put(src_tobeRemain, index_tobeRemain);

                    }

                    neighbors = temp_neighbors;
                }

                empty_slot = N;

            }

            // if there are still have some slots for coming edges, keep sampling
            if (empty_slot > 0) {
      
                p = k / t;
                double randomValue = random.nextDouble();
                if (randomValue < p) {
                    // sample the coming edge
                   
                    int insertIndex = delete_index[empty_slot - 1];
                    reservoir[0][insertIndex] = src;
                    reservoir[1][insertIndex] = dst;

                    p_and_round[0][insertIndex] = p;
                    p_and_round[1][insertIndex] = cur_round;


                    sample(src, dst, insertIndex);

                    empty_slot--;

                }

            } else return;
        }


    }
    /**
     * sample an edge to the subgraph, and store the reservoir index of this edge
     * then we can get sampling probability and sampling round of an edge
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
    private void sample(int src, int dst, int storeIndex) {
        neighbors.computeIfAbsent(src, k -> new Int2IntOpenHashMap()).put(dst, storeIndex);
        neighbors.computeIfAbsent(dst, k -> new Int2IntOpenHashMap()).put(src, storeIndex);
    }

    /**
     * delete an edge from the subgraph
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
    private void deleteEdge(int src, int dst) {
        Int2IntOpenHashMap map = neighbors.get(src);
        if (map != null) {
            map.remove(dst);
            if (map.isEmpty()) {
                neighbors.remove(src);
            }
        }

        map = neighbors.get(dst);
        if (map != null) {
            map.remove(src);
            if (map.isEmpty()) {
                neighbors.remove(dst);
            }
        }
    }

    /**
     * count triangles incident to the given edge
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
    private void count(int src, int dst) {

        // if this edge has a new node, there cannot be any triangles
        if(!neighbors.containsKey(src) || !neighbors.containsKey(dst)) {
            return;
        }

        // source node to neighbors
        Int2IntOpenHashMap srcMap = neighbors.get(src);

        // destination node to neighbors
        Int2IntOpenHashMap dstMap = neighbors.get(dst);

        if(srcMap.size() > dstMap.size()) {
            Int2IntOpenHashMap temp = srcMap;
            srcMap = dstMap;
            dstMap = temp;
        }

        // the sum of counts increased
        double countSum = 0;


        // find common neighborhood and update counter
        for (Int2IntMap.Entry entry : srcMap.int2IntEntrySet()) {
            int neighbor = entry.getIntKey();
            int indexDst = dstMap.get(neighbor);

            
            if (indexDst != 0) {    // we have special uses of index '0', it means we find a common neighborhood
                discoverd_triangles++;
                discoverd_triangles_per_round++;

                if (t < k + 1) {

                    countSum += 1;
                    nodeToCount.addTo(neighbor, 1); // update the local triangle count of the common neighbor
                } else {

                    int indexSrc = entry.getIntValue();

                    
                    // calculate the interval of triangle (src, dst, neighbor)
                    int srcSampleRound = (int)p_and_round[1][indexSrc];
                    int dstSampleRound = (int)p_and_round[1][indexDst];
                    interval += 2 * cur_round - srcSampleRound - dstSampleRound;

                    // calculate a triangle's probability
                    double src_p, dst_p;
                    if (srcSampleRound == cur_round) {
                        src_p = p_and_round[0][indexSrc];
                    } else {
                        src_p = p_and_round[0][indexSrc] * survive_rate_cache[srcSampleRound + 1][(int)cur_round];
                    }

                    if (dstSampleRound == cur_round) {
                        dst_p = p_and_round[0][indexDst];
                    } else {
                        dst_p = p_and_round[0][indexDst] * survive_rate_cache[dstSampleRound + 1][(int)cur_round];
                    }

                    double count = 1 / (src_p * dst_p);

                    
                    countSum += count;
                    nodeToCount.addTo(neighbor, count); // update the local triangle count of the common neighbor
                }

            }
        }


        if(countSum > 0) {
            nodeToCount.addTo(src, countSum); // update the local triangle count of the source node
            nodeToCount.addTo(dst, countSum); // update the local triangle count of the destination node
            globalTriangle += countSum; // update the global triangle count
            
        }  
    }

    /**
     * generate the random index of deleted edges of remaining edges
     */
    private void randomIndex() {

        if (alpha <= 0.5) {
            
            // generate N different random number
            IntOpenHashSet indicesToDelete = new IntOpenHashSet();
            while (indicesToDelete.size() < N) {
                indicesToDelete.add(1 + random.nextInt(k));
            }
            delete_index = indicesToDelete.toIntArray();

        } else {
            // shuffle a 0-k sequence
            List<Integer> number = numbers;
            Collections.shuffle(number, random);

            delete_index = number.subList(0, N).stream().mapToInt(i -> i).toArray();
            remain_index = number.subList(N, number.size()).stream().mapToInt(i -> i).toArray();

        }

    }

    public double getGlobalTriangle() {
        return globalTriangle;
    }

    public Int2DoubleMap getLocalTriangle() {
        return nodeToCount;
    }

    public double getAlpha() {
        return alpha;
    }

    public long getDiscoverd_triangles() {
        return discoverd_triangles;
    }

    /**
     * generate the alpha by triangle interval
     * @param aver_interval Average triangle interval of current round
     * @return Alpha for next round
     */
    private double generateAlphaByInterval(double aver_interval) {
        double z_temp = z;
        double x = Math.pow(z_temp, 1.0 / aver_interval);

        return Math.round((1 - x) * 10000) / 10000.0;
    }
    
    
    /**
     * output local triangle estimation to file
     */ 
    public void output() throws IOException {
        String fileName = "/data1/local-GREAT+.txt";        // local triangle estimation file path
        int v_count = 0;                                    // vertices of discovered triangles

        BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
        for (int i = 0; i <= maxID; i++) {
            if (i % 50000000 == 0) {
                System.out.println("writing node: " + i);
            }
            double count = nodeToCount.getOrDefault(i, 0.0);
            if (count > 0) {
                v_count++;
            }
            writer.write(i + "\t" + count + "\n");
        }
        
        System.out.println("vertices of discovered triangles:" + v_count);
        writer.close();
    }

    /**
     * caculate local triangle estimation error
     */ 
    public double computeLAPE() {
        String algorithmOutputFile = "/data1/local-GREAT+.txt";      // local triangle estimation file path
        String groundTruthFile = "/data1/groundTruthFile.txt";     // local triangle groundtruth file path
  
        double averageLAPE = 0;
        try (
                BufferedReader groundTruthReader = new BufferedReader(new FileReader(groundTruthFile));
                BufferedReader algorithmOutputReader = new BufferedReader(new FileReader(algorithmOutputFile))
        ) {
            String groundTruthLine; 
            String algorithmOutputLine;  

            double totalLAPE = 0;
            int vertexCount = 0;

            // calculate LAPE
            while ((groundTruthLine = groundTruthReader.readLine()) != null &&
                    (algorithmOutputLine = algorithmOutputReader.readLine()) != null) {
                String[] groundTruthParts = groundTruthLine.split("\\s+");
                String[] algorithmOutputParts = algorithmOutputLine.split("\\s+");

                
                double groundTruthValue = Double.parseDouble(groundTruthParts[1]); 
                double algorithmOutputValue = Double.parseDouble(algorithmOutputParts[1]);

                double lape = Math.abs(algorithmOutputValue - groundTruthValue) / (groundTruthValue + 1);
                totalLAPE += lape;

                vertexCount++;
            }
            System.out.println("vertex num: " + vertexCount);

            averageLAPE = totalLAPE / vertexCount;
            System.out.println("Average Local Absolute Percentage Error (LAPE): " + String.format("%4f", averageLAPE));


        } catch (IOException e) {
            e.printStackTrace();
        }
        return averageLAPE;
    }
}
