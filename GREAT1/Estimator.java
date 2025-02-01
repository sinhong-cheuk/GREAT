import it.unimi.dsi.fastutil.doubles.DoubleOpenHashSet;
import it.unimi.dsi.fastutil.ints.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import java.io.*;



public class Estimator {

    private Int2ObjectOpenHashMap<IntOpenHashSet> neighbors = new Int2ObjectOpenHashMap<>();  // graph composed of the sampled edges

    private Int2DoubleOpenHashMap nodeToCount = new Int2DoubleOpenHashMap(); // local triangle counts
    private double globalTriangle = 0;                           // global triangles
    private int maxID = -1;                                      // vertex ID

    private int[][] reservoir;
    private int k;                                               // size of the reservoir
    private double t = 0;                                        // number of streaming edges processed so far
    private int next_slot_index = 0;                             // for top-k edges in reservoir sampling
    private int empty_slot = 0;                                  // current empty slot of the reservoir
    private int[] delete_index;                                  // position of the deleted edges in reservoir
    private int[] remain_index;                                  // position of the remaining edges in reservoir
    
    private long discoverd_triangles = 0;

    private double cur_round = 1;                                // current computation round of sampling scheme
    private double cur_sample_p = 1;                             // current sampling probability
    private double weight = 1;                                   // probability of a triangle
    private double alpha;                                        // edge discard rate
    private double survive_rate;                                 // survive probability of an edge, 1-alpha
    private double[] survive_rate_array = new double[10000];     // cache of the (1-alpha)^r, accelerate the calculation of the probability
    private List<Integer> numbers = new ArrayList<>();           // for generating the random index
    
    private final Random random = new Random();

    private int N;                                               // number of deleted edges



    public Estimator(int sizeOfReservoir, double alpha) {
        this.reservoir = new int[2][sizeOfReservoir];
        this.k = sizeOfReservoir;
        this.alpha = alpha;
        this.N = (int) (k * alpha);

        for (int i = 0; i < k; i++) {
            numbers.add(i);
        }

        delete_index = new int[N];
        remain_index = new int[k - N];
    }


    /**
     * process the given edge
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
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


        count(src, dst); //count triangles involved


        if (t < k + 1) {
            // top-k edges in reservoir sampling
            reservoir[0][next_slot_index] = src;
            reservoir[1][next_slot_index] = dst;

            sample(src, dst);
            next_slot_index++;

        } else {
            // now the sampling probability turns to be (1 - alpha) ^ r
            if (empty_slot == 0) {
                // reservoir doesn't have empty slot now, we need to randomly discard N edges
                // update triangles probability
                cur_round++;
                cur_sample_p = cur_sample_p * (1 - alpha);
                weight = cur_sample_p * cur_sample_p;


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

                    Int2ObjectOpenHashMap<IntOpenHashSet> temp_neighbors = new Int2ObjectOpenHashMap<>();

                    int remaining = k - N;

                    for (int i = 0; i < remaining; i++) {
                        int index_tobeRemain = remain_index[i];

                        int src_tobeRemain = reservoir[0][index_tobeRemain];
                        int dst_tobeRemain = reservoir[1][index_tobeRemain];

                        temp_neighbors.computeIfAbsent(src_tobeRemain, k -> new IntOpenHashSet())
                                .add(dst_tobeRemain);

                        temp_neighbors.computeIfAbsent(dst_tobeRemain, k -> new IntOpenHashSet())
                                .add(src_tobeRemain);

                    }

                    neighbors = temp_neighbors;
                }


                empty_slot = N;

            }

            // if there are still have some slots for coming edges, keep sampling
            if (empty_slot > 0) {

                double randomValue = random.nextDouble();
                if (randomValue < cur_sample_p) {

                    // sample the coming edge in the empty slot
                    int insertIndex = delete_index[empty_slot - 1];
                    reservoir[0][insertIndex] = src;
                    reservoir[1][insertIndex] = dst;

                    sample(src, dst);
                    empty_slot--;

                }

            } else return;
        }


    }
    
    /**
     * sample an edge to the subgraph
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
    private void sample(int src, int dst) {

        if(!neighbors.containsKey(src)) {
            neighbors.put(src, new IntOpenHashSet());
        }
        neighbors.get(src).add(dst);

        if(!neighbors.containsKey(dst)) {
            neighbors.put(dst, new IntOpenHashSet());
        }
        neighbors.get(dst).add(src);
    }
    
    /**
     * delete an edge from the subgraph
     * @param src source node of the given edge
     * @param dst destination node of the given edge
     */
    private void deleteEdge(int src, int dst) {
        IntOpenHashSet map = neighbors.get(src);
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
        IntOpenHashSet srcMap = neighbors.get(src);

        // destination node to neighbors
        IntOpenHashSet dstMap = neighbors.get(dst);

        if(srcMap.size() > dstMap.size()) {
            IntOpenHashSet temp = srcMap;
            srcMap = dstMap;
            dstMap = temp;
        }

        // the sum of counts increased
        double countSum = 0;
        // find common neighborhood and update counter
        for(int neighbor : srcMap) {
            if (dstMap.contains(neighbor)) {
                discoverd_triangles++;              // algorithm detect a triangle

                double count = 1 / weight;

                countSum += count;
                nodeToCount.addTo(neighbor, count); // update the local triangle count of the common neighbor
            }
        }

        if(countSum > 0) {
            nodeToCount.addTo(src, countSum); // update the local triangle count of the source node
            nodeToCount.addTo(dst, countSum); // update the local triangle count of the destination node
            globalTriangle += countSum;       // update the global triangle count
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
                indicesToDelete.add(random.nextInt(k));
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
     * output local triangle estimation to file
     */
    public void output() throws IOException {
        String fileName = "/data1/local-GREAT1.txt";                  // local triangle estimation file path
        int v_count = 0;                                              // vertices of discovered triangles

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
    public void computeLAPE() {
        String algorithmOutputFile = "/data1/local-GREAT1.txt";      // local triangle estimation file path
        String groundTruthFile = "/data1/local-groundTruthFile.txt";           // local triangle groundtruth file path

        try (
                BufferedReader groundTruthReader = new BufferedReader(new FileReader(groundTruthFile));
                BufferedReader algorithmOutputReader = new BufferedReader(new FileReader(algorithmOutputFile))
        ) {
            String groundTruthLine;
            String algorithmOutputLine;

            double totalLAPE = 0;
            int vertexCount = 0;

            // calculating LAPE
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

            double averageLAPE = totalLAPE / vertexCount;
            System.out.println("Average Local Absolute Percentage Error (LAPE): " + String.format("%4f", averageLAPE));

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
