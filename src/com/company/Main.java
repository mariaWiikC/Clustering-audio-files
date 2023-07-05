package com.company;

import com.tambapps.fft4j.FastFouriers;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Main
{
    // apparently i'm getting empty clusters, aaaaaaaaa, find a way to fix this?
    public static void main(String[] args) throws UnsupportedAudioFileException, IOException
    {
        final File folder = new File("src/allSongs");
        final List<File> songList = Arrays.asList(folder.listFiles());

        List<List<Double>> allFrequencies = new ArrayList<>();
        List<String> nameSongs = new ArrayList<>();
        AudioInputStream audioInputStream;

        String[][] songKeys = {{"A#", "Gm"}, {"E"}, {"Ab"}, {"Cm", "F"}, {"Am"}, {"E"}, {"C", "F"}, {"A", "Fm"},
                {"D", "Em"}, {"Am", "C"}, {"C", "D"}, {"Eb"}, {"D"}, {"C", "Bb"}, {"G", "D"}, {"Ab", "Eb"},
                {"F#m", "Em"}, {"C"}, {"C"}, {"A"}, {"D", "C"}, {"F"}, {"Fm", "Eb"}, {"A", "Bb"}};

        // do I keep the subdominants? -> deleted them
        String[][] relatedKeys = {{"C", "Am", "F", "G", "Dm", "Em", "Cm"},
                {"G", "Em", "C", "D", "Am", "Bm", "Gm"},
                {"D", "Bm", "Cbm", "G", "A", "Em", "Fbm", "F#m", "Gbm", "Dm"},
                {"A", "F#m", "Gbm", "D", "E", "Fb", "Bm", "Cbm", "C#m", "Dbm", "Am"},
                {"E", "Fb", "C#m", "Dbm", "A", "B", "Cb", "F#m", "Gbm", "G#m", "Abm", "Em", "Fbm"},
                {"B", "Cb", "G#m", "Abm", "E", "Fb", "F#", "Gb", "C#m", "Dbm", "Dbm", "D#m", "Ebm", "Bm", "Cbm"},
                {"F#", "Gb", "D#m", "Ebm", "B", "Cb", "C#", "Db", "G#m", "Abm", "A#m", "Bbm", "F#m", "Gbm"},
                {"C#", "Db", "A#m", "Bbm", "F#", "Gb", "G#", "Ab", "D#m", "Ebm", "E#m", "Fm", "C#m", "Dbm"},
                {"G#", "Ab", "E#m", "Fm", "C#", "Db", "D#", "Eb", "A#m", "Bbm", "B#m", "Cm", "G#m", "Abm"},
                {"D#", "Eb", "B#m", "Cm", "G#", "Ab", "A#", "Bb", "E#m", "Fm", "Gm", "D#m", "Ebm"},
                {"A#", "Bb", "Gm", "D#", "Eb", "E#", "F", "B#m", "Cm", "Dm", "A#m", "Bbm"},
                {"F", "Dm", "Bb", "C", "Gm", "Am", "Fm"}};

        // for when I'm testing individual files
        /*
        audioInputStream = AudioSystem.getAudioInputStream(new File(("src/allSongs/middle-c.wav")));
        findFrequencies(audioInputStream, allFrequencies);
         */

        // FFT-ing all the songs
        for (File file : songList)
        {
            audioInputStream = AudioSystem.getAudioInputStream(file);
            if (findFrequencies(audioInputStream, allFrequencies))
                nameSongs.add(file.getName());
        }

        // preparation for K-means clustering
        // decide on number of clusters I will have
        int dimension = 5, k = 3;

         System.out.println(checkClusters(relatedKeys,
                kMeansClusterManhattanD(dimension, allFrequencies, nameSongs, k, songKeys)));

        // preparation for Hierarchical clustering
        // System.out.println(checkClusters(relatedKeys,
           //     hierarchicalClusteringManhattanD(allFrequencies, nameSongs, songKeys, k)));
    }

    static ArrayList<ArrayList<Integer>> checkClusters(String[][] relatedKeys,
                                                       ArrayList<ArrayList<ArrayList<String>>> keys)
    {
        ArrayList<ArrayList<Integer>> rNumRelKey = new ArrayList<>();
        for (int i = 0; i < keys.size(); i++)
        {
            rNumRelKey.add(new ArrayList<>());
        }

        for (int i = 0; i < keys.size(); i++)
        {
            // the only problem here is that it's getting two keys from the same song and taking it as two different
            // song keys, which might lead to a RelatedKeys row that is not the correct one
            ArrayList<String> clusterNoRep = new ArrayList<>();
            for (int j = 0; j < keys.get(i).size(); j++)
            {
                for (int l = 0; l < keys.get(i).get(j).size(); l++)
                {
                    String toAdd = keys.get(i).get(j).get(l);
                    if (!clusterNoRep.contains(toAdd))
                        clusterNoRep.add(toAdd);
                }
            }

            // check which row from related keys is the most similar
            ArrayList<Integer> similarityCounter = new ArrayList<>();
            for (int r = 0; r < relatedKeys.length; r++)
            {
                int numSimilar = 0;
                for (int c = 0; c < relatedKeys[r].length; c++)
                {
                    if (clusterNoRep.contains(relatedKeys[r][c]))
                        numSimilar++;
                }
                similarityCounter.add(numSimilar);
            }
            // what if there is more than one row with the same number of similarities? -> then ig it doesn't matter?
            // because I get the same result?
            int maxInd = findIndMaxValue(similarityCounter);


            // go through the keys of the cluster I'm at, if the row contains the key, 1 right, if not, 1 wrong.
            int correctCluster = 0, wrong = 0;
            for (int j = 0; j < keys.get(i).size(); j++)
            {
                boolean oneW = false, twoW = false;
                for (int l = 0; l < keys.get(i).get(j).size(); l++)
                {
                    // search in the row
                    // I'll just try to find a previous version where this was working, everything will be ok
                    for (int c = 0; c < relatedKeys[maxInd].length; c++)
                    {
                        if (keys.get(i).get(j).get(l).equals(relatedKeys[maxInd][c]))
                        {
                            correctCluster++;
                            if (oneW)
                                twoW = true;
                            else
                                oneW = true;
                        }
                    }
                }
                if (keys.get(i).get(j).size() == 2 && twoW)
                    wrong++;
            }
            // for these calculations, I'm getting the right number of data points, but I think the calculation
            // for the most similar row is wrong
            int total = correctCluster - wrong;
            // row of the right key group
            rNumRelKey.get(i).add(maxInd);
            // points in the right cluster
            rNumRelKey.get(i).add(total);
            //points in the wrong cluster -> smt wrong with this num
            rNumRelKey.get(i).add((keys.get(i).size() - total));
        }
        return rNumRelKey;
    }

    static int findIndMaxValue(ArrayList<Integer> arr)
    {
        int maxValue = arr.get(0), index = 0;
        for (int i = 0; i < arr.size(); i++)
        {
            if (arr.get(i) > maxValue)
            {
                index = i;
                maxValue = arr.get(i);
            }
        }
        return index;
    }

    static double[] buildCentroidK(int dimension)
    {
        double[] centroid = new double[dimension];
        int[] tens = {100};
        for (int i = 0; i < dimension; i++)
            centroid[i] = Math.random() * getRandom(tens);

        return centroid;
    }

    static boolean findFrequencies(AudioInputStream audioInputStream, List<List<Double>> allFrequencies) throws IOException
    {
        AudioFormat audioFormat = audioInputStream.getFormat();
        float sampleRate = audioFormat.getSampleRate();
        // WHAT NUMBER DO I USE FOR n -> justify!!!!!!!!! - i think it has to be a power of 2?
        // this is the num of samples i collect per second? too small, increase
        int n = 8192;
        byte[] array = new byte[n];

        int sample = audioInputStream.read(array);

        // input here
        double[] inputRe = new double[array.length];

        for (int i = 0; i < array.length; i++)
            inputRe[i] = array[i];

        // filling with 0's
        double[] inputIm = new double[array.length];
        for (int i = 0; i < array.length; i++)
            inputIm[i] = 0;

        double[] outputRe = new double[inputRe.length];
        double[] outputIm = new double[inputRe.length];

        FastFouriers.BASIC.transform(inputRe, inputIm, outputRe, outputIm);

        // calculating modulus = magnitudes
        double[] mod2 = new double[outputRe.length];
        for (int i = 0; i < outputRe.length; i++)
        {
            mod2[i] = Math.sqrt((outputRe[i] * outputRe[i]) + (outputIm[i] * outputIm[i]));
        }

        // this is me cutting the array in half
        double[] mod = new double[mod2.length / 2];
        for (int i = 0; i < mod2.length / 2; i++)
        {
            mod[i] = mod2[i];
        }

        double maxMod = 0;
        double iMaxMod = 0;
        for (int m = 0; m < mod.length; m++)
        {
            if (mod[m] > maxMod || mod[m] >= 14)
            {
                maxMod = mod[m];
                iMaxMod = m;
            }
        }

        // find all the peaks
        ArrayList<Integer> peaks = isPeak(mod);

        // AUTOMATIZAR ISSO
        // getting 3 highest peaks and their index -> all I need is the index
        double maxMod2 = 0, maxMod3 = 0, maxMod4 = 0, maxMod5 = 0, maxMod6 = 0;
        double iMaxMod2 = 0, iMaxMod3 = 0, iMaxMod4 = 0, iMaxMod5 = 0, iMaxMod6 = 0;
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod2 && mod[p] < maxMod)
            {
                maxMod2 = mod[p];
                iMaxMod2 = p;
            }
        }
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod3 && mod[p] < maxMod2)
            {
                maxMod3 = mod[p];
                iMaxMod3 = p;
            }
        }
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod4 && mod[p] < maxMod3)
            {
                maxMod4 = mod[p];
                iMaxMod4 = p;
            }
        }
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod5 && mod[p] < maxMod4)
            {
                maxMod5 = mod[p];
                iMaxMod5 = p;
            }
        }
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod6 && mod[p] < maxMod5)
            {
                maxMod6 = mod[p];
                iMaxMod6 = p;
            }
        }

        // printing the 3 frequencies = data that will be used for clustering
        //System.out.println(iMaxMod + " " + iMaxMod2 + " " + iMaxMod3);
        // adding the frequencies to the 2d array to automatically calculate distances and cluster
        if (iMaxMod != 0 && iMaxMod2 != 0 && iMaxMod3 != 0)
        {
            List<Double> fToAdd = new ArrayList<>();
            // System.out.println(iMaxMod2 + ", " + iMaxMod3 + ", " + iMaxMod4 + ", " + iMaxMod5);

            double actualFrequency = actualFrequency(iMaxMod, sampleRate, n);
            double actualFrequency2 = actualFrequency(iMaxMod2, sampleRate, n);
            double actualFrequency3 = actualFrequency(iMaxMod3, sampleRate, n);
            double actualFrequency4 = actualFrequency(iMaxMod4, sampleRate, n);
            double actualFrequency5 = actualFrequency(iMaxMod5, sampleRate, n);
            double actualFrequency6 = actualFrequency(iMaxMod6, sampleRate, n);

            // printing the frequencies just to check, with songs I know

            // System.out.println(actualFrequency);
            /*
            System.out.println(actualFrequency2);
            System.out.println(actualFrequency3);
            System.out.println(actualFrequency4);
            System.out.println(actualFrequency5);
            System.out.println(actualFrequency6);
             */

            // add values to list so I can clusterrr
            fToAdd.add(actualFrequency2);
            fToAdd.add(actualFrequency3);
            fToAdd.add(actualFrequency4);
            fToAdd.add(actualFrequency5);
            fToAdd.add(actualFrequency6);

            allFrequencies.add(fToAdd);
            return true;
        }
        else
            return false;
    }

    static double actualFrequency(double iMaxMod, float sampleRate, double n)
    {
        double actualFrequency = iMaxMod * sampleRate / n;
        // previous value was to keep the big numbers, now it's to try to approximate everything
        while (actualFrequency > 32.70)
        {
            actualFrequency = actualFrequency / 2;
        }
        return actualFrequency;
    }

    // this is not working -> supposed to substitute the repeated loops up there
    static double predFrequencies(ArrayList<Integer> peaks, double[] mod, double maxMod,
                                  double maxMod2, double iMaxMod2)
    {
        for (Integer p : peaks)
        {
            if (mod[p] > maxMod2 && mod[p] < maxMod)
            {
                maxMod2 = mod[p];
                iMaxMod2 = p;
            }
        }
        return iMaxMod2;
    }

    //method for finding peaks
    static ArrayList<Integer> isPeak(double[] arr)
    {
        ArrayList<Integer> iPeaks = new ArrayList<>();
        for (int i = 0; i < arr.length; i++)
        {
            if (i > 0 && arr[i] > arr[i - 1] && i < arr.length - 1 && arr[i] > arr[i + 1] && i > 15)
                iPeaks.add(i);
        }
        return iPeaks;
    }

    public static int getRandom(int[] array)
    {
        int rnd = new Random().nextInt(array.length);
        return array[rnd];
    }

    static boolean search(List<String> arr, String key)
    {
        for (String a : arr)
        {
            if (a.equals(key))
                return true;
        }
        return false;
    }

    static double[] findMinD(List<Double> Ds)
    {
        int indexMin = 0;
        double min = Ds.get(0);
        for (int i = 0; i < Ds.size(); i++)
        {
            if (Ds.get(i) < min)
            {
                min = Ds.get(i);
                indexMin = i;
            }
        }
        return new double[]{Ds.get(indexMin), indexMin};
    }

    // i have 5 "dimensions" now, not just X and Y
    static ArrayList<ArrayList<ArrayList<String>>> kMeansClusterManhattanD(int dimension,
                                                                           List<List<Double>> allFrequencies,
                                                                           List<String> nameSongs, int k,
                                                                           String[][] songKeys)
    {
        ArrayList<double[]> centroids = new ArrayList<>();
        for (int i = 0; i < k; i++)
            centroids.add(buildCentroidK(dimension));

        System.out.print("Initial position centroid 1: ");
        for (double c : centroids.get(0))
            System.out.print(c + " ");
        System.out.println();
        System.out.print("Initial position centroid 2: ");
        for (double c : centroids.get(1))
            System.out.print(c + " ");
        System.out.println();
        System.out.print("Initial position centroid 3: ");
        for (double c : centroids.get(2))
            System.out.print(c + " ");
        System.out.println();

        // What should be max value for iterations ?
        int counter = 0, maxIterations = 1000;

        // for loop to check in which cluster each point goes
        // calculate distance from both centroids, go to the closest one
        // recalculate centroids and distances until the points are basically equally divided
        boolean proceed = true;

        ArrayList<Double> Ds = new ArrayList<>();
        for (int i = 0; i < k; i++)
            Ds.add(0.0);

        ArrayList<ArrayList<String>> clusters = new ArrayList<>();
        for (int i = 0; i < k; i++)
            clusters.add(new ArrayList<>());

        ArrayList<ArrayList<Double>> sums = new ArrayList<>();
        for (int i = 0; i < k; i++)
        {
            sums.add(new ArrayList<>());
            for (int j = 0; j < dimension; j++)
                sums.get(i).add(0.0);
        }

        while (proceed)
        {
            for (int i = 0; i < allFrequencies.size(); i++)
            {
                for (int j = 0; j < k; j++)
                    Ds.set(j, calcManhD(centroids.get(j), allFrequencies, i));

                clusters.get((int) findMinD(Ds)[1]).add(nameSongs.get(i));
            }

            if (counter > maxIterations)
                proceed = false;

            else
            {
                counter++;
                for (ArrayList<Double> sum : sums)
                    for (int j = 0; j < dimension; j++)
                        sum.set(j, 0.0);

                for (int i = 0; i < clusters.size(); i++)
                {
                    //recalculate centroids
                    for (String key1 : clusters.get(i))
                    {
                        for (int l = 0; l < nameSongs.size(); l++)
                        {
                            if (key1.equals(nameSongs.get(l)))
                                for (int j = 0; j < dimension; j++)
                                    sums.get(i).set(j, sums.get(i).get(j) + allFrequencies.get(l).get(j));

                            if (!clusters.get(i).isEmpty())
                                for (int j = 0; j < dimension; j++)
                                    centroids.get(i)[j] = sums.get(i).get(j) / clusters.get(i).size();
                        }
                    }
                }

                for (ArrayList<String> c : clusters)
                    c.clear();
            }
        }

        // find position of song in the array and use the index to find corresponding key from key arr
        ArrayList<ArrayList<ArrayList<String>>> keys = new ArrayList<>();
        for (ArrayList<String> cluster : clusters)
        {
            ArrayList<ArrayList<String>> key = new ArrayList<>();
            keys.add(key);
            for (String c : cluster)
            {
                // c + " k: " +
                int indexC = nameSongs.indexOf(c);
                ArrayList<String> toA = new ArrayList<>();
                key.add(toA);
                for (int i = 0; i < songKeys[indexC].length; i++)
                {
                    System.out.print(songKeys[indexC][i] + ", ");
                    toA.add(songKeys[indexC][i]);
                }
                System.out.print("; ");
            }
            System.out.println();
        }

        return keys;
    }

    static ArrayList<ArrayList<ArrayList<String>>> hierarchicalClusteringManhattanD(List<List<Double>> allFrequencies,
                                                 List<String> nameSongs, String[][] songKeys, int totalNumClusters)
    {
        // choose a better number for smallestDistance
        double sumX, sumY, sumZ, sumA, sumB, smallestDistance;
        int indexSong, index2 = 0;
        // What should be max value for iterations ?
        int counter = 0, maxIterations = 100;

        // each list inside a list corresponds to a song, the values in the list are the distance
        // to other songs
        List<List<Double>> distance = new ArrayList<>();
        for (int i = 0; i < allFrequencies.size(); i++)
        {
            List<Double> d = new ArrayList<>();
            for (int j = 0; j < allFrequencies.size(); j++)
                d.add(0.0);
            distance.add(d);
        }

        List<List<Double>> centroids = new ArrayList<>();
        for (List<Double> song : allFrequencies)
        {
            List<Double> c = new ArrayList<>();
            c.add(song.get(0));
            c.add(song.get(1));
            c.add(song.get(2));
            c.add(song.get(3));
            c.add(song.get(4));
            centroids.add(c);
        }

        //creating clusters, one for each data point
        List<List<String>> clusters = new ArrayList<>();
        for (int i = 0; i < allFrequencies.size(); i++)
        {
            List<String> cluster = new ArrayList<>();
            cluster.add(nameSongs.get(i));
            clusters.add(cluster);
        }

        //calculating proximity matrix - OK
        for (int i = 0; i < centroids.size(); i++)
        {
            for (int j = 0; j < centroids.size(); j++)
            {
                distance.get(i).set(j, Math.sqrt(
                        Math.abs((centroids.get(i).get(0) - centroids.get(j).get(0))) +
                                Math.abs((centroids.get(i).get(1) - centroids.get(j).get(1))) +
                                Math.abs((centroids.get(i).get(2) - centroids.get(j).get(2))) +
                                Math.abs((centroids.get(i).get(3) - centroids.get(j).get(3))) +
                                Math.abs((centroids.get(i).get(4) - centroids.get(j).get(4)))));
            }
        }

        // is my while in the right place? or must it be after this for loo
        // i dunno, smt is not right
        while (clusters.size() > totalNumClusters && counter < maxIterations)
        {
            for (int i = 0; i < clusters.size(); i++)
            {
                smallestDistance = 1000000000;
                for (int j = 0; j < clusters.size(); j++)
                {
                    // this gets thw two clusters with smallest distance
                    if (distance.get(i).get(j) < smallestDistance && distance.get(i).get(j) != 0)
                    {
                        smallestDistance = distance.get(i).get(j);
                        index2 = j;
                    }
                }

                // move an entire cluster
                // find cluster with that has the smallest distance, pass it to the new cluster
                // and delete it
                for (String s : clusters.get(index2))
                    clusters.get(i).add(s);
                clusters.get(index2).clear();
                clusters.remove(index2);
                centroids.remove(index2);
                distance.remove(index2);

                if (checkDone(clusters, totalNumClusters)) break;

                // clearing distance list to update it
                for (List<Double> doubles : distance) doubles.clear();

                // now i gotta calculate the centroid of each new cluster - OK
                for (int j = 0; j < clusters.size(); j++)
                {
                    sumX = 0;
                    sumY = 0;
                    sumZ = 0;
                    sumA = 0;
                    sumB = 0;
                    for (String s : clusters.get(j))
                    {
                        indexSong = nameSongs.indexOf(s);
                        sumX = sumX + allFrequencies.get(indexSong).get(0);
                        sumY = sumY + allFrequencies.get(indexSong).get(1);
                        sumZ = sumZ + allFrequencies.get(indexSong).get(2);
                        sumA = sumA + allFrequencies.get(indexSong).get(3);
                        sumB = sumB + allFrequencies.get(indexSong).get(4);
                    }
                    sumX = sumX / clusters.get(j).size();
                    sumY = sumY / clusters.get(j).size();
                    sumZ = sumZ / clusters.get(j).size();
                    sumA = sumA / clusters.get(j).size();
                    sumB = sumB / clusters.get(j).size();

                    // what's going on here?!!!
                    centroids.get(j).set(0, sumX);
                    centroids.get(j).set(1, sumY);
                    centroids.get(j).set(2, sumZ);
                    centroids.get(j).set(3, sumA);
                    centroids.get(j).set(4, sumB);
                }

                // now I need to calculate the distance between the centroids and join more clusters
                // loop starts againnn
                for (int p = 0; p < centroids.size(); p++)
                {
                    for (int j = 0; j < centroids.size(); j++)
                    {
                        distance.get(p).add(j, Math.sqrt(
                                Math.abs((centroids.get(p).get(0) - centroids.get(j).get(0))) +
                                        Math.abs((centroids.get(p).get(1) - centroids.get(j).get(1))) +
                                        Math.abs((centroids.get(p).get(2) - centroids.get(j).get(2))) +
                                        Math.abs((centroids.get(p).get(3) - centroids.get(j).get(3))) +
                                        Math.abs((centroids.get(p).get(4) - centroids.get(j).get(4)))));
                    }
                }
            }
            counter++;
        }

        for (List c : clusters)
            System.out.println(c);
        System.out.println(counter);

        // creating the list used as a parameter to check "goodness" of clusters
        ArrayList<ArrayList<ArrayList<String>>> keys = new ArrayList<>();
        for (List<String> cluster : clusters)
        {
            ArrayList<ArrayList<String>> key = new ArrayList<>();
            keys.add(key);
            for (String c : cluster)
            {
                // c + " k: " +
                int indexC = nameSongs.indexOf(c);
                ArrayList<String> toA = new ArrayList<>();
                key.add(toA);
                for (int i = 0; i < songKeys[indexC].length; i++)
                {
                    System.out.print(songKeys[indexC][i] + ", ");
                    toA.add(songKeys[indexC][i]);
                }
                System.out.print("; ");
            }
            System.out.println();
        }

        return keys;
    }

    static boolean checkDone(List<List<String>> clusters, int totalNumClusters)
    {
        return clusters.size() <= totalNumClusters;
    }

    static double calcManhD(double[] centroid1, List<List<Double>> allFrequencies, int i)
    {
        return Math.sqrt(Math.abs((centroid1[0] - allFrequencies.get(i).get(0)))
                + Math.abs((centroid1[1] - allFrequencies.get(i).get(1)))
                + Math.abs((centroid1[2] - allFrequencies.get(i).get(2)))
                + Math.abs((centroid1[3] - allFrequencies.get(i).get(3)))
                + Math.abs((centroid1[4] - allFrequencies.get(i).get(4))));
    }
}