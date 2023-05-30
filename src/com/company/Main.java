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
    public static void main(String[] args) throws UnsupportedAudioFileException, IOException
    {
        final File folder = new File("src/allSongs");
        final List<File> songList = Arrays.asList(folder.listFiles());

        List<List<Double>> allFrequencies = new ArrayList<>();
        List<String> nameSongs = new ArrayList<>();
        AudioInputStream audioInputStream;

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
        int dimension = 5;
        double[] centroid1 = new double[dimension];
        double[] centroid2 = new double[dimension];

        kMeansClusterManhattanD(centroid1, centroid2, dimension, allFrequencies, nameSongs);

        // preparation for Hierarchical clustering
        //hierarchicalClusteringManhattanD(allFrequencies, nameSongs);
    }

    static void buildCentroidK(double[] centroid, int dimension)
    {
        int[] tens = {100, 1000};
        for (int i = 0; i < dimension; i++)
            centroid[i] = Math.random() * getRandom(tens);
    }

    static boolean findFrequencies(AudioInputStream audioInputStream, List<List<Double>> allFrequencies) throws IOException
    {
        AudioFormat audioFormat = audioInputStream.getFormat();
        float sampleRate = audioFormat.getSampleRate();
        // WHAT NUMBER DO I USE FOR n -> justify!!!!!!!!! - i think it has to be a power of 2?
        int n = 4096;
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
        while (actualFrequency > 14080)
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

    // i have 5 "dimensions" now, not just X and Y
    // this is not workingggg, I have to change, probs smt with the distance calculations
    static void kMeansClusterManhattanD(double[] centroid1, double[] centroid2, int dimension,
                                        List<List<Double>> allFrequencies, List<String> nameSongs)
    {
        buildCentroidK(centroid1, dimension);
        buildCentroidK(centroid2, dimension);

        ArrayList<String> cluster1 = new ArrayList<>();
        ArrayList<String> cluster2 = new ArrayList<>();

        System.out.print("Initial position centroid 1: ");
        for (double c : centroid1)
            System.out.print(c + " ");
        System.out.println();
        System.out.print("Initial position centroid 2: ");
        for (double c : centroid2)
            System.out.print(c + " ");
        System.out.println();

        // What should be max value for iterations ?
        int counter = 0, maxIterations = 5;

        // for loop to check in which cluster each point goes
        // calculate distance from both centroids, go to the closest one
        // recalculate centroids and distances until the points are basically equally divided
        boolean proceed = true;
        double d1, d2;

        double sumXc1, sumYc1, sumZc1, sumAc1, sumBc1;
        double sumXc2, sumYc2, sumZc2, sumAc2, sumBc2;

        while (proceed)
        {
            for (int i = 0; i < allFrequencies.size(); i++)
            {
                d1 = calcManhD(centroid1, allFrequencies, i);
                d2 = calcManhD(centroid2, allFrequencies, i);

                if (d1 < d2)
                    cluster1.add(nameSongs.get(i));
                else
                    cluster2.add(nameSongs.get(i));
            }

            double half = allFrequencies.size() / 2.0;
            if (cluster1.size() == Math.floor(half) || cluster1.size() == Math.ceil(half)
                    || counter > maxIterations)
            {
                proceed = false;
            }

            else
            {
                counter++;
                sumXc1 = 0;
                sumYc1 = 0;
                sumZc1 = 0;
                sumAc1 = 0;
                sumBc1 = 0;
                sumXc2 = 0;
                sumYc2 = 0;
                sumZc2 = 0;
                sumAc2 = 0;
                sumBc2 = 0;

                //recalculate centroids
                for (String key1 : cluster1)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key1.equals(nameSongs.get(l)))
                        {
                            sumXc1 = sumXc1 + allFrequencies.get(l).get(0);
                            sumYc1 = sumYc1 + allFrequencies.get(l).get(1);
                            sumZc1 = sumZc1 + allFrequencies.get(l).get(2);
                            sumAc1 = sumAc1 + allFrequencies.get(l).get(3);
                            sumBc1 = sumBc1 + allFrequencies.get(l).get(4);
                        }
                    }
                }

                if (!cluster1.isEmpty())
                {
                    centroid1[0] = sumXc1 / cluster1.size();
                    centroid1[1] = sumYc1 / cluster1.size();
                    centroid1[2] = sumZc1 / cluster1.size();
                    centroid1[3] = sumAc1 / cluster1.size();
                    centroid1[4] = sumBc1 / cluster1.size();
                }


                for (String key2 : cluster2)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key2.equals(nameSongs.get(l)))
                        {
                            sumXc2 = sumXc2 + allFrequencies.get(l).get(0);
                            sumYc2 = sumYc2 + allFrequencies.get(l).get(1);
                            sumZc2 = sumZc2 + allFrequencies.get(l).get(2);
                            sumAc2 = sumAc2 + allFrequencies.get(l).get(3);
                            sumBc2 = sumBc2 + allFrequencies.get(l).get(4);
                        }
                    }
                }

                if (!cluster2.isEmpty())
                {
                    centroid2[0] = sumXc2 / cluster2.size();
                    centroid2[1] = sumYc2 / cluster2.size();
                    centroid2[2] = sumZc2 / cluster2.size();
                    centroid2[3] = sumAc2 / cluster2.size();
                    centroid2[4] = sumBc2 / cluster2.size();
                }

                cluster1.clear();
                cluster2.clear();
            }
        }

        for (String c : cluster1)
        {
            System.out.print(c + ", ");
        }
        System.out.println();
        for (String q : cluster2)
        {
            System.out.print(q + ", ");
        }
        System.out.println();
    }

    static void hierarchicalClusteringManhattanD(List<List<Double>> allFrequencies,
                                                 List<String> nameSongs)
    {
        // choose a better number for smallestDistance
        double sumX, sumY, sumZ, sumA, sumB, smallestDistance;
        int indexSong, index2 = 0, k;
        // What should be max value for iterations ?
        int counter = 0, maxIterations = 100, totalNumClusters = 2;

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

                // the problem might be that i'm moving songs, not clusters!!!!!!
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
                /* OG
                if (!search(clusters.get(i), nameSongs.get(i)))
                    clusters.get(i).add(nameSongs.get(i));

                if (!search(clusters.get(i), nameSongs.get(index2)))
                    clusters.get(i).add(nameSongs.get(index2));
                 */

                // delete repeated elements in other clusters
                /*
                k = 0;
                while (k < clusters.size())
                {
                    if (k != i)
                    {
                        if (search(clusters.get(k), nameSongs.get(i)))
                            // changed from this: clusters.get(k).remove(nameSongs.get(i));
                            clusters.get(k).clear();

                        if (search(clusters.get(k), nameSongs.get(index2)))
                            // changed from this: clusters.get(k).remove(nameSongs.get(index2));
                            clusters.get(k).clear();
                    }
                    k++;
                }

                 */

                /*
                //go through each cluster, if empty, delete it
                for (int j = 0; j < clusters.size(); j++)
                {
                    if (clusters.get(j).isEmpty())
                    {
                        clusters.remove(j);
                        centroids.remove(j);
                        distance.remove(j);
                    }
                }
                 */

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
    }

    static boolean checkDone(List<List<String>> clusters, int totalNumClusters)
    {
        if (clusters.size() <= totalNumClusters)
            return true;
        else return false;
    }

    static double calcManhD(double[] centroid1, List<List<Double>> allFrequencies, int i)
    {
        return Math.sqrt(Math.abs((centroid1[0] - allFrequencies.get(i).get(0)))
                + Math.abs((centroid1[1] - allFrequencies.get(i).get(1)))
                + Math.abs((centroid1[2] - allFrequencies.get(i).get(2)))
                + Math.abs((centroid1[3] - allFrequencies.get(i).get(3)))
                + Math.abs((centroid1[4] - allFrequencies.get(i).get(4))));
    }

    // clustering alorithms with euclidean distance
    /*
    static void kMeansCluster (double[] centroid1, double[] centroid2, List<List<Double>> allFrequencies,
                               List<String> nameSongs)
    {
        ArrayList<String> cluster1 = new ArrayList<>();
        ArrayList<String> cluster2 = new ArrayList<>();
        // What should be max value for iterations ?
        int counter = 0, maxIterations = 20;

        // for loop to check in which cluster each point goes
        // calculate distance from both centroids, go to the closest one
        // recalculate centroids and distances until the points are basically equally divided
        boolean proceed = true;
        double d1, d2;

        double sumXc1;
        double sumYc1;

        double sumXc2;
        double sumYc2;

        while (proceed)
        {
            for (int i = 0; i < allFrequencies.size(); i++)
            {
                d1 = Math.sqrt(Math.pow((centroid1[0] - allFrequencies.get(i).get(0)), 2) + Math.pow((centroid1[1] - allFrequencies.get(i).get(1)), 2));
                d2 = Math.sqrt(Math.pow((centroid2[0] - allFrequencies.get(i).get(0)), 2) + Math.pow((centroid2[1] - allFrequencies.get(i).get(1)), 2));


                if (d1 < d2)
                {
                    cluster1.add(nameSongs.get(i));
                }
                else
                {
                    cluster2.add(nameSongs.get(i));
                }
            }

            double half = allFrequencies.size() / 2.0;
            if (cluster1.size() == Math.floor(half) || cluster1.size() == Math.ceil(half)
                    || counter > maxIterations)
                proceed = false;


            else
            {
                counter++;
                sumXc1 = 0;
                sumYc1 = 0;
                sumXc2 = 0;
                sumYc2 = 0;
                //recalculate centroids
                for (String key1 : cluster1)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key1.equals(nameSongs.get(l)))
                        {
                            sumXc1 = sumXc1 + allFrequencies.get(l).get(0);
                            sumYc1 = sumYc1 + allFrequencies.get(l).get(1);
                        }
                    }
                }

                if (!cluster1.isEmpty())
                {
                    centroid1[0] = sumXc1 / cluster1.size();
                    centroid1[1] = sumYc1 / cluster1.size();
                }


                for (String key2 : cluster2)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key2.equals(nameSongs.get(l)))
                        {
                            sumXc2 = sumXc2 + allFrequencies.get(l).get(0);
                            sumYc2 = sumYc2 + allFrequencies.get(l).get(1);
                        }
                    }
                }

                if (!cluster2.isEmpty())
                {
                    centroid2[0] = sumXc2 / cluster2.size();
                    centroid2[1] = sumYc2 / cluster2.size();
                }

                cluster1.clear();
                cluster2.clear();
            }
        }

        for (String c : cluster1)
        {
            System.out.print(c + ", ");
        }
        System.out.println();
        for (String q : cluster2)
        {
            System.out.print(q + ", ");
        }
        System.out.println();
    }


    static void hierarchicalClustering(List<List<Double>> allFrequencies, List<String> nameSongs)
    {
        // choose a better number for smallestDistance
        double sumX, sumY, smallestDistance;
        int indexSong, index2 = 0, k;
        // What should be max value for iterations ?
        int counter = 0, maxIterations = 600;

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
                        Math.pow((centroids.get(i).get(0) - centroids.get(j).get(0)), 2) +
                                Math.pow((centroids.get(i).get(1) - centroids.get(j).get(1)), 2)));
            }
        }

        while (clusters.size() > 2 && counter < maxIterations)
        {
            counter++;
            for (int i = 0; i < clusters.size(); i++)
            {
                smallestDistance = 100000000;
                for (int j = 0; j < clusters.size(); j++)
                {
                    if (distance.get(i).get(j) < smallestDistance && distance.get(i).get(j) != 0)
                    {
                        smallestDistance = distance.get(i).get(j);
                        index2 = j;
                    }
                }

                if (!search(clusters.get(i), nameSongs.get(i)))
                    clusters.get(i).add(nameSongs.get(i));

                if (!search(clusters.get(i), nameSongs.get(index2)))
                    clusters.get(i).add(nameSongs.get(index2));

                // delete repeated elements in other clusters
                k = 0;
                while (k < clusters.size())
                {
                    if (k != i)
                    {
                        if (search(clusters.get(k), nameSongs.get(i)))
                            clusters.get(k).remove(nameSongs.get(i));

                        if (search(clusters.get(k), nameSongs.get(index2)))
                            clusters.get(k).remove(nameSongs.get(index2));
                    }
                    k++;
                }

                //go through each cluster, if empty, delete it
                for (int j = 0; j < clusters.size(); j++)
                {
                    if (clusters.get(j).isEmpty())
                        clusters.remove(j);
                }

                // clearing distance list to update it
                for (List<Double> doubles : distance) doubles.clear();

                // now i gotta calculate the centroid of each new cluster - OK
                for (int j = 0; j < clusters.size(); j++)
                {
                    sumX = 0;
                    sumY = 0;
                    for (String s : clusters.get(j))
                    {
                        indexSong = nameSongs.indexOf(s);
                        sumX = sumX + allFrequencies.get(indexSong).get(0);
                        sumY = sumY + allFrequencies.get(indexSong).get(1);
                    }
                    sumX = sumX / clusters.get(j).size();
                    sumY = sumY / clusters.get(j).size();

                    // what's going on here?!!!
                    centroids.get(j).add(sumX);
                    centroids.get(j).add(sumY);
                }

                // now I need to calculate the distance between the centroids and join more clusters
                // loop starts againnn
                for (int p = 0; p < centroids.size(); p++)
                {
                    for (int j = 0; j < centroids.size(); j++)
                    {
                        distance.get(p).add(j, Math.sqrt(
                                Math.pow((centroids.get(p).get(0) - centroids.get(j).get(0)), 2) +
                                        Math.pow((centroids.get(p).get(1) - centroids.get(j).get(1)), 2)));
                    }
                }
            }
        }

        for (List c : clusters)
            System.out.print(c);
        System.out.println();
    }
     */

    // K-MEANS WITH 3 CLUSTERS FIX THIS
    static void k3MeansClusterManhattanD(double[] centroid1, double[] centroid2, double[] centroid3,
                                         List<List<Double>> allFrequencies, List<String> nameSongs)
    {
        ArrayList<String> cluster1 = new ArrayList<>();
        ArrayList<String> cluster2 = new ArrayList<>();
        ArrayList<String> cluster3 = new ArrayList<>();

        // What should be max value for iterations ?
        int counter = 0, maxIterations = 6000;

        // for loop to check in which cluster each point goes
        // calculate distance from both centroids, go to the closest one
        // recalculate centroids and distances until the points are basically equally divided
        boolean proceed = true;
        double d1, d2, d3;

        double sumXc1, sumYc1, sumZc1, sumAc1, sumBc1;
        double sumXc2, sumYc2, sumZc2, sumAc2, sumBc2;
        double sumXc3, sumYc3, sumZc3, sumAc3, sumBc3;

        while (proceed)
        {
            for (int i = 0; i < allFrequencies.size(); i++)
            {
                d1 = calcManhD(centroid1, allFrequencies, i);
                d2 = calcManhD(centroid2, allFrequencies, i);
                d3 = calcManhD(centroid3, allFrequencies, i);

                if (d1 < d2 && d1 < d3)
                    cluster1.add(nameSongs.get(i));
                if (d2 < d1 && d2 < d3)
                    cluster2.add(nameSongs.get(i));
                if (d3 < d1 && d3 < d1)
                    cluster3.add(nameSongs.get(i));
            }

            double third = allFrequencies.size() / 3.0;
            if (cluster1.size() == Math.floor(third) && cluster2.size() == Math.ceil(third)
                    || cluster1.size() == Math.ceil(third) && cluster1.size() == Math.floor(third)
                    || counter > maxIterations)
            {
                proceed = false;
            }
            else
            {
                counter++;
                sumXc1 = 0;
                sumYc1 = 0;
                sumZc1 = 0;
                sumAc1 = 0;
                sumBc1 = 0;
                sumXc2 = 0;
                sumYc2 = 0;
                sumZc2 = 0;
                sumAc2 = 0;
                sumBc2 = 0;
                sumXc3 = 0;
                sumYc3 = 0;
                sumZc3 = 0;
                sumAc3 = 0;
                sumBc3 = 0;

                //recalculate centroids
                for (String key1 : cluster1)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key1.equals(nameSongs.get(l)))
                        {
                            sumXc1 = sumXc1 + allFrequencies.get(l).get(0);
                            sumYc1 = sumYc1 + allFrequencies.get(l).get(1);
                            sumZc1 = sumZc1 + allFrequencies.get(l).get(2);
                            sumAc1 = sumAc1 + allFrequencies.get(l).get(3);
                            sumBc1 = sumBc1 + allFrequencies.get(l).get(4);
                        }
                    }
                }

                if (!cluster1.isEmpty())
                {
                    centroid1[0] = sumXc1 / cluster1.size();
                    centroid1[1] = sumYc1 / cluster1.size();
                    centroid1[2] = sumZc1 / cluster1.size();
                    centroid1[3] = sumAc1 / cluster1.size();
                    centroid1[4] = sumBc1 / cluster1.size();
                }

                for (String key2 : cluster2)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key2.equals(nameSongs.get(l)))
                        {
                            sumXc2 = sumXc2 + allFrequencies.get(l).get(0);
                            sumYc2 = sumYc2 + allFrequencies.get(l).get(1);
                            sumZc2 = sumZc2 + allFrequencies.get(l).get(2);
                            sumAc2 = sumAc2 + allFrequencies.get(l).get(3);
                            sumBc2 = sumBc2 + allFrequencies.get(l).get(4);
                        }
                    }
                }

                if (!cluster2.isEmpty())
                {
                    centroid2[0] = sumXc2 / cluster2.size();
                    centroid2[1] = sumYc2 / cluster2.size();
                    centroid2[2] = sumZc2 / cluster2.size();
                    centroid2[3] = sumAc2 / cluster2.size();
                    centroid2[4] = sumBc2 / cluster2.size();
                }

                for (String key3 : cluster3)
                {
                    for (int l = 0; l < nameSongs.size(); l++)
                    {
                        if (key3.equals(nameSongs.get(l)))
                        {
                            sumXc3 = sumXc3 + allFrequencies.get(l).get(0);
                            sumYc3 = sumYc3 + allFrequencies.get(l).get(1);
                            sumZc3 = sumZc3 + allFrequencies.get(l).get(2);
                            sumAc3 = sumAc3 + allFrequencies.get(l).get(3);
                            sumBc3 = sumBc3 + allFrequencies.get(l).get(4);
                        }
                    }
                }

                if (!cluster3.isEmpty())
                {
                    centroid2[0] = sumXc3 / cluster3.size();
                    centroid2[1] = sumYc3 / cluster3.size();
                    centroid2[2] = sumZc3 / cluster3.size();
                    centroid2[3] = sumAc3 / cluster3.size();
                    centroid2[4] = sumBc3 / cluster3.size();
                }

                cluster1.clear();
                cluster2.clear();
                cluster3.clear();
            }
        }

        for (String c : cluster1)
        {
            System.out.print(c + ", ");
        }
        System.out.println();
        for (String q : cluster2)
        {
            System.out.print(q + ", ");
        }
        System.out.println();
        for (String q : cluster3)
        {
            System.out.print(q + ", ");
        }
        System.out.println();
    }
}