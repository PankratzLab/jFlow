package org.genvisis.one.JL;

// package one.JL;
//
// import htsjdk.samtools.SAMSequenceDictionary;
// import htsjdk.samtools.SAMSequenceRecord;
//
// import java.io.UnsupportedEncodingException;
// import java.util.ArrayList;
// import java.util.HashMap;
// import java.util.List;
// import java.util.Map;
// import java.util.concurrent.Callable;
//
// import common.Logger;
// import common.WorkerTrain;
// import common.WorkerTrain.Producer;
// import seq.manage.ReferenceGenome;
//
/// **
// * @author Kitty Shan entrop, for not sure
// */
// public class ShannonEntropy {
//
// private static ArrayList<String> getWords(ReferenceGenome referenceGenome, int wordSize, Logger
// log) {
// ArrayList<String> words = new ArrayList<String>(300000000);
// SAMSequenceDictionary samSequenceDictionary =
// referenceGenome.getIndexedFastaSequenceFile().getSequenceDictionary();
//
// for (SAMSequenceRecord samSequenceRecord : samSequenceDictionary.getSequences()) {
// int totalLen = 0;
// byte[] seq =
// referenceGenome.getIndexedFastaSequenceFile().getSequence(samSequenceRecord.getSequenceName()).getBases();
// while (totalLen < samSequenceRecord.getSequenceLength()) {
// String word = "";
// for (int i = totalLen; i < wordSize; i++) {
// try {
// word += new String(new byte[] { seq[i] }, "UTF-8").toUpperCase();
// } catch (UnsupportedEncodingException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// }
// words.add(word);
// totalLen++;
// }
// }
// return words;
// }
//
// private static class ShannonProducer implements Producer<Double> {
// private ReferenceGenome referenceGenome;
// private int[] wordSizes;
// private Logger log;
// private int index;
//
// public ShannonProducer(ReferenceGenome referenceGenome, int[] wordSizes, Logger log) {
// super();
// this.referenceGenome = referenceGenome;
// this.wordSizes = wordSizes;
// this.log = log;
// this.index = 0;
//
// }
//
// @Override
// public boolean hasNext() {
//
// // TODO Auto-generated method stub
// return index < wordSizes.length;
// }
//
// @Override
// public Callable<Double> next() {
// // TODO Auto-generated method stub
// ShannonWorker worker = new ShannonWorker(referenceGenome, wordSizes[index], log);
// index++;
// return worker;
//
// }
//
// @Override
// public void shutdown() {
// // TODO Auto-generated method stub
//
// }
//
// @Override
// public void remove() {
// // TODO Auto-generated method stub
//
// }
//
// }
//
// private static class ShannonWorker implements Callable<Double> {
// private ReferenceGenome referenceGenome;
// private int wordSize;
// private Logger log;
//
// public ShannonWorker(ReferenceGenome referenceGenome, int wordSize, Logger log) {
// super();
// this.referenceGenome = referenceGenome;
// this.wordSize = wordSize;
// this.log = log;
// }
//
// @Override
// public Double call() throws Exception {
// ArrayList<String> wordsG = getWords(referenceGenome, wordSize, log);
// return calculateShannonEntropy(wordsG);
// }
// }
//
// /**
// * pilfered from http://whaticode.com/2010/05/24/a-java-implementation-for-shannon-entropy/
// */
// private static Double calculateShannonEntropy(List<String> values) {
// Map<String, Integer> map = new HashMap<String, Integer>();
// // count the occurrences of each value
// for (String sequence : values) {
// if (!map.containsKey(sequence)) {
// map.put(sequence, 0);
// }
// map.put(sequence, map.get(sequence) + 1);
// }
//
// // calculate the entropy
// Double result = 0.0;
// for (String sequence : map.keySet()) {
// Double frequency = (double) map.get(sequence) / values.size();
// result -= frequency * (Math.log(frequency) / Math.log(2));
// }
// return result;
// }
//
//// private static double[] getShannonEntropies(ReferenceGenome referenceGenome, int[] wordSizes,
// int numThreads, Logger log) {
//// double[] se = new double[wordSizes.length];
//// ShannonProducer producer = new ShannonProducer(referenceGenome, wordSizes, log);
//// WorkerTrain<Double> train = new WorkerTrain<Double>(producer, numThreads, numThreads, log);
//// int index = 0;
//// while (train.hasNext()) {
//// se[index] = train.next();
////
//// index++;
//// }
////
//// return se;
////
//// }
//
// }
