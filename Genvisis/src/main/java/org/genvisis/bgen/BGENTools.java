package org.genvisis.bgen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeSet;
import java.util.zip.DataFormatException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.Inflater;
import java.util.zip.ZipOutputStream;

import org.apache.commons.math3.util.CombinatoricsUtils;
import org.genvisis.bgen.BGENReader.BGENRecord;
import org.genvisis.bgen.BGENReader.BGENRecordMetaData;
import org.genvisis.bgen.BGENReader.COMPRESSION;

public final class BGENTools {

	private BGENTools() {}

	public static double[][] readLayout1Record(boolean skip, RandomAccessFile in,
																						 BGENRecordMetaData r, COMPRESSION c)
																																								 throws IOException {
		byte[] read = new byte[4];

		if (skip) {
			in.skipBytes((int) r.blockLength);
			return null;
		}

		read = new byte[(int) r.blockLength];
		in.read(read);
		byte[] output;

		if (c == COMPRESSION.ZLIB) {
			Inflater inflater = new Inflater();
			inflater.setInput(read);
			ByteArrayOutputStream outputStream = new ByteArrayOutputStream(read.length);
			byte[] buffer = new byte[(int) (6 * r.N)];
			while (!inflater.finished()) {
				int count;
				try {
					count = inflater.inflate(buffer);
				} catch (DataFormatException e) {
					throw new RuntimeException(e);
				}
				outputStream.write(buffer, 0, count);
			}
			outputStream.close();
			output = outputStream.toByteArray();
		} else {
			output = read;
		}

		double[][] data = new double[output.length / 6][];
		int n;
		for (int p = 0; p < output.length; p += 6) {
			n = p / 6;
			byte[] prob1B = {output[p], output[p + 1]};
			byte[] prob2B = {output[p + 2], output[p + 3]};
			byte[] prob3B = {output[p + 4], output[p + 5]};
			int prob1, prob2, prob3;
			prob1 = BGENBitMath.unsignedShortToInt(prob1B, true);
			prob2 = BGENBitMath.unsignedShortToInt(prob2B, true);
			prob3 = BGENBitMath.unsignedShortToInt(prob3B, true);
			data[n] = new double[] {prob1 / (double) 32768, prob2 / (double) 32768,
															prob3 / (double) 32768};
		}

		return data;
	}

	public static double[][] readLayout2Record(boolean skip, RandomAccessFile in,
																						 BGENRecordMetaData r, COMPRESSION c)
																																								 throws IOException {
		byte[] read = new byte[4];
		long uncompLength;

		if (skip) {
			in.skipBytes((int) r.blockLength);
			return null;
		}

		if (c == COMPRESSION.NONE) {
			uncompLength = 6 * r.N;
		} else {
			in.read(read);
			uncompLength = BGENBitMath.unsignedIntToLong(read, true);
		}

		if (r.blockLength > Integer.MAX_VALUE) {
			System.err.println("Compressed data > Integer MAX - you're going to encounter problems.");
		}

		if (uncompLength > Integer.MAX_VALUE) {
			System.err.println("Uncompressed data > Integer MAX - you're going to encounter problems.");
		}

		byte[] block = new byte[(int) r.blockLength - (c == COMPRESSION.NONE ? 0 : 4)];
		byte[] output;
		in.read(block);

		switch (c) {
			case ZLIB:
				Inflater inflater = new Inflater();
				inflater.setInput(block);
				ByteArrayOutputStream outputStream = new ByteArrayOutputStream(block.length);
				byte[] buffer = new byte[(int) uncompLength];
				while (!inflater.finished()) {
					int count;
					try {
						count = inflater.inflate(buffer);
					} catch (DataFormatException e) {
						throw new RuntimeException(e);
					}
					outputStream.write(buffer, 0, count);
				}
				outputStream.close();
				output = outputStream.toByteArray();
				break;
			case NONE:
				output = block;
				break;
			case ZSTD:
			default:
				throw new UnsupportedOperationException("ZSTD compression is not supported.");
		}

		block = new byte[4];
		System.arraycopy(output, 0, block, 0, 4); // 0-3
		long sampleCount = BGENBitMath.unsignedIntToLong(block, true);

		block = new byte[2];
		System.arraycopy(output, 4, block, 0, 2); // 4,5
		int alleleCount = BGENBitMath.unsignedShortToInt(block, true);

		// int ploidyMin = output[6];
		// int ploidyMax = output[7];

		byte[] samplePloidy = new byte[(int) sampleCount];
		System.arraycopy(output, 8, samplePloidy, 0, (int) sampleCount);
		boolean[] missing = new boolean[samplePloidy.length];
		int[] ploidyVals = new int[samplePloidy.length];
		/*-
		 * Missingness is encoded in the most-significant bit of the eight bit byte.
		 * 
		 * Therefore, if the value is > 127, the most-significant bit has been set, and the sample is
		 * set to missing.
		 * 
		 * Ploidy is encoded in the least-significant 6 bits (leaving one bit unused between missingness
		 * and ploidy). We account for the undefined bit (if set, subtract) and then use the remaining
		 * value as the ploidy value.
		 */
		int temp;
		for (int i = 0; i < samplePloidy.length; i++) {
			temp = samplePloidy[i];
			missing[i] = false;
			if (temp >= 128) {
				missing[i] = true;
				temp -= 128;
			}
			ploidyVals[i] = temp > 63 ? temp - 64 : temp;
		}

		int phasedFlag = output[(int) (8 + sampleCount)]; // 1 == haplotypes, 0 == genotypes
		int probBits = output[(int) (9 + sampleCount)];
		int probMax = 0;
		for (int i = 0; i < probBits; i++) {
			probMax |= 1;
			if (i < probBits - 1) {
				probMax <<= 1;
			}
		}

		double[][] data = new double[(int) sampleCount][];

		byte readByte;
		int K = alleleCount;
		int byteStart = (int) (10 + sampleCount);
		int bitReadPtr = 0;

		for (int n = 0; n < sampleCount; n++) {
			boolean[] sampProb = new boolean[probBits];
			/*-
			 * let Z = ploidy[n]
			 * let K = N_alleles
			 * 
			 * NProbs = {
			 * if (phased),
			 * 					Z(K-1) 
			 * else,
			 * 					(Z+K-1)!
			 * 				  -------- 
			 * 					(K-1)!Z!
			 * }
			 */
			int Z = ploidyVals[n];
			int nProbsMin1;
			if (phasedFlag == 1) {
				nProbsMin1 = Z * (K - 1);
			} else {
				nProbsMin1 = (int) ((CombinatoricsUtils.factorial(Z + K - 1))
										 / ((CombinatoricsUtils.factorial(K - 1)) * (CombinatoricsUtils.factorial(Z))));
				nProbsMin1 -= 1;
			}
			// last prob (or last prob per haplotype for phased data) is recorded as 1 - sum(probs)
			double[] probs = new double[nProbsMin1 + (phasedFlag == 0 ? 1 : Z)];
			double sum = 0;
			for (int i = 0; i < nProbsMin1; i++) {
				for (int p = 0; p < probBits; p++) {
					readByte = output[byteStart];
					sampProb[p] = BGENBitMath.getBit(readByte, bitReadPtr);
					bitReadPtr++;
					if (bitReadPtr == 8) {
						byteStart++;
						bitReadPtr = 0;
					}
				}
				probs[i] = readLayout2Probs(sampProb, probMax);
				sum += probs[i];
				if (phasedFlag == 1) {
					// untested :(
					if ((i % K - 1) == K - 2) {
						i++;
						probs[i] = 1 - sum;
						sum = 0;
					}
				}
			}
			if (phasedFlag == 0) {
				probs[probs.length - 1] = 1 - sum;
			}

			data[n] = probs;
		}

		return data;
	}

	private static double readLayout2Probs(boolean[] probBits, int probMax) {
		int probVal = 0;
		for (int i = 0; i < probBits.length; i++) {
			probVal |= probBits[i] ? 1 : 0;
			if (i < probBits.length - 1) {
				probVal <<= 1;
			}
		}
		return ((double) (probVal)) / (double) probMax;
	}

	/**
	 * Reads a serialized genomically-ordered chromosome-mapped variant metadata Map object.
	 * 
	 * @param reader
	 * @param mapFile
	 * @throws ClassNotFoundException
	 * @throws IOException
	 */
	@SuppressWarnings("unchecked")
	public static void loadMapInfo(BGENReader reader, String mapFile) throws ClassNotFoundException,
																																	 IOException {
		InputStream in;
		ObjectInputStream ois = null;

		if (mapFile.endsWith(".gz")) {
			in = new BufferedInputStream(new GZIPInputStream(new FileInputStream(mapFile)));
		} else {
			in = new BufferedInputStream(new FileInputStream(mapFile));
		}
		ois = new ObjectInputStream(in);
		reader.chrSets = (Map<Integer, TreeSet<BGENRecordMetaData>>) ois.readObject();
		ois.close();
	}

	/**
	 * Serializes the genomically-ordered chromosome-mapped variant metadata Map to a file.
	 * 
	 * @param reader
	 * @param mapFile
	 * @throws IOException
	 */
	public static void serializeMapInfo(BGENReader reader, String mapFile) throws IOException {
		if (reader.chrSets.isEmpty()) {
			System.err.println("Error - no map information is available, nothing to do.");
			return;
		}
		ObjectOutputStream oos;
		File f = new File(mapFile);
		if (f.exists() || !f.mkdirs()) {
			System.err.println("Error - No valid path to file: " + mapFile);
			return;
		}

		if (mapFile.endsWith(".gz")) {
			oos = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(mapFile)));
		} else {
			oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(mapFile)));
		}
		oos.writeObject(reader.chrSets);
		oos.flush();
		oos.close();
	}

	public static void writeMapInfoFile(BGENReader reader, String mapFile) throws IOException {
		if (reader.chrSets.isEmpty()) {
			System.err.println("Error - no map information has been stored; written map file will be exported in current order.");
		}
		PrintWriter writer;
		if (mapFile.endsWith(".gz")) {
			writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(mapFile, false)));
		} else if (mapFile.endsWith(".zip")) {
			writer = new PrintWriter(new ZipOutputStream(new FileOutputStream(mapFile, false)));
		} else {
			writer = new PrintWriter(new BufferedWriter(new FileWriter(mapFile, false)));
		}
		writer.println("CHR\tPOS\tID\tA1\tA2\tLOC");
		if (reader.chrSets.isEmpty()) {
			reader.reset();
			for (BGENRecord rec : new BGENIterators.BGENIterable(new BGENIterators.BGENIterator(reader,
																																													false))) {
				writeBGENRecHead(rec.getMetaData(), reader, writer);
			}
		} else {
			for (Entry<Integer, TreeSet<BGENRecordMetaData>> entry : reader.chrSets.entrySet()) {
				for (BGENRecordMetaData rec : entry.getValue()) {
					writeBGENRecHead(rec, reader, writer);
				}
			}
		}
		writer.flush();
		writer.close();
	}


	private static void writeBGENRecHead(BGENRecordMetaData record, BGENReader bgenReader,
																			 PrintWriter writer) {
		if (record.getAlleles().length > 2) {
			System.err.print("Error - variant " + record.getID() + "," + record.getRsID()
											 + " has >2 alleles: {");
			for (int i = 0; i < record.getAlleles().length; i++) {
				System.err.println(record.getAlleles()[i]);
				if (i < record.getAlleles().length - 1) {
					System.err.println(",");
				}
				System.err.println("}");
			}
		}
		writer.println(record.getChr() + "\t" + record.getPos() + "\t" + record.getID() + ","
									 + record.getRsID() + "\t"
									 + record.getAlleles()[0] + "\t" + record.getAlleles()[1] + "\t" + record.ptrByt);
	}

	public static List<String> getSortedChrPos(BGENReader reader) {
		List<String> values = new ArrayList<>();
		for (Entry<Integer, TreeSet<BGENRecordMetaData>> entry : reader.chrSets.entrySet()) {
			for (BGENRecordMetaData rec : entry.getValue()) {
				values.add(rec.getChr() + ":" + rec.getPos());
			}
		}
		return values;
	}

}
