package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCOps;

import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class ConvertChargeToVCF {

	private static final String[] CHARGE_HEADER =
																							{	"SNP", "REF", "ALT", "SKATgene", "func_region",
																								"MAF_whites", "n_whites", "MAF_blacks", "n_blacks"};

	public static void processChargeMafFile(String fullPathToFile, Logger log) {
		if (!Files.headerOfFileContainsAll(fullPathToFile, CHARGE_HEADER, log)) {
			log.reportTimeError("This is designed for a specific file format with header "
													+ Array.toStr(CHARGE_HEADER));
		} else {
			VCFHeader vcfHeader = new VCFHeader();

			for (int i = 0; i < CHARGE_HEADER.length; i++) {
				VCFInfoHeaderLine vHeaderLine = new VCFInfoHeaderLine(CHARGE_HEADER[i], 1,
																															i <= 4	? VCFHeaderLineType.String
																																			: VCFHeaderLineType.Float,
																															CHARGE_HEADER[i] + " from Charge");
				vcfHeader.addMetaDataLine(vHeaderLine);
			}

			int countSkipped = 0;
			String outputVCF = fullPathToFile + ".vcf";
			ArrayList<Segment> segs = new ArrayList<Segment>(12000000);
			ArrayList<VariantContext> vcs = new ArrayList<VariantContext>(12000000);

			VariantContextWriterBuilder builder =
																					new VariantContextWriterBuilder().setOutputFile(outputVCF);
			builder.clearOptions();

			builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);
			VariantContextWriter writer = builder.build();
			writer.writeHeader(vcfHeader);
			try {
				int[] indices = ext.indexFactors(	CHARGE_HEADER, Files.getHeaderOfFile(fullPathToFile, log),
																					true, false);
				BufferedReader reader = Files.getAppropriateReader(fullPathToFile);
				reader.readLine();
				int count = 0;
				while (reader.ready()) {
					count++;
					if (count % 100000 == 0) {
						log.reportTimeInfo(count + "");
					}
					String[] line = reader.readLine().trim().split("[\\s]+");
					String ID = line[indices[0]];
					String[] tmp = ID.split(":");
					String chr = Positions.getChromosomeUCSC(Positions.chromosomeNumber(tmp[0]), true);
					int pos = Integer.parseInt(tmp[1]);
					String ref = line[indices[1]];
					String alt = line[indices[2]];
					if ((!ref.equals(".") || !alt.equals(".")) && !alt.contains(",")) {
						ArrayList<Allele> alleles = new ArrayList<Allele>();
						if (!ref.equals(".")) {
							alleles.add(Allele.create(ref, true));
						}
						alleles.add(Allele.create(alt, false));

						// alleles.add(alt);
						VariantContextBuilder vBuilder = new VariantContextBuilder();
						vBuilder.chr(chr);
						vBuilder.start(pos);
						vBuilder.stop(pos + ref.length() - 1);
						vBuilder.id(ID);
						vBuilder.alleles(alleles);
						// vBuilder.filter("PASS");
						for (int i = 0; i < CHARGE_HEADER.length; i++) {
							vBuilder.attribute(	CHARGE_HEADER[i],
																	i <= 4	? line[indices[i]]
																					: Float.parseFloat(line[indices[i]]) + "");
						}
						VariantContext vc = vBuilder.make();
						segs.add(VCOps.getSegment(vc));
						vcs.add(vc);
						// writer.add(vBuilder.make());
					} else {
						countSkipped++;
					}

				}

				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + fullPathToFile + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + fullPathToFile + "\"");
				return;
			}

			log.reportTimeInfo("Sorting and writing variants");

			Segment[] segsA = segs.toArray(new Segment[segs.size()]);
			Segment[] segsASorted = segs.toArray(new Segment[segs.size()]);
			Arrays.sort(segsASorted);
			log.reportTimeInfo("Finished sorting");

			int[] order = new int[segsASorted.length];
			Arrays.fill(order, -1);
			for (int i = 0; i < segsA.length; i++) {
				int[] indices = Segment.binarySearchForAllOverLappingIndices(segsA[i], segsASorted);
				boolean found = false;
				for (int indice : indices) {
					if (segsA[i].equals(segsASorted[indice])) {
						order[i] = indice;
						if (found) {
							log.reportTimeError("multiple founds");
						}
						found = true;
					}
				}
				if (!found) {
					log.reportTimeError("Could not find");
				}

			}
			for (int element : order) {
				writer.add(vcs.get(element));
			}
			writer.close();
			if (countSkipped > 0) {
				log.reportTimeWarning("Skipped "	+ countSkipped
															+ " variants due to ambigous alt alleles, or to lack of ref/alt alleles");
			}

			try {
				File indexFile = Tribble.indexFile(new File(outputVCF));

				Index index = IndexFactory.createLinearIndex(new File(outputVCF), new VCFCodec());
				LittleEndianOutputStream stream =
																				new LittleEndianOutputStream(new FileOutputStream(indexFile));
				index.write(stream);
				stream.close();
			} catch (IOException e) {
				log.reportError("Error - could not create index file for" + fullPathToFile + ".vcf.gz");
			}
			// VCFFileReader reader = new VCFFileReader(fullPathToFile + ".vcf.gz", true);
			// for (VariantContext vc : reader) {
			// System.out.println("hi");
			// }
			// VCFHeaderLine info = new VCFHeaderLine(key, value)
			// vcfHeader.addMetaDataLine(headerLine);
		}

	}

	public static void convert(String fullPathToFile) {
		processChargeMafFile(fullPathToFile, new Logger());
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "ConvertChargeToVCF.dat";
		String fullPathToFile =
													"D:/data/Project_Tsai_Project_021/Variants/charge_fibrinogen_mafs_and_counts.xln";
		// String fullPathToFile = null;
		// String logfile = null;

		String usage = "\n"	+ "jlDev.ConvertChargeToVCF requires 0-1 arguments\n"
										+ "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("file=")) {
				fullPathToFile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			convert(fullPathToFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
