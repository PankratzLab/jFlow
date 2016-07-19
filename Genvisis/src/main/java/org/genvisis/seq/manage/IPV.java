//package seq.manage;
//
//import java.io.File;
//import java.io.FileWriter;
//import java.io.PrintWriter;
//import java.util.ArrayList;
//
//import htsjdk.variant.variantcontext.VariantContext;
//import htsjdk.variant.vcf.VCFFileReader;
//import common.Array;
//import common.Logger;
//import common.ext;
//import seq.analysis.VCFSimpleTally;
//import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
//import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
//import seq.qc.FilterNGS.VariantContextFilter;
//
///**
// * @author lane0212
// *
// *         Generate input for IP-V, see what happens
// */
//public class IPV {
//	// print $convertedTsv "Variant Alleles"."\t"."Coding Start"."\t"."Validation Status"."\t"."Strand"."\t"."Gene Location(Transcript Strand)"."\t"."TranscriptID"."\t"."polyphen2"."\t"."sift"."\t"."rsID"."\t"."Variation Source"."\n";
//
//	private static final String[] HEADER = new String[] {
//			"Variant Alleles",
//			"Coding Start",
//			"Validation Status",
//			"Strand",
//			"Gene Location(Transcript Strand)",
//			"TranscriptID",
//			"polyphen2",
//			"sift",
//			"rsID",
//			"Variation Source",
//			"UCSC" };
//
//	private static final String[] ANNO = new String[] { "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_TRANSCRIPT_ID", "LJB23_Polyphen2_HDIV_score", "LJB23_SIFT_score", "snp138", "AAChange.refGene" };
//
//	private static void generate(String[] vcfs, String outputTSV, String[] geneNames, VariantContextFilter filter, Logger log) {
//		System.out.println("Don't use this");
//		System.exit(status);
//		String outDir = ext.parseDirectoryOfFile(outputTSV);
//		new File(outDir).mkdirs();
//
//		try {
//			PrintWriter writer = new PrintWriter(new FileWriter(outputTSV));
//			writer.println(Array.toStr(HEADER));
//			for (int i = 0; i < vcfs.length; i++) {
//
//				VCFFileReader reader = new VCFFileReader(vcfs[i], true);
//				ArrayList<String> ncbiAcc = new ArrayList<String>();
//				for (VariantContext vc : reader) {
//					if (filter.filter(vc).passed() && ext.indexOfStr(VCOps.getSNP_EFFGeneName(vc), geneNames) >= 0) {
//						if (vc.isBiallelic()) {
//							ArrayList<String> out = new ArrayList<String>();
//							String variantAlleles = vc.getReference().getBaseString() + "/" + vc.getAlternateAllele(0).getBaseString();
//							out.add(variantAlleles);
//							String[] annos = VCOps.getAnnotationsFor(ANNO, vc, ".");
//							int AAPOS = Integer.parseInt(annos[0].replaceAll("[^\\d.]", "")) * 3;
//							AAPOS = AAPOS - 3;
//							System.out.println(annos[annos.length - 1] + "");
//
//							String[] cc = annos[annos.length - 1].split(":");
//							String transcript = annos[1];
//							// for (int j = 0; j < cc.length; j++) {
//							// if(ext.indexOfStr(cc[j], array))
//							// }
//							// int add = 0;
//							// if (cc.length > 1) {
//							// for (int j = 0; j < cc[0].length(); j++) {
//							// if (!(cc[0].charAt(j) + "").equals(cc[1].charAt(j) + "")) {
//							// add = j;
//							// break;
//							// }
//							// }
//							// }
//							// AAPOS += add;
//							out.add("" + AAPOS);
//							out.add("N/A");
//							out.add("+");
//							out.add("+");
//							out.add(annos[1]);
//							out.add(annos[2]);
//							out.add(annos[3]);
//							out.add(annos[4]);
//
//							out.add("dbsnp");
//							out.add("miss");
//							out.add(VCOps.getSegment(vc).getUCSClocation());
//
//							writer.println(Array.toStr(Array.toStringArray(out)));
//
//						} else {
//							reader.close();
//							writer.close();
//							throw new IllegalArgumentException("Tri-allelic not supported");
//						}
//					}
//				}
//				reader.close();
//				break;
//			}
//			writer.close();
//		} catch (Exception e) {
//			log.reportError("Error writing to " + outputTSV);
//			log.reportException(e);
//		}
//
//	}
//
//	private static void test() {
//		String dir = "/panfs/roc/groups/5/pankrat2/lanej/test/IPV_TEST/";
//		String[] vcfs = new String[] { dir + "OSTEO_OFF_INHERIT.final.vcf.gz", dir + "OSTEO_OFF_INHERIT_CONTROL.final.vcf.gz" };
//		String dirOut = "/panfs/roc/groups/5/pankrat2/lanej/test/IPV_TEST/ipv/circos-p/Input/";
//
//		Logger log = new Logger(dir + ".log");
//		VARIANT_FILTER_BOOLEAN[] bFilters = new VARIANT_FILTER_BOOLEAN[] { VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER };
//		String[][] jexl = VCFSimpleTally.getG1000ESPFilter(0.01);
//		String[] genes = new String[] { "RECQL4" };
//		VariantContextFilter filter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, bFilters, jexl[0], jexl[1], log);
//		for (int i = 0; i < genes.length; i++) {
//			String geneout = dirOut + genes[i] + "/";
//			new File(geneout).mkdirs();
//			String output = geneout + genes[i] + "_input.txt";
//			generate(vcfs, output, new String[] { genes[i] }, filter, log);
//		}
//	}
//
//	public static void main(String[] args) {
//		test();
//	}
//
//}
