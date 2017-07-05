/**
 * 
 */
package org.genvisis.seq.analysis.mutMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import org.genvisis.seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import org.genvisis.seq.qc.FilterNGS.VariantContextFilter;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 *
 * Class to format {@link VariantContext} to be used as input to
 * http://www.cbioportal.org/mutation_mapper.jsp
 */
public class VCMut {

	private VariantContext vc;
	private String hugoSymbol;// HUGO symbol for the gene
	private String proteinChange;// Amino acid change V600E
	private String mutationType;// Translational effect of variant allele Missense_Mutation,
															// Nonsense_Mutation, etc.
	private String chromosome;// Chromosome number X, Y, M, 1, 2, etc.
	private int startPosition;// Lowest numeric position of the reported variant on the genomic
														// reference sequence 666
	private int endPosition;// Highest numeric position of the reported variant on the genomic
													// reference sequence 667
	private Allele referenceAllele;// The plus strand reference allele at this position A
	// private List<Allele> variantAlleles; // The plus strand reference allele at this position A


	private boolean hasFullInfo;



	/**
	 * See headers at http://www.cbioportal.org/mutation_mapper.jsp
	 * 
	 * @param vc
	 * @param hugoSymbol
	 * @param proteinChange
	 * @param mutationType
	 * @param chromosome
	 * @param startPosition
	 * @param endPosition
	 * @param referenceAllele
	 * @param hasFullInfo has hugo symbol, etc
	 * 
	 */
	public VCMut(VariantContext vc, String hugoSymbol, String proteinChange, String mutationType,
							 String chromosome, int startPosition, int endPosition, Allele referenceAllele,
							 boolean hasFullInfo) {
		super();
		this.vc = vc;
		this.hugoSymbol = hugoSymbol;
		this.proteinChange = proteinChange;
		this.mutationType = mutationType;
		this.chromosome = chromosome;
		this.startPosition = startPosition;
		this.endPosition = endPosition;
		this.referenceAllele = referenceAllele;
		this.hasFullInfo = hasFullInfo;
	}

	/**
	 * @param other
	 * @param includeVC copy the {@link VariantContext}
	 */
	public VCMut(VCMut other, boolean includeVC) {
		super();
		this.vc = includeVC ? other.vc : null;
		this.hugoSymbol = other.hugoSymbol;
		this.proteinChange = other.proteinChange;
		this.mutationType = other.mutationType;
		this.chromosome = other.chromosome;
		this.startPosition = other.startPosition;
		this.endPosition = other.endPosition;
		this.referenceAllele = other.referenceAllele;
		this.hasFullInfo = other.hasFullInfo;
	}



	public VariantContext getVc() {
		return vc;
	}



	public String getHugoSymbol() {
		return hugoSymbol;
	}



	public String getProteinChange() {
		return proteinChange;
	}



	public String getMutationType() {
		return mutationType;
	}



	public String getChromosome() {
		return chromosome;
	}



	public int getStartPosition() {
		return startPosition;
	}



	public int getEndPosition() {
		return endPosition;
	}



	public Allele getReferenceAllele() {
		return referenceAllele;
	}



	public boolean isHasFullInfo() {
		return hasFullInfo;
	}


	/**
	 * Where to extract info from the {@link VariantContext}
	 *
	 */
	public enum PARSE_METHOD {

		SNP_EFF("SNPEFF_GENE_NAME", "SNPEFF_AMINO_ACID_CHANGE", "SNPEFF_EFFECT", "SNPEFF_IMPACT");

		/**
		 * Gene change annotation
		 */
		private String geneFlag;
		/**
		 * AA change annotation
		 */
		private String aaFlag;
		/**
		 * mutation_Type annotation
		 */
		private String effFlag;

		private String impactFlag;


		private PARSE_METHOD(String geneFlag, String aaFlag, String effFlag, String impactFlag) {
			this.geneFlag = geneFlag;
			this.aaFlag = aaFlag;
			this.effFlag = effFlag;
			this.impactFlag = impactFlag;
		}
	}


	/**
	 * By sample
	 *
	 */
	public static class MutInd {
		private String sampleID;
		private Allele variantAllele;
		private VCMut mut;

		/**
		 * @param sampleID
		 * @param variantAllele allele this sample has
		 * @param mut {@link VCMut} for variant level annotation
		 */
		public MutInd(String sampleID, Allele variantAllele, VCMut mut) {
			super();
			this.sampleID = sampleID;
			this.variantAllele = variantAllele;
			this.mut = mut;
		}


		public List<String> getMutIndMapFormat() {
			List<String> data = new ArrayList<>();
			data.add(sampleID);
			data.add(variantAllele.getDisplayString());
			data.addAll(mut.getMutMapFormat());

			return data;
		}

		public static List<String> getMutIndMapHeader() {
			List<String> data = new ArrayList<>();
			data.add("Sample_ID");
			data.add("Variant_Allele");
			data.addAll(getMutMapHeader());
			return data;
		}


	}

	// Hugo_Symbol Sample_ID Protein_Change Mutation_Type Chromosome Start_Position End_Position
	// Reference_Allele Variant_Allele

	private static List<String> getMutMapHeader() {
		List<String> data = new ArrayList<>();
		data.add("Hugo_Symbol");
		data.add("Protein_Change");
		data.add("Mutation_Type");
		data.add("Chromosome");
		data.add("Start_Position");
		data.add("End_Position");
		data.add("Reference_Allele");
		return data;
	}

	private List<String> getMutMapFormat() {
		List<String> data = new ArrayList<>();
		data.add(hugoSymbol);
		data.add(proteinChange);
		data.add(mutationType);
		data.add(Positions.getChromosomeUCSC(Positions.chromosomeNumber(chromosome),
																				 false));
		data.add(Integer.toString(startPosition));
		data.add(Integer.toString(endPosition));
		data.add(referenceAllele.getDisplayString());

		return data;
	}

	/**
	 * @param inds individuals to use, or an empty set for all
	 * @param filter filters to be applied to the variant
	 * @param gq genotype quality
	 * @return
	 */
	public List<MutInd> parseToInds(Set<String> inds, VariantContextFilter filter, int gq) {
		List<MutInd> mutInds = new ArrayList<>();
		if (filter.filter(vc).passed()) {
			GenotypesContext gc = vc.getGenotypes();
			for (Genotype g : gc) {
				if (inds.contains(g.getSampleName()) && g.isCalled() && !g.isHomRef()) {
					if (g.getGQ() > gq) {
						List<Allele> als = g.getAlleles();
						for (Allele a : als) {
							if (a.isNonReference()) {
								MutInd mutInd = new MutInd(g.getSampleName(), a, new VCMut(this, false));
								mutInds.add(mutInd);
								break;
							}
						}
					}
				}
			}
		}
		return mutInds;
	}


	/**
	 * @param vc
	 * @param parseMethod
	 * @param mutMap
	 * @param log
	 * @return a {@link VCMut} parsed from {@link VariantContext} according to {@link PARSE_METHOD}
	 */
	public static VCMut fromVariantContext(VariantContext vc, PARSE_METHOD parseMethod,
																				 Map<String, String> mutMap) {
		String chromosome = vc.getContig();
		int startPosition = vc.getStart();
		int endPosition = vc.getEnd();
		Allele referenceAllele = vc.getReference();

		String hugoSymbol = VCOps.getAnnotationsFor(new String[] {parseMethod.geneFlag}, vc,
																								".")[0];
		String proteinChange = VCOps.getAnnotationsFor(new String[] {parseMethod.aaFlag}, vc,
																									 ".")[0];
		String mutationAnn = VCOps.getAnnotationsFor(new String[] {parseMethod.effFlag}, vc,
																								 ".")[0];
		String mutationImpact = VCOps.getAnnotationsFor(new String[] {parseMethod.impactFlag}, vc,
																										".")[0];
		boolean hasFullInfo = false;
		String mutationType = "NA";
		if (mutMap.containsKey(mutationAnn)) {
			hasFullInfo = true;
			mutationType = mutMap.get(mutationAnn);
		} else if (!".".equals(mutationImpact)) {
			throw new IllegalStateException("should really be able to parse " + mutationAnn + "\t"
																			+ mutationImpact);
		}


		return new VCMut(vc, hugoSymbol, proteinChange, mutationType, chromosome, startPosition,
										 endPosition, referenceAllele, hasFullInfo);
	}



	/**
	 * @param log
	 * @return the SNPEFF EFFECT to Sequence Ontology mapping
	 */
	public static Map<String, String> getSNPEFFOntologyMap(Logger log) {



		HashMap<String, String> effMap = new HashMap<>();

		BufferedReader br;
		try {
			br = new BufferedReader(new InputStreamReader(loadMutMap()));
			String sCurrentLine;

			while ((sCurrentLine = br.readLine()) != null) {
				if (!sCurrentLine.startsWith("#")) {
					String[] info = sCurrentLine.trim().split("\t");
					effMap.put(info[1], info[0]);
				}
			}
			br.close();
		} catch (IOException e) {
			log.reportException(e);
		}

		return effMap;

	}

	private static InputStream loadMutMap() {
		return VCMut.class.getResourceAsStream("SNPEFF_MutMap.txt");
		// return new File("/home/tsaim/lane0212/SNPEFF_MutMap.txt");
	}



	@Override
	public String toString() {
		return "VCMut [ hugoSymbol=" + hugoSymbol + ", proteinChange=" + proteinChange
					 + ", mutationType=" + mutationType + ", chromosome=" + chromosome + ", startPosition="
					 + startPosition + ", endPosition=" + endPosition + ", referenceAllele=" + referenceAllele
					 + ", hasFullInfo=" + hasFullInfo + "]";
	}


	public static void run(String vcf, String outputDir, String vpopFile,
												 double maf, boolean removeFilter, Segment seg) {
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "mutMap.log");

		VcfPopulation vpop = VcfPopulation.load(vpopFile,
																						POPULATION_TYPE.ANY, log);
		vpop.report();
		VariantContextFilter filter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {},
																													 removeFilter ? new VARIANT_FILTER_BOOLEAN[] {VARIANT_FILTER_BOOLEAN.FAILURE_FILTER}
																																				: new VARIANT_FILTER_BOOLEAN[] {},
																													 new String[] {"Freq"},
																													 new String[] {FilterNGS.getPopFreqFilterString(maf)},
																													 log);

		Map<String, String> effmap = VCMut.getSNPEFFOntologyMap(new Logger());


		for (String pop : vpop.getSuperPop().keySet()) {

			String output = outputDir + VCFOps.getAppropriateRoot(vcf, true) + "_pop_" + pop + "_maf_"
											+ maf
											+ "_filt_" + removeFilter + ".txt";
			log.reportTimeInfo("converting population " + pop + ", reporting to " + output);
			PrintWriter writer = Files.getAppropriateWriter(output);
			writer.println(ArrayUtils.toStr(MutInd.getMutIndMapHeader()));
			VCFFileReader reader = new VCFFileReader(new File(vcf));
			CloseableIterator<VariantContext> iter = reader.query(seg.getChromosomeUCSC(),
																														seg.getStart(),
																														seg.getStop());
			while (iter.hasNext()) {
				VariantContext vc = iter.next();
				VCMut vcMut = VCMut.fromVariantContext(vc, PARSE_METHOD.SNP_EFF, effmap);
				if (vcMut.isHasFullInfo()) {
					List<MutInd> inds = vcMut.parseToInds(vpop.getSuperPop().get(pop), filter, -4);
					for (MutInd ind : inds) {
						writer.println(ArrayUtils.toStr(ind.getMutIndMapFormat()));
					}
				}

			}
			iter.close();
			reader.close();
			writer.close();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CLI c = new CLI(VCMut.class);
		c.addArgWithDefault(CLI.ARG_VCF, CLI.DESC_VCF, "a.vcf.gz");
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, "out/");
		c.addArgWithDefault("vpop", "VCF population file", "a.vpop");
		c.addArgWithDefault("maf", "maf threshold (from 1000G ExAC, etc), set to 1.2 for no filter",
												"0.01");
		c.addFlag("noFilt", "do not remove filtered variants (FILTER field)");
		c.addArgWithDefault("seg", "segment to limit search space,", "chr11:108056992-108276393");
		c.parseWithExit(args);

		run(c.get(CLI.ARG_VCF), c.get(CLI.ARG_OUTDIR), c.get("vpop"), c.getD("maf"), !c.has("noFilt"),
				new Segment(c.get("seg")));


	}
}
