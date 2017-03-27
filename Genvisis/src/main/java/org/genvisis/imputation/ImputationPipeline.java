package org.genvisis.imputation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.Qc;
import org.genvisis.seq.manage.ReferenceGenome;

import com.google.common.collect.ImmutableSet;


public class ImputationPipeline {
	
	Project proj;
	Set<String> dropMarkers = new HashSet<String>();
	Set<String> dropSamples = new HashSet<String>();
	Map<String, Marker> prepMarkers = new HashMap<String, Marker>();
	Set<String> prepMarkersNames = new HashSet<String>();
	
	Set<String> markersToExport;
	
	public ImputationPipeline(Project proj, String referenceFile) {
		this.proj = proj;
		ImputationPrep prep = new ImputationPrep(proj, referenceFile);
		Set<Marker> markers = prep.getMatchingMarkers();
		for (Marker m : markers) {
			prepMarkersNames.add(m.getName());
			prepMarkers.put(m.getName(), m);
		}
	}
	
	public void loadDefaultDropFiles(String plinkDir) {
		String dir = plinkDir + Qc.QC_SUBDIR + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
		String mark = dir + FurtherAnalysisQc.MARKER_QC_DROPS;
		String samp = dir + FurtherAnalysisQc.SAMPLE_QC_DROPS;
		setMarkersToDropFile(mark);
		setSamplesToDropFile(samp);
	}
	
	public void setSamplesToDropFile(String samplesToDropFile) {
		if (!Files.exists(samplesToDropFile)) {
			proj.getLog().reportTimeWarning("Sample drop file doesn't exist: " + samplesToDropFile);
			return;
		}
		dropSamples = HashVec.loadFileToHashSet(samplesToDropFile, false);
	}

	public void setMarkersToDropFile(String markersToDropFile) {
		if (!Files.exists(markersToDropFile)) {
			proj.getLog().reportTimeWarning("Marker drop file doesn't exist: " + markersToDropFile);
			return;
		}
		dropMarkers = HashVec.loadFileToHashSet(markersToDropFile, false);
	}
	
	private Set<String> getMarkersToExport() {
		if (markersToExport == null) {
  		markersToExport = new HashSet<String>(Arrays.asList(proj.getMarkerNames()));
  		markersToExport.removeAll(dropMarkers);
  		markersToExport.retainAll(prepMarkersNames);
		}
		return markersToExport;
	}
	
	private Set<String> getChrMarkers(int chr) {
		Set<String> markersToExport = getMarkersToExport();
		Set<String> chrMarkers = new HashSet<String>();
		for (String m : markersToExport) {
			if (prepMarkers.get(m).getChr() == (byte) chr) {
				chrMarkers.add(m);
			}
		}
		return chrMarkers;
	}
	
	private ArrayList<String> getMarkersSortedNoDupes(int chr) {
		Set<String> chrMarkers = getChrMarkers(chr);

		int[] pos = new int[chrMarkers.size()];
		String[] mkr = chrMarkers.toArray(new String[chrMarkers.size()]);
		for (int i = 0; i < mkr.length; i++) {
			pos[i] = prepMarkers.get(mkr[i]).getPosition();
		}
		
		int[] indices = Sort.getSortedIndices(pos);
		
		ArrayList<String> mkrs = new ArrayList<String>();
		for (int i = 0; i < indices.length; i++) {
			if (i == 0 || pos[indices[i]] != pos[indices[i-1]]) { // skip if prev (in sorted array) was same position
				mkrs.add(mkr[indices[i]]);
			}
		}
		return mkrs;
	}

	public void exportToPlink(String plinkDirAndRoot) {
		// TODO (??) Only alphanumeric characters in FID/IID
		String[] writtenDNAs = PlinkData.createFamFile(proj, plinkDirAndRoot, dropSamples);
		if (writtenDNAs == null) {
			// TODO error
			return;
		}
		int[] indicesOfTargetSamplesInProj = PlinkData.getIndicesOfTargetSamplesInProj(proj, writtenDNAs, proj.getLog());
		
		// UNUSED - could potentially apply
		String clusterFilterFileName = null;
		
		float gcThreshold = 0; /* this is unused in plink export code */ // proj.GC_THRESHOLD.getValue().floatValue();
		
		for (int chr = 1; chr < 23; chr++) {
			ArrayList<String> mkrs = getMarkersSortedNoDupes(chr);
			
  		String[] targetMarkers = mkrs.toArray(new String[mkrs.size()]);
  		int[] indicesOfTargetMarkersInProj = new int[targetMarkers.length];
  		HashMap<String, Byte> chrsOfTargetMarkers = new HashMap<String, Byte>();
  		HashMap<String, Integer> posOfTargetMarkers = new HashMap<String, Integer>();
  		PlinkData.getIndicesOfTargetMarkers(proj, targetMarkers, indicesOfTargetMarkersInProj,
  															chrsOfTargetMarkers, posOfTargetMarkers, proj.getLog());
  		
  		String dirAndRoot = plinkDirAndRoot + "_chr" + chr;
  		boolean success = PlinkData.createBedFileSnpMajor10KperCycle(proj, ImmutableSet.copyOf(targetMarkers),
  																											 chrsOfTargetMarkers, posOfTargetMarkers,
  																											 indicesOfTargetSamplesInProj,
  																											 clusterFilterFileName, gcThreshold,
  																											 dirAndRoot, proj.getLog());
  		
  		if (success) {
  			PrintWriter refWriter = Files.getAppropriateWriter(dirAndRoot + "_alleles.ref");
				for (String s : targetMarkers) {
					refWriter.println(s + "\t" + prepMarkers.get(s).getRef());	
				}
  			refWriter.flush();
  			refWriter.close();
  			Files.copyFile(plinkDirAndRoot + ".fam", dirAndRoot + ".fam");
  		}
		}
		
	}
	
	public void exportToVCF(String vcfDirAndRoot) {
		SampleData sd = proj.getSampleData(0, false);
		String[] allSamples = proj.getSamples();
		List<String> idsToInclude = new ArrayList<String>();
		Map<String, Integer> idIndexMap = new HashMap<String, Integer>();
		for (int i = 0; i < allSamples.length; i++) {
			String s = allSamples[i];
			if (!dropSamples.contains(sd.lookup(s)[1])) {
				idsToInclude.add(s);
			}
			idIndexMap.put(s, i);
		}

		// UNUSED - could potentially apply
		ClusterFilterCollection clusterFilters = null;
		
		float gcThreshold = 0; /* this is unused in plink export code, so we won't use it here, either */ //proj.GC_THRESHOLD.getValue().floatValue();
		
		
		for (int chr = 1; chr < 23; chr++) {
			ArrayList<String> mkrs = getMarkersSortedNoDupes(chr);
			
			String fileOut = vcfDirAndRoot + "_chr" + chr + ".vcf.gz";
			
			VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(fileOut);
			builder.clearOptions();
			builder.setOption(Options.INDEX_ON_THE_FLY);
			
			HashSet<VCFHeaderLine> lines = new HashSet<VCFHeaderLine>();
			VCFFormatHeaderLine format = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "GT");
			lines.add(format);
			
			VCFHeader vcfHeader = new VCFHeader(lines, idsToInclude);

			SAMSequenceDictionary samSequenceDictionary = new ReferenceGenome(Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog()).getFASTA().getAbsolute(), proj.getLog())
																															.getIndexedFastaSequenceFile().getSequenceDictionary();

			builder.setReferenceDictionary(samSequenceDictionary);
			vcfHeader.setSequenceDictionary(samSequenceDictionary);
			
			VariantContextWriter writer = builder.build();
			vcfHeader.hasGenotypingData();
			writer.writeHeader(vcfHeader);

			
			MDL mdl = new MDL(proj, proj.getMarkerSet(), mkrs.toArray(new String[mkrs.size()]));

			while (mdl.hasNext()) {
				MarkerData markerData = mdl.next();
				
				VariantContextBuilder builderVc = new VariantContextBuilder();
				builderVc.chr("chr" + chr);
				Marker mkr = prepMarkers.get(markerData.getMarkerName());
				ArrayList<Allele> a = new ArrayList<Allele>();
				Allele aR = Allele.create(Character.toString(mkr.getRef()), true);
				Allele aA = Allele.create(Character.toString(mkr.getAlt()), false);
				a.add(aR);
				a.add(aA);
				builderVc.alleles(a);
				builderVc.start(mkr.getPosition());
				builderVc.stop(mkr.getPosition());
				builderVc.id(mkr.getName());

				Collection<Genotype> genos = new ArrayList<Genotype>();
				byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilters, markerData.getMarkerName(), gcThreshold, proj.getLog());
				for (int k = 0; k < idsToInclude.size(); k++) {
					int idInd = idIndexMap.get(idsToInclude.get(k));
					String id = idsToInclude.get(k);
					Genotype g;
					List<Allele> all;
					byte indGeno = genotypes[idInd];
//					0 for A/A, 1 for A/B, 2 for B/B, and -1 for null
					boolean aIsRef = mkr.getRef() == mkr.getA();
					switch(indGeno) {
						case 0:
							all = aIsRef ? Arrays.asList(aR, aR) : Arrays.asList(aA, aA);
							break;
						case 1:
							all = aIsRef ? Arrays.asList(aR, aA) : Arrays.asList(aR, aA);
							break;
						case 2:
							all = aIsRef ? Arrays.asList(aA, aA) : Arrays.asList(aR, aR);
							break;
						case -1:
						default:
							all = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
							break;
					}
					g = GenotypeBuilder.create(id, all);
					
					genos.add(g);
				}
				
				builderVc.genotypes(genos);
				writer.add(builderVc.make());
				
			}
			mdl.shutdown();
			mdl = null;
			
			writer.close();
			writer = null;
			System.gc();
		}
		
		
		
	}
	
}
