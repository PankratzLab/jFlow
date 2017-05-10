package org.genvisis.one.JL.aging.bam;

import java.util.List;

import org.genvisis.common.Logger;

import htsjdk.samtools.SAMRecord;

// mtDNA CN, telomere length, chr X/Y loss will implement this interface ... goal is to read bam file once
/**
 * For use in analyzing a .bam file, focusing on a single read at a time
 *
 */
public interface BamAnalysis {

	/**
	 * @return base file name pattern for reporting
	 */
	public String getRootOutputFile();

	/**
	 * Prepare any primary information from the bam file
	 * 
	 * @param bamFile
	 * @param log
	 */
	public void init(String bamFile, Logger log);

	/**
	 * @param samRecord
	 *            analyze this read
	 * @param log
	 */
	public void analyze(SAMRecord samRecord, Logger log);

	/**
	 * @param samRecord
	 *            do any filtering on this record here
	 * @param log
	 * @return whether to use the read or not
	 */
	public boolean shouldAnalyze(SAMRecord samRecord, Logger log);

	/**
	 * Summarize the run, called at end of processing
	 */
	public void summarize();

	/**
	 * @return the summary output to be written/passed on
	 */
	public List<String> getOutput();

}
