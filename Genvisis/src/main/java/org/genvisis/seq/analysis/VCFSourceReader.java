package org.genvisis.seq.analysis;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * This is a complete recapitulation of the {@link VCFFileReader} to allow creation with non-file
 * sources, as they are not strictly required elsewhere in the API.
 * <p>
 * TODO Ideally to be submitted upstream, or at least a common interface could be provided and used
 * in Genvisis API (as this class must itself be a VCFFileReader).
 * </p>
 */
public class VCFSourceReader extends VCFFileReader implements Closeable, Iterable<VariantContext> {

	private final FeatureReader<VariantContext> reader;

	public VCFSourceReader(final String file, final boolean requireIndex) {
		super(new File(file), requireIndex);// unfortunately, I think passing null to this constructor
																				// will fail...so maybe we should just use the
																				// actual constructors in most places

		reader = AbstractFeatureReader.getFeatureReader(file,
																										file.endsWith(".bcf")	? (FeatureCodec) new BCF2Codec()
																																					: new VCFCodec(),
																										requireIndex);
	}

	/** Returns the VCFHeader associated with this VCF/BCF file. */
	@Override
	public VCFHeader getFileHeader() {
		return (VCFHeader) reader.getHeader();
	}

	/** Returns an iterator over all records in this VCF/BCF file. */
	@Override
	public CloseableIterator<VariantContext> iterator() {
		try {
			return reader.iterator();
		} catch (final IOException ioe) {
			throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
		}
	}

	/** Queries for records within the region specified. */
	@Override
	public CloseableIterator<VariantContext> query(	final String chrom, final int start,
																									final int end) {
		try {
			return reader.query(chrom, start, end);
		} catch (final IOException ioe) {
			throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
		}
	}

	@Override
	public void close() {
		try {
			reader.close();
		} catch (final IOException ioe) {
			throw new TribbleException("Could not close a variant context feature reader.", ioe);
		}
	}
}
