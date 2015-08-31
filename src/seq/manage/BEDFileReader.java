package seq.manage;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212<br>
 *         Trying to make an indexed bed file reader, mimicking {@link VCFFileReader}
 */
public class BEDFileReader implements Closeable, Iterable<BEDFeature> {
	private final FeatureReader<BEDFeature> reader;

	public BEDFileReader(final String file, final boolean requireIndex) {
		this.reader = AbstractFeatureReader.getFeatureReader(file, new BEDCodec(), requireIndex);
	}

	/** Queries for records within the region specified. */
	public CloseableIterator<BEDFeature> query(final String chrom, final int start, final int end) {
		try {
			return reader.query(chrom, start, end);
		} catch (final IOException ioe) {
			throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
		}
	}

	public void close() {
		try {
			this.reader.close();
		} catch (final IOException ioe) {
			throw new TribbleException("Could not close a bed context feature reader.", ioe);
		}
	}

	/** Returns an iterator over all records in this VCF/BCF file. */
	public CloseableIterator<BEDFeature> iterator() {
		try {
			return reader.iterator();
		} catch (final IOException ioe) {
			throw new TribbleException("Could not create an iterator from a feature reader.", ioe);
		}
	}

}
