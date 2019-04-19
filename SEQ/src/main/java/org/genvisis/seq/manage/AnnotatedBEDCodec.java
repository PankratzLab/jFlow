package org.genvisis.seq.manage;

import java.util.regex.Pattern;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;

public class AnnotatedBEDCodec extends AsciiFeatureCodec<BEDFeature> {

  /** Default extension for BED files. */
  public static final String BED_EXTENSION = ".bed";

  private static final Pattern SPLIT_PATTERN = Pattern.compile("\\t|( +)");
  private final int startOffsetValue;

  /**
   * Calls {@link #BEDCodec(StartOffset)} with an argument of {@code StartOffset.ONE}
   */
  public AnnotatedBEDCodec() {
    this(StartOffset.ONE);
  }

  /**
   * BED format is 0-based, but Tribble is 1-based. Set desired start position at either ZERO or ONE
   */
  public AnnotatedBEDCodec(final StartOffset startOffset) {
    super(BEDFeature.class);
    this.startOffsetValue = startOffset.value();
  }

  public BEDFeature decodeLoc(String line) {
    return decode(line);
  }

  @Override
  public BEDFeature decode(String line) {

    if (line.trim().isEmpty()) {
      return null;
    }
    // discard header lines in case the caller hasn't called readHeader
    if (isBEDHeaderLine(line)) {
      this.readHeaderLine(line);
      return null;
    }

    String[] tokens = SPLIT_PATTERN.split(line, -1);
    return decode(tokens);
  }

  /**
   * The BED codec doesn't retain the actual header, but we need to parse through it and advance to
   * the beginning of the first feature. This is especially true if we're indexing, since we want to
   * underlying stream offset to be established correctly, but is also the case for when we're
   * simply iterating to satisfy a feature query (otherwise the feature reader can terminate
   * prematurely if the header is large).
   * 
   * @param lineIterator
   * @return Always null. The BEDCodec currently doesn't model or preserve the BED header.
   */
  @Override
  public Object readActualHeader(final LineIterator lineIterator) {
    while (lineIterator.hasNext()) {
      // Only peek, since we don't want to actually consume a line of input unless its a header
      // line.
      // This prevents us from advancing past the first feature.
      final String nextLine = lineIterator.peek();
      if (isBEDHeaderLine(nextLine)) {
        // advance the iterator and consume the line (which is a no-op)
        this.readHeaderLine(lineIterator.next());
      } else {
        return null; // break out when we've seen the end of the header
      }
    }

    return null;
  }

  // Return true if the candidateLine looks like a BED header line.
  private boolean isBEDHeaderLine(final String candidateLine) {
    return candidateLine.startsWith("#") || candidateLine.startsWith("track")
           || candidateLine.startsWith("browser");
  }

  public BEDFeature decode(String[] tokens) {
    int tokenCount = tokens.length;

    // The first 3 columns are non optional for BED. We will relax this
    // and only require 2.

    if (tokenCount < 2) {
      return null;
    }

    String chr = tokens[0];

    // The BED format uses a first-base-is-zero convention, Tribble features use 1 => add 1.
    int start = Integer.parseInt(tokens[1]) + startOffsetValue;

    int end = start;
    if (tokenCount > 2) {
      end = Integer.parseInt(tokens[2]);
    }

    AnnotatedBEDFeature feature = new AnnotatedBEDFeature(chr, start, end);

    // Name
    if (tokenCount > 3) {
      String name = tokens[3].replaceAll("\"", "");
      feature.setName(name);
    }

    // The rest of the columns are values unrelated to the BED format, but associated with this
    // entry in the BED file

    for (int i = 4; i < tokens.length; i++) {
      feature.addAnnotation(tokens[i]);
    }

    return feature;
  }

  protected boolean readHeaderLine(String line) {
    // We don't parse BED header
    return false;
  }

  @Override
  public boolean canDecode(final String path) {
    final String toDecode;
    if (AbstractFeatureReader.hasBlockCompressedExtension(path)) {
      toDecode = path.substring(0, path.lastIndexOf("."));
    } else {
      toDecode = path;
    }
    return toDecode.toLowerCase().endsWith(BED_EXTENSION);
  }

  public int getStartOffset() {
    return this.startOffsetValue;
  }

  /**
   * Indicate whether co-ordinates or 0-based or 1-based.
   * <p/>
   * Tribble uses 1-based, BED files use 0. e.g.: start_position = bedline_start_position -
   * startIndex.value()
   */
  public enum StartOffset {
    ZERO(0), ONE(1);

    private int start;

    private StartOffset(int start) {
      this.start = start;
    }

    public int value() {
      return this.start;
    }
  }

  @Override
  public TabixFormat getTabixFormat() {
    return TabixFormat.BED;
  }
}
