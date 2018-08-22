package org.genvisis.common.parsing;

public class ByteWrapperColumn extends NumberWrapperColumn<Byte> {

  public ByteWrapperColumn(IndexedFileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public ByteWrapperColumn(IndexedFileColumn<?> base, boolean dieOnParseFailure) {
    super(base, Byte::parseByte, dieOnParseFailure);
  }

}
