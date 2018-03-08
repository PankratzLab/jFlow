package org.genvisis.common.parsing;

public class ByteWrapperColumn extends NumberWrapperColumn<Byte> {

  public ByteWrapperColumn(FileColumn<?> base) {
    this(base, base.dieOnParseFailure());
  }

  public ByteWrapperColumn(FileColumn<?> base, boolean dieOnParseFailure) {
    super(base, Byte::parseByte, dieOnParseFailure);
  }

}
