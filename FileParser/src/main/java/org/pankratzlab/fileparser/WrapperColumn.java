package org.pankratzlab.fileparser;

import java.util.function.Function;

/**
 * A {@link CachedFileColumn} that wraps a base {@link FileColumn} of one type and determines a
 * calculated value in another type
 *
 * @param <W> The wrapping type, to be returned by {@link #getValue(String[])}
 * @param <B> The base type, that of the supplied base {@link FileColumn}
 */
public abstract class WrapperColumn<W, B> extends CachedFileColumn<W>
                                   implements IndexedFileColumn<W> {

  private final IndexedFileColumn<? extends B> base;

  /**
   * @param base
   */
  public WrapperColumn(IndexedFileColumn<? extends B> base) {
    this(base, base.dieOnParseFailure());
  }

  /**
   * @param base
   * @param dieOnParseFailure
   */
  public WrapperColumn(IndexedFileColumn<? extends B> base, boolean dieOnParseFailure) {
    this(base.getName(), base, dieOnParseFailure);
  }

  protected WrapperColumn(String name, IndexedFileColumn<? extends B> base,
                          boolean dieOnParseFailure) {
    super(name, dieOnParseFailure);
    this.base = base;
  }

  @Override
  public final void initialize(FileParser parser) {
    base.initialize(parser);
  }

  @Override
  public final W calculateValue(String[] line) throws ParseFailureException {
    return calculateFromBase(getBaseValue(line));
  }

  @Override
  public int getIndex() {
    return base.getIndex();
  }

  private final B getBaseValue(String[] line) throws ParseFailureException {
    return base.getValue(line);
  }

  /**
   * @param base the result of {@link FileColumn#getValue(String[])} for the base
   * @return value to return for {@link #getValue(String[])}
   */
  protected abstract W calculateFromBase(B base) throws ParseFailureException;

  public static <W, B> WrapperColumn<W, B> of(IndexedFileColumn<? extends B> base,
                                              final Function<B, W> calculateFunction) {
    return new WrapperColumn<W, B>(base) {

      @Override
      protected W calculateFromBase(B base) {
        return calculateFunction.apply(base);
      }
    };
  }

}
