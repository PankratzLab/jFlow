package org.pankratzlab.common;

import java.io.Serializable;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;

public class Functionals {

  private Functionals() {}

  public static interface SerializableSupplier<T> extends Supplier<T>, Serializable {};

  public static <T> SerializableSupplier<T> serializableSupplier(final SerializableSupplier<T> supplier) {
    return supplier;
  }

  public static interface SerializableConsumer<T> extends Consumer<T>, Serializable {};

  public static <T> SerializableConsumer<T> serializableConsumer(final SerializableConsumer<T> consumer) {
    return consumer;
  }

  public static interface SerializableFunction<F, T> extends Function<F, T>, Serializable {};

  public static <F, T> SerializableFunction<F, T> serializableFunction(final SerializableFunction<F, T> function) {
    return function;
  }

}
