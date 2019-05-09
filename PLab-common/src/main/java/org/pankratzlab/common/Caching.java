package org.pankratzlab.common;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.lang.ref.SoftReference;
import java.util.function.Supplier;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;

public class Caching {

  private Caching() {}

  /**
   * Returns a supplier which caches in a {@link SoftReference} the instance retrieved during the
   * first call to {@code get()} and returns that value on subsequent calls to {@code get()} until
   * the referent is Garbage Collected. See:
   * <a href="http://en.wikipedia.org/wiki/Memoization">memoization</a>
   *
   * <p>
   * The returned supplier is thread-safe and will only call the delegate once unless the referent
   * has been Garbage Collected or manually cleared. The supplier's serialized form does not contain
   * the cached value, which will be recalculated when {@code get()} is called on the reserialized
   * instance.
   *
   * <p>
   * When the underlying delegate throws an exception then this memoizing supplier will keep
   * delegating calls until it returns valid data.
   *
   * <p>
   * If {@code delegate} is an instance created by an earlier call to
   * {@code memoizeWithSoftReference}, it is returned directly.
   */
  public static <T> SoftRefMemoizingSupplier<T> memoizeWithSoftRef(Supplier<T> delegate) {
    if (delegate instanceof NonSerializableSoftRefMemoizingSupplier
        || delegate instanceof SerializableSoftRefMemoizingSupplier) {
      return (SoftRefMemoizingSupplier<T>) delegate;
    }
    return delegate instanceof Serializable ? new SerializableSoftRefMemoizingSupplier<>(delegate)
                                            : new NonSerializableSoftRefMemoizingSupplier<>(delegate);
  }

  public static interface SoftRefMemoizingSupplier<T> extends Supplier<T> {
    public void clear();
  }

  private static class SerializableSoftRefMemoizingSupplier<T>
                                                           implements SoftRefMemoizingSupplier<T>,
                                                           Serializable {
    final Supplier<T> delegate;
    transient volatile SoftReference<T> value = new SoftReference<>(null);

    SerializableSoftRefMemoizingSupplier(Supplier<T> delegate) {
      this.delegate = Preconditions.checkNotNull(delegate);
    }

    @Override
    public T get() {
      T cachedVal = value.get();
      if (cachedVal == null) {
        synchronized (this) {
          cachedVal = value.get();
          if (cachedVal == null) {
            T t = delegate.get();
            value = new SoftReference<>(t);
            return t;
          }
        }
      }
      return cachedVal;
    }

    @Override
    public String toString() {
      return "Caching.memoizeWithSoftRef(" + delegate + ")";
    }

    private static final long serialVersionUID = 0;

    private void readObject(ObjectInputStream inputStream) throws IOException,
                                                           ClassNotFoundException {
      inputStream.defaultReadObject();
      value = new SoftReference<>(null);
    }

    @Override
    public void clear() {
      synchronized (this) {
        value = new SoftReference<>(null);
      }
    }
  }

  @VisibleForTesting
  static class NonSerializableSoftRefMemoizingSupplier<T> implements SoftRefMemoizingSupplier<T> {
    final Supplier<T> delegate;
    volatile SoftReference<T> value = new SoftReference<>(null);

    NonSerializableSoftRefMemoizingSupplier(Supplier<T> delegate) {
      this.delegate = Preconditions.checkNotNull(delegate);
    }

    @Override
    public T get() {
      T cachedVal = value.get();
      if (cachedVal == null) {
        synchronized (this) {
          cachedVal = value.get();
          if (cachedVal == null) {
            T t = delegate.get();
            value = new SoftReference<>(t);
            return t;
          }
        }
      }
      return cachedVal;
    }

    @Override
    public String toString() {
      return "Caching.memoizeWithSoftRef(" + delegate + ")";
    }

    @Override
    public void clear() {
      synchronized (this) {
        value = new SoftReference<>(null);
      }
    }
  }

}
