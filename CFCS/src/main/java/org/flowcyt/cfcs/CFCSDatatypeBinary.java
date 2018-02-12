package org.flowcyt.cfcs;
// CFCSDatatypeBinary.java

/*
 * ------------------------------------------------------------------------- *\ This software and
 * documentation are provided 'as is' and Tree Star, Inc., its contractors and partners specifically
 * disclaim all other warranties, expressed or implied, including but not limited to implied
 * warranties of merchantability and fitness for a particular purpose, or during any particular date
 * range. By using this software, you are agreeing to these limits of liability, and to hold Tree
 * Star harmless for any information, accurate or erroneous, that might be generated by the program.
 * This software is intended for research use only. Christopher Lane <cdl@best.classes> for Tree
 * Star 1/14/2002 Copyright 2002 \*
 * -------------------------------------------------------------------------
 */

import java.io.DataOutputStream;
import java.io.IOException;

public class CFCSDatatypeBinary extends CFCSAbstractDatatype {

  private final static int[] SUBDATATYPES = {UNDEFINED, BYTE, SHORT, UNDEFINED, INTEGER};

  // --------------------------------------------------------------------

  /* friendly */
  static CFCSDatatype getDataSubtype(final CFCSParameters parameters) {
    final int count = parameters.getCount();

    if (count == 0) return new CFCSDatatypeBinary(); // For creation, use most general type

    int max = Integer.MIN_VALUE, min = Integer.MAX_VALUE;

    for (int i = 0; i < count; i++) {
      CFCSParameter parameter = parameters.getParameter(i);

      int bytes = (parameter.getFieldSize() + 7) / 8;

      if (bytes > max) max = bytes;
      if (bytes < min) min = bytes;
    }

    if (min == max && min < SUBDATATYPES.length) {
      switch (SUBDATATYPES[min]) { // standard fixed size type
        case CFCSDatatypeBinary.BYTE:
          return new CFCSDatatypeByte();
        case CFCSDatatypeBinary.SHORT:
          return new CFCSDatatypeShort();
        case CFCSDatatypeBinary.INTEGER:
          return new CFCSDatatypeInteger();
        default:
          throw new CFCSError(CFCSSystemError);
      }
    } else if (min < SUBDATATYPES.length && max < SUBDATATYPES.length) {
      return new CFCSDatatypeBinary();
    }

    return new CFCSDatatypeOversize();
  }

  // --------------------------------------------------------------------

  public void readData(final int index, final byte[][] cinch) {
    final int[] row = new int[cinch.length];

    final boolean maskable = sizing.isRangeMaskable();
    final int[] masks = (maskable) ? sizing.getRangeMasks() : null;

    for (int cell = 0; cell < cinch.length; cell++) {
      int bytes = cinch[cell].length;

      if (bytes > INTEGER) throw new CFCSError(CFCSNotImplemented); // fix

      int datum = (cinch[cell][0] & 0xff);

      for (int i = 1; i < bytes; i++) {
        datum = (datum << 8) | (cinch[cell][i] & 0xff);
      }

      if (maskable) datum &= masks[cell];

      row[cell] = datum;
    }

    if (index == data.size()) data.add(row);
    else data.set(index, row);
  }

  // --------------------------------------------------------------------

  public void writeData(final int index, final int count, final DataOutputStream stream) {
    final int[] row = (int[]) data.get(index);

    int bytes = 0;
    final int variability = sizing.getSizeVariability();
    final int[] sizes = (variability != CFCSDataSizing.FIXED) ? sizing.getByteSizes() : null;

    if (variability == CFCSDataSizing.ROW) bytes = sizes[index];
    else if (variability == CFCSDataSizing.FIXED) bytes = sizing.getByteSize();

    try {
      for (int i = 0; i < count; i++) {
        int datum = row[i]; // Sign extension should be OK here

        if (variability == CFCSDataSizing.COLUMN) bytes = sizes[i];

        switch (bytes) {
          case BYTE:
            stream.write(datum);
            break;
          case SHORT:
            stream.writeShort(datum);
            break;
          case INTEGER:
            stream.writeInt(datum);
            break;
          default:
            throw new CFCSError(CFCSNotImplemented);
        }
      }
    } catch (IOException exception) {
      throw new CFCSError(CFCSIOError, exception);
    }
  }

  // --------------------------------------------------------------------
  // Variations on getData(index, array) for byte, short, int, float & double

  public void getData(final int index, final byte[] array) throws CFCSDataSizeError {
    final int[] row = (int[]) data.get(index);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientSpace, array.length);
    }

    for (int i = 0; i < row.length; i++) {
      int datum = row[i];

      if (datum > Byte.MAX_VALUE || datum < Byte.MIN_VALUE) {
        throw new CFCSDataSizeError();
      }

      array[i] = (byte) datum;
    }
  }

  public void getData(final int index, final short[] array) throws CFCSDataSizeError {
    final int[] row = (int[]) data.get(index);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientSpace, array.length);
    }

    for (int i = 0; i < row.length; i++) {
      int datum = row[i];

      if (datum > Short.MAX_VALUE || datum < Short.MIN_VALUE) {
        throw new CFCSDataSizeError();
      }

      array[i] = (short) datum;
    }
  }

  public void getData(final int index, final int[] array) throws CFCSDataSizeError {
    final int[] row = (int[]) data.get(index);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientSpace, array.length);
    }

    System.arraycopy(row, 0, array, 0, row.length);
  }

  public final void getData(final int index, final float[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  public final void getData(final int index, final double[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  // --------------------------------------------------------------------

  public void setData(final int index, final byte[] array) {
    final int[] row = (int[]) data.get(index);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientData, array.length);
    }

    for (int i = 0; i < row.length; i++)
      row[i] = (array[i] & 0xff);
  }

  public void setData(final int index, final short[] array) {
    final int[] row = (int[]) data.get(index);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientData, array.length);
    }

    for (int i = 0; i < row.length; i++)
      row[i] = (array[i] & 0xffff);
  }

  public void setData(final int index, final int[] array) {
    final int[] row = (int[]) data.get(index);

    System.arraycopy(array, 0, row, 0, row.length);

    if (array.length < row.length) {
      throw new CFCSError(CFCSInsufficientData, array.length);
    }
  }

  public final void setData(final int index, final float[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  public final void setData(final int index, final double[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  // --------------------------------------------------------------------

  public void addData(final byte[] array) {
    final int[] row = new int[array.length];

    for (int i = 0; i < row.length; i++)
      row[i] = (array[i] & 0xff);

    data.add(row);
  }

  public void addData(final short[] array) {
    final int[] row = new int[array.length];

    for (int i = 0; i < row.length; i++)
      row[i] = (array[i] & 0xff);

    data.add(row);
  }

  public void addData(final int[] array) {
    final int[] row = new int[array.length];

    System.arraycopy(array, 0, row, 0, row.length);

    data.add(row);
  }

  public final void addData(final float[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  public final void addData(final double[] array) throws CFCSDataChannelError {
    throw new CFCSDataChannelError();
  }

  // --------------------------------------------------------------------

}
