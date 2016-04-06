package org.flowcyt.cfcs;
// CFCSListModeData.java

/* ------------------------------------------------------------------------- *\
This software and documentation are provided 'as is' and Tree Star, Inc., its
contractors and partners specifically disclaim all other warranties, expressed
or implied, including but not limited to implied warranties of merchantability
and fitness for a particular purpose, or during any particular date range.

By using this software, you are agreeing to these limits of liability, and to
hold Tree Star harmless for any information, accurate or erroneous, that might
be generated by the program.  This software is intended for research use only.

Christopher Lane <cdl@best.classes> for Tree Star  1/16/2002      Copyright 2002
\* ------------------------------------------------------------------------- */


import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import common.ext;

public final class CFCSListModeData extends CFCSAbstractData implements CFCSErrorCodes
{

    // --------------------------------------------------------------------

    /* friendly */
    CFCSListModeData(final CFCSDatatype datatype, final CFCSKeywords keywords)
    {
        super(LISTMODE, datatype, keywords);

        sizing = new CFCSDataSizing(keywords)
        {
            public int getSizeVariability()
            {
                return COLUMN;
            }

            public boolean isRangeMaskable()
            {
                return true;
            }
        };

        datatype.setSizingObject(sizing);
    }

    // --------------------------------------------------------------------

    public final byte[][][] cinchBytes(final byte[] bytes)
    {
        int[] sizes = null;
        int size = 0;
	    final int nEvents = getCount();
	    final int nParameters = parameters.getCount();

	    final boolean variable = (sizing.getSizeVariability() != CFCSDataSizing.FIXED);

        final byte[][][] cinched = new byte[nEvents][nParameters][];

        if (variable)
            sizes = sizing.getByteSizes();
        else
            size = sizing.getByteSize();

        for (int event = 0, position = 0; event < nEvents; event++)
        {
            for (int parameter = 0, used = 0; parameter < nParameters; parameter++, position += used)
            {
                if (variable)
                    size = sizes[parameter];

                if (size == CFCSParameter.FIELDSIZE_FREE_FORMAT)
                {
                    char character;
                    StringBuffer buffer = new StringBuffer();

                    used = 0;

                    // JS, Feb 24, 2006
                    try {
                    	while (isDelimiter(character = (char) bytes[position + used++]))
                    	{
                    	}
                    	do
                    	{
                    		buffer.append(character);
                    	}
                    	while (!isDelimiter(character = (char) bytes[position + used++]));
                    } catch(java.lang.ArrayIndexOutOfBoundsException e) {
                    }
                    
                    cinched[event][parameter] = (buffer.toString()).getBytes();
                }
                else
                {
                    byte[] cinch = cinched[event][parameter] = new byte[size];

                    System.arraycopy(bytes, position, cinch, 0, used = size);
                }
            }
        }

        return cinched;
    }

    // --------------------------------------------------------------------
    
    private volatile int loadCount = 0;
    public boolean isLoaded() {
        return loadCount == 2;
    }
    
    public final void setBytes(final byte[] bytes)
    {
        super.setBytes(bytes);
        
        ((ArrayList<Object>)(((CFCSAbstractDatatype)datatype).data)).ensureCapacity((int) (getCount() * 1.76));
        for (int i = 0; i < getCount(); i++) {
            ((ArrayList<Object>)(((CFCSAbstractDatatype)datatype).data)).add(new float[0]);
        }
        
        // this thread:
//        new Thread(new Runnable() {
//            @Override
//            public void run() {
//                for (int event = 0, events = getCount(); event < events; event++) {
//                    try {
//                        datatype.readData(event, cinched[event]);
//                    } catch (IndexOutOfBoundsException exception) {
//                        throw new CFCSError(CFCSNoSuchEvent, event);
//                    }
//                }
//                loadCount += 1;
//            }
//        }).start();
        
        new Thread(new Runnable() {
            @Override
            public void run() {
                for (int event = 0, events = getCount() / 2; event < events; event++)
                {
                    try
                    {
                        datatype.readData(event, cinched[event]);
                    }
                    catch (IndexOutOfBoundsException exception)
                    {
                        throw new CFCSError(CFCSNoSuchEvent, event);
                    }
                }
                loadCount += 1;
            }
        }).start();
        new Thread(new Runnable() {
            @Override
            public void run() {
                for (int event = getCount() / 2, events = getCount(); event < events; event++)
                {
                    try
                    {
                        datatype.readData(event, cinched[event]);
                    }
                    catch (IndexOutOfBoundsException exception)
                    {
                        throw new CFCSError(CFCSNoSuchEvent, event);
                    }
                }
                loadCount += 1;
            }
        }).start();
        

//        int step = getCount() / 3;
//        
//        new Thread(new Runnable() {
//            @Override
//            public void run() {
//                for (int event = 0, events = step; event < events; event++)
//                {
//                    try
//                    {
//                        datatype.readData(event, cinched[event]);
//                    }
//                    catch (IndexOutOfBoundsException exception)
//                    {
//                        throw new CFCSError(CFCSNoSuchEvent, event);
//                    }
//                }
//                loadCount += 1;
//            }
//        }).start();
//        new Thread(new Runnable() {
//            @Override
//            public void run() {
//                for (int event = step, events = step * 2; event < events; event++)
//                {
//                    try
//                    {
//                        datatype.readData(event, cinched[event]);
//                    }
//                    catch (IndexOutOfBoundsException exception)
//                    {
//                        throw new CFCSError(CFCSNoSuchEvent, event);
//                    }
//                }
//                loadCount += 1;
//            }
//        }).start();
//        new Thread(new Runnable() {
//            @Override
//            public void run() {
//                for (int event = step * 2, events = getCount(); event < events; event++)
//                {
//                    try
//                    {
//                        datatype.readData(event, cinched[event]);
//                    }
//                    catch (IndexOutOfBoundsException exception)
//                    {
//                        throw new CFCSError(CFCSNoSuchEvent, event);
//                    }
//                }
//                loadCount += 1;
//            }
//        }).start();
        
    }
    
    // --------------------------------------------------------------------

    public final byte[] getBytes()
    {
        final int count = parameters.getCount();

        final ByteArrayOutputStream buffer = new ByteArrayOutputStream();

        try
        {
            DataOutputStream stream = new DataOutputStream(buffer);

            for (int event = 0, events = getCount(); event < events; event++)
            {
                try
                {
                    datatype.writeData(event, count, stream);
                }
                catch (IndexOutOfBoundsException exception)
                {
                    throw new CFCSError(CFCSNoSuchEvent, event);
                }
            }

            stream.close();
        }
        catch (IOException exception)
        {
            throw new CFCSError(CFCSIOError, exception);
        }

        byte[] bytes = buffer.toByteArray();

        return (sizing.isPackedData()) ? packBytes(bytes) : bytes;
    }

    // --------------------------------------------------------------------
    // classespute how many bits to truncate to stuff large data into little
    // containers.  Based on the size of containers and range of parameters.

    private int classesputeShift(final int have, final int need)
    {
        int shift, range = Integer.MIN_VALUE;
        final int limit = (1 << (need * BITSPERBYTE));

        for (int i = 0, count = parameters.getCount(); i < count; i++)
        {
            range = Math.max((parameters.getParameter(i)).getRange(), range);
        }

        // Data may have been adjusted once already to fit container
        range = Math.min(range, (1 << (have * BITSPERBYTE)));

        for (shift = 0; range > limit; range >>= 2, shift++)
        {
        }

        return shift;
    }

    // --------------------------------------------------------------------

    protected final byte[] packBytes(final byte[] bytes)
    {
        final int[] sizes = sizing.getSizes();
        final int[] byteSizes = sizing.getByteSizes();
        final int[] masks = sizing.getFieldMasks();

        final ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        int current = 0;
	    int valid = 0;
	    final int count = parameters.getCount();
	    long accumulator = 0L;

        try
        {
            DataOutputStream stream = new DataOutputStream(buffer);

            for (int event = 0, events = getCount(); event < events; event++)
            {
                for (int parameter = 0; parameter < count; parameter++)
                {
                    int datum = 0, size = sizes[parameter];

                    for (int i = 0, nBytes = byteSizes[parameter]; i < nBytes; i++)
                    {
                        datum = (datum << BITSPERBYTE) | (bytes[current++] & 0xff);
                    }

                    accumulator = (accumulator << size) | (datum & masks[parameter]);

                    valid += size;

                    while (valid > BITSPERBYTE)
                    {
                        valid -= BITSPERBYTE;
                        stream.write((int) ((accumulator >>> valid) & 0xff));
                        accumulator &= (1 << valid) - 1;
                    }
                }
            }

            if (valid > 0)
                stream.write((int) (accumulator << (BITSPERBYTE - valid)));

            stream.close();
        }
        catch (IOException exception)
        {
            throw new CFCSError(CFCSIOError, exception);
        }

        return buffer.toByteArray();
    }

    // --------------------------------------------------------------------

    public final byte[] unpackBytes(final byte[] bytes)
    {
        final int[] swap = keywords.getByteSwapArray();

        boolean swapped = false;

        for (int i = 0; i < swap.length; i++)
        {
            if ((swapped = (swap[i] != JAVA_BYTESWAP_ARRAY[i])) == true)
                break;
        }

        if (swapped)
            return unpackSwappedBytes(bytes);

        return unpackNativeBytes(bytes);
    }

    // --------------------------------------------------------------------

    protected static final byte[] unpackSwappedBytes(final byte[] bytes)
    {
        throw new CFCSError(CFCSNotImplemented);
    }

    // --------------------------------------------------------------------

    protected final byte[] unpackNativeBytes(final byte[] bytes)
    {
        final int[] sizes = sizing.getSizes();
        final int[] byteSizes = sizing.getByteSizes();
        final int[] masks = sizing.getFieldMasks();

        final ByteArrayOutputStream buffer = new ByteArrayOutputStream();
        int current = 0;
	    int valid = 0;
	    int accumulator = 0;
	    final int count = parameters.getCount();

	    try
        {
            DataOutputStream stream = new DataOutputStream(buffer);

            for (int event = 0, events = getCount(); event < events; event++)
            {
                for (int parameter = 0; parameter < count; parameter++)
                {
                    int size = sizes[parameter];

                    for (; valid < size; valid += BITSPERBYTE)
                    {
                        accumulator = (accumulator << BITSPERBYTE) | (bytes[current++] & 0xff);
                    }

                    valid -= size;

                    int datum = (accumulator >>> valid) & masks[parameter];

                    accumulator &= (1 << valid) - 1;

                    switch (byteSizes[parameter])
                    {
                    case CFCSDatatypeBinary.BYTE:
                        stream.write(datum);
                        break;
                    case CFCSDatatypeBinary.SHORT:
                        stream.writeShort(datum);
                        break;
                    case CFCSDatatypeBinary.INTEGER:
                        stream.writeInt(datum);
                        break;
                    default :
                        throw new CFCSError(CFCSNotImplemented);
                    }
                }
            }

            stream.close();
        }
        catch (IOException exception)
        {
            throw new CFCSError(CFCSIOError, exception);
        }

        return buffer.toByteArray();
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int index, final byte[] array)
    {
        try
        {
            datatype.getData(index, array);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
        catch (CFCSDataSizeError exception)
        {
            short[] event = new short[array.length];
            int shift = classesputeShift(CFCSDatatype.SHORT, CFCSDatatype.BYTE);

            getEvent(index, event);

            for (int i = 0; i < event.length; i++)
                array[i] = (byte) (event[i] >>> shift);
        }
    }

    public final void setEvent(final int index, final byte[] array)
    {
        try
        {
            datatype.setData(index, array);
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            setEvent(index, event);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }
    }

    public final void addEvent(final byte[] array)
    {
        try
        {
            datatype.addData(array);
            addCount(1); // Moved from below
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            addEvent(event);
        }

        // JS: Moved into the try section, otherwise the count gets added multiple times if exception occurs
        // addCount(1);
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int index, final short[] array)
    {
        try
        {
            datatype.getData(index, array);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
        catch (CFCSDataSizeError exception)
        {
            int[] event = new int[array.length];
            int shift = classesputeShift(CFCSDatatype.INTEGER, CFCSDatatype.SHORT);

            getEvent(index, event);

            for (int i = 0; i < event.length; i++)
                array[i] = (short) (event[i] >>> shift);
        }
    }

    public final void setEvent(final int index, final short[] array)
    {
        try
        {
            datatype.setData(index, array);
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            setEvent(index, event);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }
    }

    public final void addEvent(final short[] array)
    {
        try
        {
            datatype.addData(array);
            addCount(1); // Moved from below
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            addEvent(event);
        }

        // JS: Moved into the try section, otherwise the count gets added multiple times if exception occurs
        // addCount(1);
        
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int index, final int[] array)
    {
        try
        {
            datatype.getData(index, array);
            for (int i = 0; i < array.length; i++)
            {
                array[i] = (int) (parameters.getParameter(i)).channelToScale(array[i]);
            }
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            getEvent(index, event);

            for (int i = 0; i < event.length; i++)
            {
                array[i] = (parameters.getParameter(i)).scaleToChannel(event[i]);
            }
        }
        catch (CFCSDataSizeError exception)
        {
            throw new CFCSError(CFCSNotImplemented);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
    }

    public final void setEvent(final int index, final int[] array)
    {
        try
        {
            datatype.setData(index, array);
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            setEvent(index, event);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }
    }

    public final void addEvent(final int[] array)
    {
        try
        {
            datatype.addData(array);
            addCount(1); // Moved from below
        }
        catch (CFCSDataScaleError exception)
        {
            double event[] = new double[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).channelToScale(array[i]);
            }

            addEvent(event);
        }

        // JS: Moved into the try section, otherwise the count gets added multiple times if exception occurs
        // addCount(1);
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int index, final float[] array)
    {
        try
        {
            datatype.getData(index, array);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            try {
				datatype.getData(index, event);
			} catch (CFCSDataSizeError e) {
				throw new CFCSError(CFCSNotImplemented);
			}

            for (int i = 0; i < event.length; i++)
            {
                array[i] = (float) (parameters.getParameter(i)).channelToScale(event[i]);
            }
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
    }

    public final void setEvent(final int index, final float[] array)
    {
        try
        {
            datatype.setData(index, array);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).scaleToChannel(array[i]);
            }

            setEvent(index, event);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }
    }

    public final void addEvent(final float[] array)
    {
        try
        {
            datatype.addData(array);
            addCount(1); // Moved from below
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).scaleToChannel(array[i]);
            }

            addEvent(event);
        }

        // JS: Moved into the try section, otherwise the count gets added multiple times if exception occurs
        // addCount(1);
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int index, final double[] array)
    {
        try
        {
            datatype.getData(index, array);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            try {
				datatype.getData(index, event);
			} catch (CFCSDataSizeError e) {
				throw new CFCSError(CFCSNotImplemented);
			}

            for (int i = 0; i < event.length; i++)
            {
                array[i] = (parameters.getParameter(i)).channelToScale(event[i]);
            }
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
    }

    public final void getEventAsInTheFile(final int index, final double[] array)
    {
    	try
        {
            datatype.getData(index, array);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];
            try {
				datatype.getData(index, event);
			} catch (CFCSDataSizeError e) {
				throw new CFCSError(CFCSNotImplemented);
			}

            for (int i = 0; i < event.length; i++)
            {
                array[i] = (double)event[i];
            }
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchEvent, index);
        }
    }
    
    public final void setEvent(final int index, final double[] array)
    {
        try
        {
            datatype.setData(index, array);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).scaleToChannel(array[i]);
            }

            setEvent(index, event);
        }
        catch (IndexOutOfBoundsException exception)
        {
            throw new CFCSError(CFCSNoSuchParameter, index);
        }
    }

    public final void addEvent(final double[] array)
    {
        try
        {
            datatype.addData(array);
            addCount(1);
        }
        catch (CFCSDataChannelError exception)
        {
            int event[] = new int[array.length];

            for (int i = 0; i < event.length; i++)
            {
                event[i] = (parameters.getParameter(i)).scaleToChannel(array[i]);
            }

            addEvent(event);
        }

        // JS: Moved into the try section, otherwise the count gets added multiple times if exception occurs
        // addCount(1);
    }

    // --------------------------------------------------------------------

    public final void getEvents(final int startIndex, final int endIndex, final int[][] table)
    {
        for (int index = startIndex; index <= endIndex; index++)
        {
            getEvent(index, table[index]);
        }
    }

    public final void setEvents(final int startIndex, final int endIndex, final int[][] table)
    {
        for (int index = startIndex; index <= endIndex; index++)
        {
            setEvent(index, table[index]);
        }
    }

    // --------------------------------------------------------------------

    public final void getEvent(final int event, final int nBytes, final byte[][] table)
    {
        final int count = parameters.getCount();

        if (table.length < count)
        {
            throw new CFCSError(CFCSInsufficientSpace, table.length);
        }

        for (int parameter = 0; parameter < count; parameter++)
        {
            byte[] bytes = table[parameter];

            if (bytes.length < nBytes)
            {
                throw new CFCSError(CFCSInsufficientSpace, bytes.length);
            }

            byte[] cinch = cinched[event][parameter];

            int length = cinch.length;
            int limit = nBytes - length;

            for (int i = 0; i < limit; i++)
                bytes[i] = 0;

            for (int i = 0; i < length; i++)
                bytes[limit + i] = cinch[i];
        }
    }

    // --------------------------------------------------------------------

    public static final void setEvent(final int index, final int nBytes, final byte[][] table)
    {
        throw new CFCSError(CFCSNotImplemented);
    }

    // --------------------------------------------------------------------

}
