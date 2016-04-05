package org.flowcyt.cfcs;
// CFCSKeywords.java

/* ------------------------------------------------------------------------- *\
This software and documentation are provided 'as is' and Tree Star, Inc., its
contractors and partners specifically disclaim all other warranties, expressed
or implied, including but not limited to implied warranties of merchantability
and fitness for a particular purpose, or during any particular date range.

By using this software, you are agreeing to these limits of liability, and to
hold Tree Star harmless for any information, accurate or erroneous, that might
be generated by the program.  This software is intended for research use only.

Christopher Lane <cdl@best.classes> for Tree Star  1/18/2002      Copyright 2002
\* ------------------------------------------------------------------------- */


import java.io.UnsupportedEncodingException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public final class CFCSKeywords implements CFCSErrorCodes
{

    /* friendly */
    static final String DATATYPE_KEYWORD = "$DATATYPE";
    /* friendly */
    static final String GATE_KEYWORD = "$GATE";
    /* friendly */
    static final String MODE_KEYWORD = "$MODE";
    /* friendly */
    static final String PARAMETER_KEYWORD = "$PAR";

    /* friendly */
    static final String BYTEORDER_KEYWORD = "$BYTEORD";
    /* friendly */
    static final String FILE_KEYWORD = "$FIL";
    /* friendly */
    static final String NEXTDATA_KEYWORD = "$NEXTDATA";
    /* friendly */
    static final String TOTAL_KEYWORD = "$TOT";
    
    /* friendly */
    // JS, added for FCS3.1 support
    // Eventually, we could add support for "all" FCS keywords
    static final String SAMPLE_VOLUME_KEYWORD = "$VOL";
    static final String PLATEID_KEYWORD = "$PLATEID";
    static final String PLATENAME_KEYWORD = "$PLATENAME";
    static final String WELLID_KEYWORD = "$WELLID";
    static final String SPILLOVER_KEYWORD = "$SPILLOVER";
    static final String SPILLOVER_ALTERNATIVE_KEYWORD = "SPILL"; // This is not in the standard but BD uses it in their FCS3.0 files
    static final String ORIGINALITY_KEYWORD = "$ORIGINALLITY";
    static final String LAST_MODIFIER_KEYWORD = "$LAST_MODIFIER";
    static final String LAST_MODIFIED_KEYWORD = "$LAST_MODIFIED";
    
    public static enum OriginalityEnum {Original, NonDataModified, Appended, DataModified, Undefined };
        
    /* friendly */
    static final Set<String> SYSTEM_KEYWORDS = new HashSet<String>();
    /* friendly */
    static final Set<String> NONTRANSFERABLE_KEYWORDS = new HashSet<String>();

    /* friendly */
    static final String[] SEGMENT_ROOTS = {"HEADER", "STEXT", "DATA", "ANALYSIS"};
    /* friendly */
    static final String SEGMENT_BEGIN_PREFIX = "$BEGIN";
    /* friendly */
    static final String SEGMENT_END_PREFIX = "$END";

    static
    {
        SYSTEM_KEYWORDS.add(DATATYPE_KEYWORD);
        SYSTEM_KEYWORDS.add(GATE_KEYWORD);
        SYSTEM_KEYWORDS.add(MODE_KEYWORD);
        SYSTEM_KEYWORDS.add(PARAMETER_KEYWORD);

        NONTRANSFERABLE_KEYWORDS.add(BYTEORDER_KEYWORD);
        NONTRANSFERABLE_KEYWORDS.add(FILE_KEYWORD);
        NONTRANSFERABLE_KEYWORDS.add(NEXTDATA_KEYWORD);
        NONTRANSFERABLE_KEYWORDS.add(TOTAL_KEYWORD);

        for (int i = CFCSDataSet.TEXT; i < CFCSDataSet.OTHER_START; i++)
        {
            NONTRANSFERABLE_KEYWORDS.add(SEGMENT_BEGIN_PREFIX + SEGMENT_ROOTS[i]);
            NONTRANSFERABLE_KEYWORDS.add(SEGMENT_END_PREFIX + SEGMENT_ROOTS[i]);
        }

        SYSTEM_KEYWORDS.addAll(NONTRANSFERABLE_KEYWORDS);
    }

    private static final char DEFAULT_DELIMITER = '/';
    private static final char DELIMITER_GUARD = '_'; // can't be the same as delimiter
    private static final char[] delimiters = {
        DEFAULT_DELIMITER, '#', '%', '*', '|', '!', '&', '+', ';', '\\', '~',
    };
    private static final char DELIMITER_UNDEFINED = 0;

    private char delimiter = DELIMITER_UNDEFINED;

    private static final int NUMBER = 0, LETTER = 1;

    public static final Object[][] DATATYPE_LOOKUP_TABLE = {
        {new Integer(CFCSDatatype.ASCII), "A"},
        {new Integer(CFCSDatatype.FLOAT), "F"},
        {new Integer(CFCSDatatype.DOUBLE), "D"},
        {new Integer(CFCSDatatype.BINARY), "I"},
    };

    private static final Object[][] MODE_LOOKUP_TABLE = {
        {new Integer(CFCSData.LISTMODE), "L"},
        {new Integer(CFCSData.CORRELATED), "C"},
        {new Integer(CFCSData.UNCORRELATED), "U"},
    };

    private final CFCSKeywordTable keywords = new CFCSKeywordTable();

    // --------------------------------------------------------------------

    private static final /* inner */ class CFCSKeywordTable
    {
        private final List<CFCSKeyword> list = new LinkedList<CFCSKeyword>();
        private final Map<String, CFCSKeyword> map = new HashMap<String, CFCSKeyword>();

        // ------------------------------------------------------------

        final int size()
        {
            return list.size();
        }

        // ------------------------------------------------------------

        final boolean containsKeyword(final String name)
        {
            return map.containsKey(name.toUpperCase());
        }

        // ------------------------------------------------------------

        final CFCSKeyword getKeyword(final String name)
        {
            return map.get(name.toUpperCase());
        }

        final CFCSKeyword getKeyword(final int index)
        {
            return list.get(index);
        }

        // ------------------------------------------------------------

        final void addKeyword(final CFCSKeyword keyword)
        {
            map.put((keyword.getKeywordName()).toUpperCase(), keyword);
            list.add(keyword);
        }

        // ------------------------------------------------------------

        final void deleteKeyword(final String name)
        {
            list.remove(map.remove(name.toUpperCase()));
        }
    }

    // --------------------------------------------------------------------
    // Initial methods that follow actually know where keywords are stored.
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Returns the total number of keywords currently defined classesbined
    // between TEXT, STEXT and ANALYSIS (and possibly OTHER segments.)

    public final int getCount()
    {
        return keywords.size();
    }

    // --------------------------------------------------------------------
    // Returns the keyword named "name".  (Note: the instance returned
    // is a "copy" of the actual keyword).

    public final CFCSKeyword getKeyword(final String name)
    {
        if (name == null || name.length() == 0)
        { /* 3.2.9 */
            throw new CFCSError(CFCSIllegalName);
        }

        if (keywords.containsKeyword(name))
            return (keywords.getKeyword(name)).copy();
        /* else */ throw new CFCSError(CFCSKeywordNotFound, name);
    }

    // --------------------------------------------------------------------
    // Returns the keyword at index "index".  (Note: the instance returned
    // is a "copy" of the actual keyword).

    public final CFCSKeyword getKeyword(final int index)
    {
        if (index < 0 || index > keywords.size() - 1)
        {
            throw new CFCSError(CFCSBadKeywordIndex, index);
        }

        return keywords.getKeyword(index).copy();
    }

    // --------------------------------------------------------------------
    // Adds the keyword to the current list.  If the keyword already exists,
    // then it is replaced by that specified in the parameter (i.e.,
    // "replaceKeyword" functionality.)

    public final void addKeyword(final CFCSKeyword keyword)
    {
        final String name = keyword.getKeywordName();

        if (SYSTEM_KEYWORDS.contains(name))
        {
            throw new CFCSError(CFCSCannotModifyValue, name);
        }
        else if (CFCSParameters.isParameter(name) || CFCSGatingParameters.isParameter(name))
        {
            throw new CFCSError(CFCSCannotModifyValue, name);
        }
        else if (!keywords.containsKeyword(name))
        {

            String value = keyword.getKeywordValue();
            int segIdx = keyword.getKeywordSource();

            keywords.addKeyword(new CFCSKeyword(name, value, segIdx));
        }
        else
            replaceKeyword(keyword);
    }
    
    public static Boolean canAdd(final CFCSKeyword keyword) {
    	return canAdd(keyword.getKeywordName());
    }
    
    public static Boolean canAdd(final String keywordName) {
    	if(
    		SYSTEM_KEYWORDS.contains(keywordName) || 
    		CFCSParameters.isParameter(keywordName) || 
    		CFCSGatingParameters.isParameter(keywordName)) 
    		return false;
    	else return true;
    }	

    // --------------------------------------------------------------------
    // Deletes the keyword corresponding to "name".
    // If the keyword doesn't exist, then an exception is thrown.

    public final void deleteKeyword(final String name)
    {
        if (name == null || name.length() == 0)
        { /* 3.2.9 */
            throw new CFCSError(CFCSIllegalName);
        }
        else if (SYSTEM_KEYWORDS.contains(name))
        {
            throw new CFCSError(CFCSCannotModifyValue, name);
        }
        else if (keywords.containsKeyword(name))
        {
            keywords.deleteKeyword(name);
        }
        else
            throw new CFCSError(CFCSKeywordNotFound, name);
    }

    // --------------------------------------------------------------------
    // Replaces the keyword with that specified by the parameter.  If the
    // keyword doesn't exist, then an exception is thrown.

    public final void replaceKeyword(final CFCSKeyword keyword)
    {
        final String name = keyword.getKeywordName();

        if (SYSTEM_KEYWORDS.contains(name))
        {
            throw new CFCSError(CFCSCannotModifyValue, name);
        }
        else if (keywords.containsKeyword(name))
        {
            CFCSKeyword existing = keywords.getKeyword(name);
            int segIdx = keyword.getKeywordSource();

            if (segIdx != existing.getKeywordSource())
            {
                throw new CFCSError(CFCSCannotModifySource, segIdx);
            }
            /* else */ existing.setKeywordValue(keyword.getKeywordValue());
        }
        else
            throw new CFCSError(CFCSKeywordNotFound, name);
    }

    // --------------------------------------------------------------------
    // Friendly methods that bypass safety checks and avoid copying.
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Fast, non-copying keyword lookup routines for friends who need to
    // loop over lots of keywords quickly.

    /* friendly */
    final CFCSKeyword getSystemKeyword(final String name)
    {

        if (keywords.containsKeyword(name))
            return keywords.getKeyword(name);

        return null;
    }

    /* friendly */
    final CFCSKeyword getSystemKeyword(final int index)
    {
        return keywords.getKeyword(index);
    }

    // --------------------------------------------------------------------
    // Adds the keyword to the current list.  If the keyword already exists,
    // then replaceSystemKeyword() is called instead.  Doesn't copy data and
    // allows any keyword to be set.

    /* friendly */
    final void addSystemKeyword(final CFCSKeyword keyword)
    {
        if (!keywords.containsKeyword(keyword.getKeywordName()))
        {
            keywords.addKeyword(keyword);
        }
        else
            replaceSystemKeyword(keyword);
    }
    
    public void setFileKeyword(String value) {
    	CFCSKeyword keyword = new CFCSKeyword();
    	keyword.setKeywordName(CFCSKeywords.FILE_KEYWORD);
    	keyword.setKeywordValue(value);
    	keyword.setKeywordSource(CFCSDataSet.TEXT);
    	replaceSystemKeyword(keyword);
    }

    // --------------------------------------------------------------------
    // Replaces the keyword with that specified by the parameter.  If the
    // keyword doesn't exist, then an exception is thrown.  Allows any
    // value to be reset.

    /* friendly */
    final void replaceSystemKeyword(final CFCSKeyword keyword)
    {
        final String name = keyword.getKeywordName();

        if (keywords.containsKeyword(name))
        {
            CFCSKeyword existing = keywords.getKeyword(name);
            int segIdx = keyword.getKeywordSource();

            if (segIdx != existing.getKeywordSource())
            {
                throw new CFCSError(CFCSCannotModifySource, segIdx);
            }
            /* else */ existing.setKeywordValue(keyword.getKeywordValue());
        }
        else
            throw new CFCSError(CFCSKeywordNotFound, name);
    }

    // --------------------------------------------------------------------
    // Methods below don't know where keywords are stored and use the above.
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // Escape (double) delimiters found inside keyword name or value and
    // put a guard character in front of delimiter found at beginning of
    // keyword name or value.  If none found, returns original string.

    private static String escape(final String string, final char delimiter)
    {
        final int length = string.length();
	    int start = 0;
	    int end;
	    final StringBuffer buffer = new StringBuffer();

        if (string.charAt(0) == delimiter)
            buffer.append(DELIMITER_GUARD);

        while (start < length && (end = string.indexOf(delimiter, start)) >= 0)
        {
            (buffer.append(string.substring(start, ++end))).append(delimiter);
            start = end;
        }

        if (start == 0)
            return string;
        /* else */ if (start < length)
            buffer.append(string.substring(start, length));

        return buffer.toString();
    }

    // --------------------------------------------------------------------
    // Pick delimiter that doesn't conflict with keyword names and values if possible.
    // Save it, as we're required to reuse the same one for TEXT, ANALYSIS and STEXT

    private char getDelimiter()
    {
        if (delimiter != DELIMITER_UNDEFINED)
            return delimiter;

        int count = getCount();

        outer: for (int j = 0; j < delimiters.length; j++)
        {
            char character = delimiters[j];

            for (int i = 0; i < count; i++)
            {
                CFCSKeyword keyword = getSystemKeyword(i);

                if ((keyword.getKeywordName()).indexOf(character) > -1 ||
                        (keyword.getKeywordValue()).indexOf(character) > -1)
                {
                    continue outer;
                }
            }

            return (delimiter = character);
        }

        return (delimiter = DEFAULT_DELIMITER);
    }

    // --------------------------------------------------------------------
    // Load data into the keywords from bytes found in a file.

    /* friendly */
    final void setBytes(final byte[] buffer, final int segIdx)
    {
        final char delimiter = (char) buffer[0];
        boolean lastWasDelimiter = false;
        boolean readingKey = true;

        StringBuffer key = new StringBuffer();
        StringBuffer value = new StringBuffer();

        final int length = buffer.length;

        for (int i = 1; i < length; i++)
        {
            char ch = (char) buffer[i];

            if (ch == delimiter)
            {
                if ((lastWasDelimiter = !lastWasDelimiter) == false)
                { // double delimiter == escaped delimiter
                    ((readingKey) ? key : value).append(ch); // put a delimiter in the buffer
                }
            }
            else
            {
                if (lastWasDelimiter)
                {
                    if ((readingKey = !readingKey) == true)
                    {
                    	// Added UTF8 support
                    	// JS, 2008-08-17
                    	// Getting the chars from the String buffer, converting these to bytes and the bytes to
                    	// a UTF8 encoded String
                    	char chars[] = new char[value.length()];
                    	byte bytes[] = new byte[chars.length];
                    	value.getChars(0, value.length(), chars, 0);
                    	for(int j = 0; j < chars.length; j++) bytes[j] = (byte)chars[j];
                    	// addSystemKeyword(new CFCSKeyword(new String(key), new String(value), segIdx));
                    	try { 
                    		addSystemKeyword(new CFCSKeyword(new String(key), new String(bytes, "UTF8"), segIdx));
						} catch (UnsupportedEncodingException e) {
							e.printStackTrace();
						}
                    	
                        key = new StringBuffer(); // clear out the buffers
                        value = new StringBuffer();
                    }
                    lastWasDelimiter = false;
                }

                ((readingKey) ? key : value).append(ch);
            }
        }

        if (lastWasDelimiter && !readingKey)
        {	
        	// check for encountering a delimiter at the very end
        	
        	// Added UTF8 support
        	// JS, 2008-08-17
        	// Getting the chars from the String buffer, converting these to bytes and the bytes to
        	// a UTF8 encoded String
        	char chars[] = new char[value.length()];
        	byte bytes[] = new byte[chars.length];
        	value.getChars(0, value.length(), chars, 0);
        	for(int j = 0; j < chars.length; j++) bytes[j] = (byte)chars[j];
        	// addSystemKeyword(new CFCSKeyword(new String(key), new String(value), segIdx));
        	try { 
        		addSystemKeyword(new CFCSKeyword(new String(key), new String(bytes, "UTF8"), segIdx));
			} catch (UnsupportedEncodingException e) {
				e.printStackTrace();
			}
        }
    }

    // --------------------------------------------------------------------
    // Unload data from the keywords in a form that can be stored to file.
    // Returns null if there are no keywords to store in the segIdx segment.

    /* friendly */
    final byte[] getBytes(final int segIdx)
    {
        final char delimiter = getDelimiter();

        final StringBuffer buffer = new StringBuffer();

        final int count = getCount();

        for (int i = 0; i < count; i++)
        {
            CFCSKeyword keyword = getSystemKeyword(i);

            if (keyword.getKeywordSource() == segIdx)
            {
                String name = escape(keyword.getKeywordName(), delimiter);
                String value = escape(keyword.getKeywordValue(), delimiter);

                (((buffer.append(delimiter)).append(name)).append(delimiter)).append(value);
            }
        }

        if (buffer.length() == 0)
            return null;

        buffer.append(delimiter);

        return (buffer.toString()).getBytes();
    }

    // --------------------------------------------------------------------
    // Similar to getBytes() but only for the TEXT segment.  Separately
    // handles FCS-defined keywords ($XYZ) from other keywords when the
    // TEXT segment needs to overflow into an STEXT segment.

    /* friendly */
    final byte[] getTextBytes(final boolean want_fcs_defined)
    {
        final int count = getCount();
        final char delimiter = getDelimiter();
        final StringBuffer buffer = new StringBuffer();

        for (int i = 0; i < count; i++)
        {
            CFCSKeyword keyword = getSystemKeyword(i);

            if (keyword.getKeywordSource() != CFCSDataSet.TEXT)
                continue;

            boolean is_fcs_defined = keyword.isFCSDefined();

            if ((want_fcs_defined && !is_fcs_defined) || (!want_fcs_defined && is_fcs_defined))
                continue;

            String name = escape(keyword.getKeywordName(), delimiter);
            String value = escape(keyword.getKeywordValue(), delimiter);

            (((buffer.append(delimiter)).append(name)).append(delimiter)).append(value);
        }

        if (buffer.length() == 0)
            return null;

        buffer.append(delimiter);

        return (buffer.toString()).getBytes();
    }

    // --------------------------------------------------------------------
    // Front-end onto the $DATATYPE keyword, substitutes a program constant

    /* friendly */
    final int getDatatype()
    {
        String value = null;

        try
        {
            value = (getKeyword(DATATYPE_KEYWORD)).getKeywordValue();
        }
        catch (CFCSError error)
        {
            return CFCSDatatype.UNDEFINED;
        }

        for (int i = 0; i < DATATYPE_LOOKUP_TABLE.length; i++)
        {
            if (value.equalsIgnoreCase((String) DATATYPE_LOOKUP_TABLE[i][LETTER]))
            {
                int type = ((Integer) DATATYPE_LOOKUP_TABLE[i][NUMBER]).intValue();
                return type;
            }
        }

        return CFCSDatatype.UNDEFINED;
    }

    // --------------------------------------------------------------------

    /* friendly */
    final void setDatatype(final int type)
    {
        for (int i = 0; i < DATATYPE_LOOKUP_TABLE.length; i++)
        {
            if (type == ((Integer) DATATYPE_LOOKUP_TABLE[i][NUMBER]).intValue())
            {
                String value = (String) DATATYPE_LOOKUP_TABLE[i][LETTER];
                addSystemKeyword(new CFCSKeyword(DATATYPE_KEYWORD, value));
                return;
            }
        }

        throw new CFCSError(CFCSIllegalValue, type);
    }

    // --------------------------------------------------------------------
    // Front-end onto the $MODE keyword, substitutes a program constant

    /* friendly */
    final int getMode()
    {
        String value = null;

        try
        {
            value = (getKeyword(MODE_KEYWORD)).getKeywordValue();
        }
        catch (CFCSError error)
        {
            return CFCSData.UNDEFINED;
        }

        for (int i = 0; i < MODE_LOOKUP_TABLE.length; i++)
        {
            if (value.equalsIgnoreCase((String) MODE_LOOKUP_TABLE[i][LETTER]))
            {
                return ((Integer) MODE_LOOKUP_TABLE[i][NUMBER]).intValue();
            }
        }

        return CFCSData.UNDEFINED;
    }

    // --------------------------------------------------------------------

    void setMode(final int mode)
    {
        for (int i = 0; i < MODE_LOOKUP_TABLE.length; i++)
        {
            if (mode == ((Integer) MODE_LOOKUP_TABLE[i][NUMBER]).intValue())
            {
                String value = (String) MODE_LOOKUP_TABLE[i][LETTER];
                addSystemKeyword(new CFCSKeyword(MODE_KEYWORD, value));
                return;
            }
        }

        throw new CFCSError(CFCSIllegalValue, mode);
    }

    // --------------------------------------------------------------------
    // Convert from $BYTEORD string value to an array

    protected final int[] getByteSwapArray()
    {
        final CFCSKeyword keyword = getSystemKeyword(BYTEORDER_KEYWORD);

        final String string = keyword.getKeywordValue();

        final List<Integer> positions = new LinkedList<Integer>();

        for (int i = 0; i < string.length(); i++)
        {
            char character = string.charAt(i);

            if (Character.isDigit(character))
            {
                positions.add(new Integer(Character.digit(character, 10)));
            }
        }

        final int[] array = new int[positions.size()];

        for (int i = 0; i < array.length; i++)
            array[i] = positions.get(i).intValue();

        return array;
    }

    // --------------------------------------------------------------------
    // Convert from an array to a string value for $BYTEORD

    protected final void setByteSwapArray(final int[] array)
    {
        final StringBuffer buffer = new StringBuffer();

        for (int i = 0; i < array.length; i++)
        {
            buffer.append(array[i]);

            if (i < (array.length - 1))
                buffer.append(CFCSSystem.VALUE_SEPARATOR_CHAR);
        }

        addSystemKeyword(new CFCSKeyword(BYTEORDER_KEYWORD, buffer.toString()));
    }

    // --------------------------------------------------------------------

    
    /*
     * JS: Do more of these? All FCS defined keywords?
     */
    public final double getSampleVolume()
    {
        String value = null;

        try
        {
            value = (getKeyword(SAMPLE_VOLUME_KEYWORD)).getKeywordValue();
        }
        catch (CFCSError error)
        {
            return Double.NaN;
        }

        try 
        { 
        	return Double.parseDouble(value); 
        } 
        catch (Exception e) 
        { 
        	return Double.NaN; 
        }
    }
    
    /*
     * JS: Do more of these? All FCS defined keywords?
     */
    public final void setSampleVolume(double volume)
    {
    	addSystemKeyword(new CFCSKeyword(SAMPLE_VOLUME_KEYWORD, volume + ""));
    }
    
    // JS
    public final String getPlateId()
    {
    	String value = null;
    	try
        {
            value = (getKeyword(PLATEID_KEYWORD)).getKeywordValue();
            return value;
        }
        catch (CFCSError error)
        {
            return null;
        }
    }

    // JS
    public final void setPlateId(String plateId)
    {
    	addSystemKeyword(new CFCSKeyword(PLATEID_KEYWORD, plateId));
    }
    
    // JS
    public final String getPlateName()
    {
    	String value = null;
    	try
        {
            value = (getKeyword(PLATENAME_KEYWORD)).getKeywordValue();
            return value;
        }
        catch (CFCSError error)
        {
            return null;
        }
    }

    // JS
    public final void setPlateName(String plateName)
    {
    	addSystemKeyword(new CFCSKeyword(PLATENAME_KEYWORD, plateName));
    }
    
    // JS
    public final String getWellId()
    {
    	String value = null;
    	try
        {
            value = (getKeyword(WELLID_KEYWORD)).getKeywordValue();
            return value;
        }
        catch (CFCSError error)
        {
            return null;
        }
    }

    // JS
    public final void setWellId(String wellId)
    {
    	addSystemKeyword(new CFCSKeyword(WELLID_KEYWORD, wellId));
    }
    
    // JS
    public final CFCSSpillover getSpillover()
    {
    	String value = null;
    	
    	try
        {
            value = (getKeyword(SPILLOVER_KEYWORD)).getKeywordValue();
        }
        catch (CFCSError error)
        {
        	try
            {
        		value = (getKeyword(SPILLOVER_ALTERNATIVE_KEYWORD)).getKeywordValue();
            } catch (CFCSError err) 
            {
            	return null;
            }
        }
        
        return new CFCSSpillover(value);
    }
    
    // JS
    public final void setSpillover(CFCSSpillover spillover)
    {
    	addSystemKeyword(new CFCSKeyword(SPILLOVER_KEYWORD, spillover.getSpilloverKeywordValue()));
    }
    
    // JS
    public final OriginalityEnum getOriginality()
    {
    	String value = null;
    	OriginalityEnum ret = OriginalityEnum.Undefined;
    	try
        {
            value = (getKeyword(ORIGINALITY_KEYWORD)).getKeywordValue();
            if(value.compareToIgnoreCase("Original") == 0) ret = OriginalityEnum.Original;
            else if(value.compareToIgnoreCase("NonDataModified") == 0) ret = OriginalityEnum.NonDataModified;
            else if(value.compareToIgnoreCase("Appended") == 0) ret = OriginalityEnum.Appended;
            else if(value.compareToIgnoreCase("DataModified") == 0) ret = OriginalityEnum.DataModified;
        }
        catch (CFCSError error) {}
        return ret;
    }

    // JS
    public final void setOriginality(OriginalityEnum originality)
    {
    	String value = null;
    	switch(originality) {
    	case Original: value = "Original"; break;
    	case NonDataModified: value = "NonDataModified"; break;
    	case Appended: value = "Appended"; break;
    	case DataModified: value = "DataModified"; break;
    	}
    	if(value != null) addSystemKeyword(new CFCSKeyword(ORIGINALITY_KEYWORD, value));
    }
    
    // JS
    public final String getLastModifier()
    {
    	String value = null;
    	try
        {
            value = (getKeyword(LAST_MODIFIER_KEYWORD)).getKeywordValue();
            return value;
        }
        catch (CFCSError error)
        {
            return null;
        }
    }

    // JS
    public final void setLastModifier(String lastModifier)
    {
    	addSystemKeyword(new CFCSKeyword(LAST_MODIFIER_KEYWORD, lastModifier));
    }
    
    // JS
    public final Date getLastModified()
    {
    	try
        {
    		String value = null;
        	Calendar cal = null;
        	value = (getKeyword(LAST_MODIFIED_KEYWORD)).getKeywordValue();
            // e.g., 25-SEP-2008 15:22:10.47
            String s[] = value.split(" ", 2);
            if(s.length != 2) throw new CFCSError(CFCSErrorCodes.CFCSIllegalValue, value);
            String sd[] = s[0].split("-", 3);
            if(sd.length != 3) throw new CFCSError(CFCSErrorCodes.CFCSIllegalValue, value);
            int year = Integer.parseInt(sd[2]);
            int month = monthStrtoNum(sd[1]);
            int date = Integer.parseInt(sd[0]);
            String sh[] = s[1].split(":", 3);
            if(sh.length != 3) throw new CFCSError(CFCSErrorCodes.CFCSIllegalValue, value);
            int hrs = Integer.parseInt(sh[0]);
            int min = Integer.parseInt(sh[1]);
            int sec = (int) Double.parseDouble(sh[2]);
            int msec = (int)(100 * (Double.parseDouble(sh[2]) - (double)sec));
            
            cal = Calendar.getInstance();
            cal.set(Calendar.YEAR, year);
            cal.set(Calendar.MONTH, month - 1);
            cal.set(Calendar.DATE, date);
            cal.set(Calendar.HOUR_OF_DAY, hrs);
            cal.set(Calendar.MINUTE, min);
            cal.set(Calendar.SECOND, sec);
            cal.set(Calendar.MILLISECOND, msec);
            
            return cal.getTime();
        }
        catch (Exception error)
        {
            return null;
        }
    }
    
    // JS
    public final void setLastModified(Date d)
    {
    	SimpleDateFormat df = new SimpleDateFormat("dd-MMM-yyyy HH:mm:ss.SS");
    	String value = df.format(d).toUpperCase();
    	if(value != null) addSystemKeyword(new CFCSKeyword(LAST_MODIFIED_KEYWORD, value));
    }
    	
    
    private static int monthStrtoNum(String month) {
    	if(month.compareToIgnoreCase("JAN") == 0) return 1;
    	else if(month.compareToIgnoreCase("FEB") == 0) return 2;
    	else if(month.compareToIgnoreCase("MAR") == 0) return 3;
    	else if(month.compareToIgnoreCase("APR") == 0) return 4;
    	else if(month.compareToIgnoreCase("MAY") == 0) return 5;
    	else if(month.compareToIgnoreCase("JUN") == 0) return 6;
    	else if(month.compareToIgnoreCase("JUL") == 0) return 7;
    	else if(month.compareToIgnoreCase("AUG") == 0) return 8;
    	else if(month.compareToIgnoreCase("SEP") == 0) return 9;
    	else if(month.compareToIgnoreCase("OCT") == 0) return 10;
    	else if(month.compareToIgnoreCase("NOV") == 0) return 11;
    	else if(month.compareToIgnoreCase("DEC") == 0) return 12;
    	else return -1;
    }
    
}
