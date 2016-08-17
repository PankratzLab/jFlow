package org.flowcyt.cfcs;
// CFCSParameter.java

/* ------------------------------------------------------------------------- *\
This software and documentation are provided 'as is' and Tree Star, Inc., its
contractors and partners specifically disclaim all other warranties, expressed
or implied, including but not limited to implied warranties of merchantability
and fitness for a particular purpose, or during any particular date range.

By using this software, you are agreeing to these limits of liability, and to
hold Tree Star harmless for any information, accurate or erroneous, that might
be generated by the program.  This software is intended for research use only.

Christopher Lane <cdl@best.classes> for Tree Star  1/14/2002      Copyright 2002
\* ------------------------------------------------------------------------- */

public final class CFCSParameter extends CFCSAbstractParameter
{

    private int fieldSize = Integer.MIN_VALUE;
    private double power = Double.NaN;
    private double gain = Double.NaN;

    private String lambda = null;
    
    public enum PreferredDisplayScale { Linear, Logarithmic, Undefined };
    private String preferredDisplay = null;
    private PreferredDisplayScale preferredDisplayScale = PreferredDisplayScale.Undefined;
    private double[] preferredDisplayArgs = null; 
    
    private double calibrationFactor = Double.NaN;;
    private String calibrationUnit = null;    
    
    private static final double LOG10 = Math.log(10);

    /* friendly */
    static final int FIELDSIZE_FREE_FORMAT = 0;
    /* friendly */
    static final String FIELDSIZE_FREE_FORMAT_STRING = "*";

    // $PnB ---------------------------------------------------------------

    public final int getFieldSize()
    {
        if (isNotSet(fieldSize))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "FieldSize");
        }

        return fieldSize;
    }

    public final void setFieldSize(final int fieldSize)
    {
        if (fieldSize != FIELDSIZE_FREE_FORMAT && fieldSize < 1)
        {
            throw new CFCSError(CFCSIllegalValue, fieldSize);
        }
        /* else */ this.fieldSize = fieldSize;
    }

    // $PnB (bis) ---------------------------------------------------------
    // CFCS API interface to FieldSize is integer but attribute isn't really

    public final String getFieldSizeString()
    {
        if (fieldSize == FIELDSIZE_FREE_FORMAT)
            return FIELDSIZE_FREE_FORMAT_STRING;

        return Integer.toString(getFieldSize());
    }
    public final void setFieldSizeString(final String fieldSize)
    {
        if (fieldSize.equalsIgnoreCase(FIELDSIZE_FREE_FORMAT_STRING))
        {
            setFieldSize(FIELDSIZE_FREE_FORMAT);
        }
        else
        {
            try
            {
                setFieldSize((new Integer(fieldSize)).intValue());
            }
            catch (NumberFormatException exception)
            {
                throw new CFCSError(CFCSBadValueConversion, exception);
            }
        }
    }
    // $PnG ---------------------------------------------------------------
    public final double getGain()
    {
        if (isNotSet(gain))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "Gain");
        }

        return gain;
    }
    public final void setGain(final double gain)
    {
        this.gain = gain;
    }
    // $PnL ---------------------------------------------------------------
    // JS 2009-08-18
    // Modified for FCS3.1 where there may be more excitation wavelengths for a single parameter 
    // Set the string of excitation wavelengths
    public final void setExcitationWavelength(final String lambda)
    {
    	this.lambda = lambda;
    }

    // JS 2009-08-18
    // Modified for FCS3.1 where there may be more excitation wavelengths for a single parameter 
    // Get the string of excitation wavelengths
    public final String getExcitationWavelength()
    {
    	if (isNotSet(lambda))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "ExcitationWavelength");
        }

        return lambda;
    }
    
    // JS 2009-08-18
    // Added to support FCS3.1 where there may be more excitation wavelengths for a single parameter 
    // Return only the first excitation wavelength if there is more.
    public final int getSingleExcitationWavelength()
    {
    	if (isNotSet(lambda))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "ExcitationWavelength");
        }

    	return Integer.parseInt(lambda.split(",")[0]);
    }
    
    // JS 2009-08-18
    // Added to support FCS3.1 where there may be more excitation wavelengths for a single parameter 
    // Return only the first excitation wavelength if there is more.
    public final int[] getExcitationWavelengthArray()
    {
    	if (isNotSet(lambda))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "ExcitationWavelength");
        }

    	String[] sa = lambda.split(",");
    	int[] ret = new int[sa.length];
    	for(int i = 0; i < sa.length; i++) ret[i] = Integer.parseInt(sa[i]);

    	return ret;
    }
	
    // JS 2009-08-18
    // Modified for FCS3.1 where there may be more excitation wavelengths for a single parameter
    // Only set a single excitation wavelength
    public final void setExcitationWavelength(final int lambda)
    {
        this.lambda = new String(new Integer(lambda).toString());
    }
    
    // $PnO ---------------------------------------------------------------
    public final double getLaserPower()
    {
        if (isNotSet(power))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "LaserPower");
        }

        return power;
    }

    public final void setLaserPower(final double power)
    {
        if (power < 0)
        {
            throw new CFCSError(CFCSIllegalValue, power);
        }
        /* else */ this.power = power;
    }
    
    // $PnD ---------------------------------------------------------------
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for the preferred display 
    public final String getPreferredDisplay()
    {
        if (isNotSet(preferredDisplay))
        {
            throw new CFCSError(CFCSUndefinedAttribute, "PreferredDisplay");
        }

        return preferredDisplay;
    }
    
    // $PnD ---------------------------------------------------------------
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for the preferred display
    // $PnD/string,f1,f2/, e.g., $P3D/Linear,0,1024/ or $P2D/Logarithmic,4,0.1/
    public final void setPreferredDisplay(final String preferredDisplay)
    {
    	String s[] = preferredDisplay.split(",");
    	if(s.length != 3) throw new CFCSError(CFCSIllegalValue, preferredDisplay);
    	if(s[0].compareToIgnoreCase("Linear") == 0) preferredDisplayScale = PreferredDisplayScale.Linear;
    	else if(s[0].compareToIgnoreCase("Logarithmic") == 0) preferredDisplayScale = PreferredDisplayScale.Logarithmic;
    	else throw new CFCSError(CFCSIllegalValue, preferredDisplay);
    	
    	preferredDisplayArgs = new double[2];
    	preferredDisplayArgs[0] = Double.parseDouble(s[1]); 
	    preferredDisplayArgs[1] = Double.parseDouble(s[2]);
    	this.preferredDisplay = preferredDisplay;
    }
    
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for the preferred display 
    public final PreferredDisplayScale getPreferredDisplayScale() {
    	return preferredDisplayScale;
    }
    
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for the preferred display 
    public final double[] getPreferredDisplayArgs() {
    	if (isNotSet(preferredDisplayArgs)) {
            throw new CFCSError(CFCSUndefinedAttribute, "PreferredDisplay");
        }
    	if (isNotSet(preferredDisplayScale)) {
    		throw new CFCSError(CFCSUndefinedAttribute, "PreferredDisplay");
    	}
    	return preferredDisplayArgs;
    }
    
    
    static final boolean isNotSet(final PreferredDisplayScale preferredDisplayScale) {
          return preferredDisplayScale == PreferredDisplayScale.Undefined;
	}
    
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for calibration 
    public final String getCalibration() {
    	 if (isNotSet(calibrationFactor)) {
             throw new CFCSError(CFCSUndefinedAttribute, "Calibration");
         }

         if(isNotSet(calibrationUnit)) return calibrationFactor + "";
         else return calibrationFactor + "," + calibrationUnit;
    }

    // JS 2009-08-18
    // Added for FCS3.1; new keyword for calibration 
    public final void setCalibration(final String calibration) {
    	String sa[] = calibration.split(",", 2);
    	try {
    		calibrationFactor = Double.parseDouble(sa[0]);
    	} catch (Exception e) {
    		throw new CFCSError(CFCSIllegalValue, calibration);
    	}
    	
    	if(sa.length == 2) this.calibrationUnit = sa[1];
    }
    
    // JS 2009-08-18
    // Added for FCS3.1; new keyword for calibration 
    public final double getCalibrationFactor() {
   	 if (isNotSet(calibrationFactor)) {
            throw new CFCSError(CFCSUndefinedAttribute, "Calibration");
        }
        else return calibrationFactor;
    }
    
 	// JS 2009-08-18
    // Added for FCS3.1; new keyword for calibration 
    public final String getCalibrationUnit() {
   	 if (isNotSet(calibrationUnit)) {
            throw new CFCSError(CFCSUndefinedAttribute, "Calibration Unit");
        }
        else return calibrationUnit;
    }
    
    
    // $PnE ---------------------------------------------------------------
    // Overide $xnE version in CFCSAbstractParameter to add error checking.
    // Zero implies linear scaling, in which case the gain must be nonzero.
	public final void setLogDecades(final double decades)
    {
        if (decades == 0.0 && isSet(gain) && gain == 0.0)
        {
            throw new CFCSError(CFCSInconsistentValue, "Gain");
        }
        /* else */ super.setLogDecades(decades);
    }

    // --------------------------------------------------------------------
    // Using scaling information, converts the value corresponding to "scale"
    // into a channel number.  Out-of-range values are converted.

    public final int scaleToChannel(double scale)
    {
        final double channel;
	    final double decades = getLogDecades();

	    if (decades > 0.0)
        { // logarithmic
            double offset = getOffset();

            if (offset > 0.0)
                scale /= offset;

            channel = (Math.log(scale) / LOG10) * (getRange() / decades);
        }
        else
        { // linear
            channel = (isSet(gain) && gain > 0.0) ? scale * gain : scale;
        }

        return Math.round((float) channel);
    }

    // --------------------------------------------------------------------
    // Using scaling information, converts a channel number into a scaled
    // value.  Out-of-range values are converted.

    public final double channelToScale(final int channel)
    {
        final double decades = getLogDecades();

        if (decades > 0.0)
        { // logarithmic
        	
            double scale = Math.pow(10.0, (channel * decades / getRange()));
            double offset = getOffset();

            if (offset > 0.0)
                return (scale * offset);
            /* else */ return scale;
        }
        else
        { // linear
            return (isSet(gain) && gain > 0.0) ? channel / gain : channel;
        }
    }

    // --------------------------------------------------------------------

}
