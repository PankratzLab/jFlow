/**
 * 
 */
package cnv.var;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import cnv.plots.CompPlot;

/**
 * This class reads in CNV files and generates a list of rectangles representing them. The list can be retrieved in different contexts based on the current display mode.
 * 
 * @author Michael Vieths
 * 
 */
public class CNVRectangles {
    private HashMap<String, ArrayList<CNVariant>> fileMap;
    private ArrayList<CNVRectangle>               cnvRectangles;
    public Color[]                                colorScheme;
    private int                                   rectangleHeight;
    private float                                 scalingFactor;
    int                                           lowestStart;
    int                                           probes;
    int                                           minSize;
    int                                           qualityScore;

    public CNVRectangles(String[] files, int[] location, int probes, int minSize, int qualityScore) {
        // Set the color scheme
        colorScheme = CompPlot.colorScheme;
        // Read the data from the CNV files

        fileMap = new HashMap<String, ArrayList<CNVariant>>();

        for (String file : files) {
            // Load the CNVs out of the files
            CNVariantHash cnvHash = CNVariantHash.load(file, CNVariantHash.CONSTRUCT_ALL, false);

            // All CNVs are loaded into an array, which is then stored in an array by file
            // Format is:
            // file1: CNV1..CNVN
            // fileN: CNV1..CNVN
            // All CNVs in each file, when they are rendered, will be a single color
            // Any CNVs where CNVariant.getCN() != 2 will be a different shade of that color
            CNVariant[] cnvs = cnvHash.getAllInRegion((byte) location[0], location[1], location[2], probes, minSize, qualityScore);
            ArrayList<CNVariant> cnvList = new ArrayList<CNVariant>(Arrays.asList(cnvs));
            fileMap.put(file, cnvList);
        }

        // Populate the hashmap and rectangles
        cnvRectangles = new ArrayList<CNVRectangle>();
        int i = 0;
        for (String key : fileMap.keySet()) {
            for (CNVariant variant : fileMap.get(key)) {
                // Set the color
                CNVRectangle cnvRect = new CNVRectangle(variant, location[1]);
                cnvRect.setCNVColor(colorScheme[i % colorScheme.length]);
                cnvRect.addCNV(variant);
                cnvRectangles.add(cnvRect);
            }
        }
    }

    /**
     * Create an empty object
     */
    public CNVRectangles() {
        cnvRectangles = new ArrayList<CNVRectangle>();
    }

    /**
     * Collapse the rectangles so that all CNVs with the same footprint (start and end) appear as a single rectangle
     * 
     * @return ArrayList of collapsed rectangles
     */
    public ArrayList<CNVRectangle> getCollapsedRectangles() {
        ArrayList<CNVRectangle> collapsedRectangles = new ArrayList<CNVRectangle>();
        HashMap<String, CNVRectangle> cnvMap = new HashMap<String, CNVRectangle>();

        /**
         * Store the rectangles in a HashMap based on their signature
         */
        for (CNVRectangle cnvRect : cnvRectangles) {
            String key = (int) cnvRect.getStartXValue() + "-" + (int) cnvRect.getStopXValue();
            if (cnvMap.containsKey(key)) {
                cnvMap.get(key).addCNV(cnvRect.getCNV());
            } else {
                cnvMap.put(key, cnvRect);
            }
        }

        /**
         * Associate all of the CNVs with the same signature into their own rectangles
         */
        int i = 0;
        for (String key : cnvMap.keySet()) {
            CNVRectangle cnvRect = cnvMap.get(key);
            int width = Math.round(((int) cnvRect.getStopXValue() - (int) cnvRect.getStartXValue()) * scalingFactor);
            int x = Math.round((int) cnvRect.getStartXValue() * scalingFactor);
            // Ensure a 1 pixel gap between CNVs
            int y = (i * rectangleHeight) + i;

            // Store the rectangle for later bounds checking
            cnvRect.setRect(x, y, width, rectangleHeight);
            collapsedRectangles.add(cnvRect);
            i++;
        }

        return collapsedRectangles;
    }

    /**
     * Return all of the rectangles with one Individual ID per line
     * 
     * @return ArrayList of all rectangles sorted by Individual ID
     */
    public ArrayList<CNVRectangle> getFullRectangles() {
        ArrayList<CNVRectangle> fullRectangles = new ArrayList<CNVRectangle>();

        /**
         * Store everything in a hashmap with a key of the individual ID so we can get all CNVs associated with that ID on the same line
         */
        HashMap<String, ArrayList<CNVRectangle>> cnvMap = new HashMap<String, ArrayList<CNVRectangle>>();
        for (CNVRectangle cnvRect : cnvRectangles) {
            String iid = cnvRect.getCNV().getIndividualID();
            if (cnvMap.containsKey(iid)) {
                cnvMap.get(iid).add(cnvRect);
            } else {
                ArrayList<CNVRectangle> cnvRects = new ArrayList<CNVRectangle>();
                cnvRects.add(cnvRect);
                cnvMap.put(iid, cnvRects);
            }
        }

        // Set the coordinates based on the line
        int i = 0;
        for (String key : cnvMap.keySet()) {
            ArrayList<CNVRectangle> cnvRects = cnvMap.get(key);
            for (CNVRectangle cnvRect : cnvRects) {
                int width = Math.round(((int) cnvRect.getStopXValue() - (int) cnvRect.getStartXValue()) * scalingFactor);
                int x = Math.round((int) cnvRect.getStartXValue() * scalingFactor);
                // Ensure a 1 pixel gap between CNVs
                int y = (i * rectangleHeight) + i;

                // Figure out what shade we should be
                cnvRect.setCNVColor(cnvRect.getCNVColor(), "Full");

                // Store the rectangle for later bounds checking
                cnvRect.setRect(x, y, width, rectangleHeight);
                fullRectangles.add(cnvRect);
            }
            i++;
        }

        return fullRectangles;
    }

    /**
     * Pack the rectangles with as few gaps as possible. Perform no sorting.
     * 
     * @return ArrayList of packed rectangles
     */
    public ArrayList<CNVRectangle> getPackedRectangles() {
        ArrayList<CNVRectangle> packedRectangles = new ArrayList<CNVRectangle>();
        clearUsed(cnvRectangles);
        lowestStart = getLowestRect(cnvRectangles);
        packedRectangles = packRectangles(cnvRectangles);
        return packedRectangles;
    }

    /**
     * Pack the provided rectangles, filling in line by line and setting the x,y coordinates
     * 
     * @param cnvRects
     * @return ArrayList of packed rectangles
     */
    private ArrayList<CNVRectangle> packRectangles(ArrayList<CNVRectangle> cnvRects) {
        ArrayList<CNVRectangle> packedRectangles = new ArrayList<CNVRectangle>();
        int i = 0;
        while (hasUnused(cnvRects)) {
            CNVRectangle cnvRect = getLeftMost(lowestStart, cnvRects);
            if (cnvRect != null) {
                // Set the Y coordinate
                int x = Math.round((int) cnvRect.getStartXValue() * scalingFactor);
                int y = (i * rectangleHeight) + i;
                int width = Math.round(((int) cnvRect.getStopXValue() - (int) cnvRect.getStartXValue()) * scalingFactor);
                cnvRect.setRect(x, y, width, rectangleHeight);

                // Figure out what shade we should be
                cnvRect.setCNVColor(cnvRect.getCNVColor(), "Pack");

                packedRectangles.add(cnvRect);

                // Mark this one as used so we don't add it again
                cnvRect.setUsed(true);
                do {
                    cnvRect = getLeftMost((int) cnvRect.getStopXValue(), cnvRects);
                    if (cnvRect == null) {
                        i++;
                        break;
                    } else {
                        x = Math.round((int) cnvRect.getStartXValue() * scalingFactor);
                        y = (i * rectangleHeight) + i;
                        width = Math.round(((int) cnvRect.getStopXValue() - (int) cnvRect.getStartXValue()) * scalingFactor);
                        cnvRect.setRect(x, y, width, rectangleHeight);

                        // Figure out what shade we should be
                        cnvRect.setCNVColor(cnvRect.getCNVColor(), "Full");

                        packedRectangles.add(cnvRect);
                        // Mark this one as used so we don't add it again
                        cnvRect.setUsed(true);
                    }
                } while (cnvRect != null);
            } else {
                break;
            }
        }

        return packedRectangles;
    }

    /**
     * Find the rectangle with the lowest X value (furthest left)
     * 
     * @param cnvRects
     * @return integer value of the lowest X value
     */
    private int getLowestRect(ArrayList<CNVRectangle> cnvRects) {
        int lowest = 0;

        for (CNVRectangle cnvRect : cnvRects) {
            int startX = (int) cnvRect.getStartXValue();
            if (startX < lowest) {
                lowest = startX;
            }
        }

        return lowest;
    }

    /**
     * Determine whether there are any rectangles we haven't looked at
     * 
     * @param cnvRects
     * @return true if any have not been used
     */
    private boolean hasUnused(ArrayList<CNVRectangle> cnvRects) {
        boolean unused = false;

        for (CNVRectangle cnvRect : cnvRects) {
            if (!cnvRect.isInUse()) {
                unused = true;
                break;
            }
        }

        return unused;
    }

    /**
     * Clear the used flag on all rectangles
     * 
     * @param cnvRects
     */
    private void clearUsed(ArrayList<CNVRectangle> cnvRects) {
        for (CNVRectangle cnvRect : cnvRects) {
            cnvRect.setUsed(false);
        }
    }

    /**
     * Get the next leftmost rectangle, starting at startX
     * 
     * @param startX
     * @param cnvRects
     * @return The leftmost rectangle
     */
    private CNVRectangle getLeftMost(int startX, ArrayList<CNVRectangle> cnvRects) {
        int lastPosition = startX + 2; // Want a 2-pixel buffer between CNVs
        if (startX == lowestStart) {
            // The first time through we'll be starting at lowestStart so we don't want the offset
            lastPosition = lowestStart;
        }

        CNVRectangle leftMostCNV = null;

        // Iterate through all of the retangles looking for any with a lower x
        for (CNVRectangle cnvRect : cnvRects) {
            if (((int) cnvRect.getStartXValue() >= lastPosition) && (!cnvRect.isInUse())) {
                if (leftMostCNV == null) {
                    leftMostCNV = cnvRect;
                } else {
                    int oldDifference = (int) leftMostCNV.getStartXValue() - lastPosition;
                    int newDifference = (int) cnvRect.getStartXValue() - lastPosition;
                    if (newDifference < oldDifference) {
                        leftMostCNV = cnvRect;
                    }
                }
            }
        }

        return leftMostCNV;
    }

    /**
     * Get the height of the rectangles
     * 
     * @return Height of the rectangles as an integer
     */
    public int getRectangleHeight() {
        return rectangleHeight;
    }

    /**
     * Set the height of the rectangles
     * 
     * @param rectangleHeight
     */
    public void setRectangleHeight(int rectangleHeight) {
        this.rectangleHeight = rectangleHeight;
    }

    /**
     * Get the scaling factor (based on the current viewing window)
     * 
     * @return scaling factor as a float
     */
    public float getScalingFactor() {
        return scalingFactor;
    }

    /**
     * Set the scaling factor (from CompPanel)
     * 
     * @param scalingFactor
     */
    public void setScalingFactor(float scalingFactor) {
        this.scalingFactor = scalingFactor;
    }

    /**
     * Get the number of probes on which we're filtering
     * 
     * @return the number of probes on which we're filtering
     */
    public int getProbes() {
        return probes;
    }

    /**
     * Set the number of probes on which to filter
     * 
     * @param The
     *            minimum number of probes
     */
    public void setProbes(int probes) {
        this.probes = probes;
    }

    /**
     * Get the smallest size of CNV on which we're filtering
     * 
     * @return the smallest size of CNV on which we're filtering
     */
    public int getMinSize() {
        return minSize;
    }

    /**
     * Set the minimum CNV size on which to filter
     * 
     * @param minSize
     *            Minimum CNV size
     */
    public void setMinSize(int minSize) {
        this.minSize = minSize;
    }

    /**
     * Get the minimum quality score on which we're filtering
     * 
     * @return the minimum quality score on which we're filtering
     */
    public int getQualityScore() {
        return qualityScore;
    }

    /**
     * Set the minimum quality score on which to filter
     * 
     * @param qualityScore
     *            The minimum quality score on which to filter
     */
    public void setQualityScore(int qualityScore) {
        this.qualityScore = qualityScore;
    }
}
