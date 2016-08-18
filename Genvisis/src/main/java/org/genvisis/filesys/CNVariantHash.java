package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JProgressBar;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.SerializedFiles;

public class CNVariantHash implements Serializable {
    public static final long                                  serialVersionUID = 1L;

    public static final int                                   CONSTRUCT_BY_IND = 1;
    public static final int                                   CONSTRUCT_ALL    = 2;

    private String                                            filename;             // The file this hash came from
    private Hashtable<String, Hashtable<String, CNVariant[]>> hashes;

    public CNVariantHash(String filename, int structureType, boolean jar, Logger log) {
        Hashtable<String, Hashtable<String, Vector<CNVariant>>> vHashes;
        Hashtable<String, Vector<CNVariant>> vHash;
        Vector<CNVariant> v;
        Hashtable<String, CNVariant[]> finalHash;
        CNVariant cnv;
        String trav, temp;
        String[] inds, chrs;
        // long time;
        BufferedReader reader;
        String[] line;
        ProgressMonitor progMonitor = new ProgressMonitor(new JProgressBar(), log);
        int count;

        setFilename(filename);
        String taskName = "CNV_CONVERSION";
        progMonitor.beginDeterminateTask(taskName, "Converting CNVs", Files.countLines(filename, 1), ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);

        // time = new Date().getTime();
        vHashes = new Hashtable<String, Hashtable<String, Vector<CNVariant>>>();
        count = 0;
        try {
            reader = Files.getReader(filename, jar, true, true);

            reader.mark(1000);
            temp = reader.readLine();
            count += temp.length();
            line = temp.trim().split("[\\s]+");
            if (!line[2].toLowerCase().equals("chr") && Positions.chromosomeNumber(line[2]) != -1) {
                reader.reset();
            }
            while (reader.ready()) {
                temp = reader.readLine();
                count += temp.length();
                progMonitor.updateTask(taskName);
                cnv = new CNVariant(temp.trim().split("[\\s]+"));
                trav = cnv.getFamilyID() + "\t" + cnv.getIndividualID();
                if (structureType == CONSTRUCT_ALL) {
                    // cnv.setFamilyID(null); // CompPanel will still need to link to the proper IDs
                    // cnv.setIndividualID(null);
                    trav = "all";
                }

                if (vHashes.containsKey(trav)) {
                    vHash = vHashes.get(trav);
                } else {
                    vHashes.put(trav, vHash = new Hashtable<String, Vector<CNVariant>>());
                }
                if (vHash.containsKey(cnv.getChr() + "")) {
                    v = vHash.get(cnv.getChr() + "");
                } else {
                    vHash.put(cnv.getChr() + "", v = new Vector<CNVariant>());
                }
                v.add(cnv);

            }
            reader.close();
        } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + filename + "\" not found in current directory");
            fnfe.printStackTrace();
        } catch (IOException ioe) {
            System.err.println("Error reading file \"" + filename + "\"");
            ioe.printStackTrace();
        }
        progMonitor.endTask(taskName);

        hashes = new Hashtable<String, Hashtable<String, CNVariant[]>>();
        // time = new Date().getTime();

        inds = HashVec.getKeys(vHashes);
        taskName = "CNV_SERIALIZATION";
        progMonitor.beginDeterminateTask(taskName, "Serializing CNVs", inds.length, ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
        for (int i = 0; i < inds.length; i++) {
        	progMonitor.updateTask(taskName);
            vHash = vHashes.get(inds[i]);
            finalHash = new Hashtable<String, CNVariant[]>();
            chrs = HashVec.getKeys(vHash);
            for (int j = 0; j < chrs.length; j++) {
                finalHash.put(chrs[j], CNVariant.toCNVariantArray(vHash.get(chrs[j])));
            }
            hashes.put(inds[i], finalHash);
        }
        progMonitor.endTask(taskName);
    }

    /**
     * @return the filename
     */
    public String getFilename() {
        return filename;
    }

    /**
     * @param filename
     *            the filename to set
     */
    public void setFilename(String filename) {
        this.filename = filename;
    }

    public Hashtable<String, CNVariant[]> getDataFor(String famID_indID) {
        if (hashes.containsKey(famID_indID)) {
            return hashes.get(famID_indID);
        } else {
            return new Hashtable<String, CNVariant[]>();
        }
    }

    public void serialize(String filename) {
        SerializedFiles.writeSerial(this, filename);
    }

    public static CNVariantHash load(String filename, int structureType, boolean jar, Logger log) {
        CNVariantHash hashes = null;
        String suffix;

        if (structureType == CONSTRUCT_BY_IND) {
            suffix = ".ind.ser";
        } else if (structureType == CONSTRUCT_ALL) {
            suffix = ".all.ser";
        } else {
            System.err.println("Error - invalid CONSTRUCT type");
            suffix = null;
        }

        boolean parse = Files.exists(filename + suffix, jar);
        
        if (parse) {
            hashes = (CNVariantHash) SerializedFiles.readSerial(filename + suffix, jar, log, false);
        } 
        if (!parse || hashes == null) {
        	if (hashes == null) {
        		log.report("Detected that CNVariantHash needs to be updated from cnv.var.CNVariantHash to filesys.CNVariantHash; reparsing...");
        	}        		
            hashes = new CNVariantHash(filename, structureType, jar, log);
            hashes.serialize(filename + suffix);
        }

        hashes.setFilename(filename);

        return hashes;
    }

    /**
     * This method will return all of the CNVs within a specified range
     * 
     * @param chr
     *            Which chromosome to use
     * @param start
     *            Starting position
     * @param stop
     *            Ending position
     * @param minProbes
     *            Minimum number of probes (markers)
     * @param minSizeKb
     *            Minimum size in kilobases
     * @param minQualityScore
     *            Minimum quality score
     * @return Sorted array of CNVs in the region
     */
    public CNVariant[] getAllInRegion(byte chr, int start, int stop, int minProbes, int minSizeKb, int minQualityScore) {
        Vector<CNVariant> inRegion = new Vector<CNVariant>();

        Enumeration<String> e = hashes.keys();
        while (e.hasMoreElements()) {
            Hashtable<String, CNVariant[]> h = hashes.get(e.nextElement());
            Enumeration<String> ee = h.keys();
            while (ee.hasMoreElements()) {
                CNVariant[] cnv = h.get(ee.nextElement());
                for (int i = 0; i < cnv.length; i++) {
                    byte myChr = cnv[i].getChr();

                    // Only look at the specified chromosome
                    if (myChr == chr) {
                        int myStart = cnv[i].getStart();
                        int myStop = cnv[i].getStop();

                        if ((myStop < start) || (myStart > stop)) {
                            // We stop before the window, or start after it
                        } else {
                            int minSizeBases = minSizeKb * 1000;
                            if ((cnv[i].getNumMarkers() >= minProbes) && (cnv[i].getSize() >= minSizeBases) && (cnv[i].getScore() >= minQualityScore)) {
                                inRegion.add(cnv[i]);
                            }
                        }
                    }
                }
            }
        }
        return CNVariant.sortCNVs(CNVariant.toCNVariantArray(inRegion));
    }

    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\data\\penncnv_1SNP.cnv";

        String usage = "\n" + "cnv.var.CNVariantHash requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("file=")) {
                filename = args[i].split("=")[1];
                numArgs--;
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            new CNVariantHash(filename, 1, false, new Logger());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
