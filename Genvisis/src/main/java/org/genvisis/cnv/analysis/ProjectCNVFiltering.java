package cnv.analysis;

import java.util.Hashtable;

import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.SampleData;

import common.Array;
import common.CNVFilter;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import common.CNVFilter.FreqFilter;

import filesys.CNVariant;
import filesys.LocusSet;

public class ProjectCNVFiltering {
    

    public static void setIndividualsToKeepFromSampleData(CNVFilter filter, Project proj) {
        String samples[] = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
        SampleData sampleData = proj.getSampleData(0, false);
        for (int i = 0; i < samples.length; i++) {
            samples[i] = sampleData.lookup(samples[i])[1];
        }
        proj.getLog().report("Info - retaining " + samples.length + (samples.length == 1 ? " sample" : " samples") + " in the cnv filter set");
        filter.setIndHash(samples);
    }


    public static void setCentromereBoundariesFromProject(CNVFilter filter, Project proj) {
//      setCentromereBoundariesFromFile(proj.getFilename(proj.MARKERSET_FILENAME));
        filter.setCentromereBoundariesFromFile(proj.MARKERSET_FILENAME.getValue());
    }
    
    public static void setCNVDefaults(CNVFilter filter, Project proj) {
        filter.setMinNumMarkers(CNVFilter.DEFAULT_MIN_NUM_MARKERS);
        filter.setMinScore(CNVFilter.DEFAULT_MIN_SCORE_THRESHOLD);
        setCentromereBoundariesFromProject(filter, proj);
        filter.setBreakupCentromeres(CNVFilter.DEFAULT_BREAK_UP_CENTROMERES);
        setIndividualsToKeepFromSampleData(filter, proj);
        filter.setCommonIn(CNVFilter.DEFAULT_COMMON_IN);
    }

    private static void loadMarkerSet(Project proj, CNVFilter filter) {
    	String markerSetFile = proj.MARKERSET_FILENAME.getValue();
    	MarkerSet markerSet = MarkerSet.load(markerSetFile, false);
    	filter.setPositions(markerSet.getPositionsByChr());
        filter.setCentromereBoundaries(Positions.determineCentromereBoundariesFromMarkerSet(markerSet.getChrs(), markerSet.getPositions(), proj.GENOME_BUILD_VERSION.getValue().getBuildInt(), proj.getLog()));
        filter.computeCentromereMidPoints();
    }
    
    /**
     * @param proj
     * @param args
     *            these are usually the args passed to a main method, if an arg is not valid we do not warn or exit, simply skip
     * @param filter
     *            a pre-existing filter, can be null
     * @param defaults
     *            set the filter to the default values as defined by {@link CNVFilter#setCNVDefaults(Project)}
     * @param log
     * @return
     */
    public static CNVFilter setupCNVFilterFromArgs(Project proj, String[] args, CNVFilter filter, boolean defaults, Logger log) {
        if (filter == null) {
            filter = new CNVFilter(log);

        }
        if (defaults) {
            setCNVDefaults(filter, proj);
        }
        filter.setCommandLineFiltersInEffect(new Hashtable<String, String>());
//      filter.setCentromereBoundariesFromFile(proj.getFilename(proj.MARKERSET_FILENAME));// resets if neccesary
        loadMarkerSet(proj, filter);
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith(CNVFilter.COMMAND_MIN_NUM_MARKERS)) {
                filter.setMinNumMarkers(ext.parseIntArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MIN_NUM_MARKERS);
            } else if (args[i].startsWith(CNVFilter.COMMAND_MAX_NUM_MARKERS)) {
                filter.setMaxNumMarkers(ext.parseIntArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MAX_NUM_MARKERS);
            } else if (args[i].startsWith(CNVFilter.COMMAND_MIN_SIZE)) {
                filter.setMinSize(ext.parseIntArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MIN_SIZE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_MAX_SIZE)) {
                filter.setMaxSize(ext.parseIntArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MAX_SIZE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_MIN_SCORE)) {
                filter.setMinScore(ext.parseDoubleArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MIN_SCORE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_MAX_SCORE)) {
                filter.setMaxScore(ext.parseDoubleArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_MAX_SCORE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_PROBLEM_REGIONS_FILE)) {
                filter.setProblemRegionsFromFile(proj.PROJECT_DIRECTORY.getValue() + ext.parseStringArg(args[i], ""));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_PROBLEM_REGIONS_FILE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_COMMON_REFERENCE)) {
                filter.setCommonReferenceFromFile(proj.PROJECT_DIRECTORY.getValue() + ext.parseStringArg(args[i], ""), CNVFilter.DEFAULT_COMMON_IN);
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_COMMON_REFERENCE);
            } else if (args[i].startsWith(CNVFilter.COMMAND_INDIVIDUALS_TO_KEEP)) {
                filter.setIndividualsToKeepFromFile(proj.PROJECT_DIRECTORY.getValue() + ext.parseStringArg(args[i], ""));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_INDIVIDUALS_TO_KEEP);
            } else if (args[i].startsWith(CNVFilter.COMMAND_BUILD)) {
                filter.setBuild(ext.parseIntArg(args[i]));
//              filter.setCentromereBoundariesFromFile(proj.getFilename(proj.MARKERSET_FILENAME));
                filter.setCentromereBoundariesFromFile(proj.MARKERSET_FILENAME.getValue());
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_BUILD);
            } else if (args[i].startsWith(CNVFilter.COMMAND_COMMON_IN)) {
                filter.setCommonIn(ext.parseBooleanArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_COMMON_IN);
            } else if (args[i].startsWith(CNVFilter.COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA)) {
                if (filter.getIndHash() == CNVFilter.NO_FILTER_INDIVIDUAL_HASH && ext.parseBooleanArg(args[i])) {
                    setIndividualsToKeepFromSampleData(filter, proj);
                } else if (!ext.parseBooleanArg(args[i])) {
                    proj.getLog().reportTimeWarning("Over-riding sample data exclusion");
                    filter.setIndHash(CNVFilter.NO_FILTER_INDIVIDUAL_HASH);
                }
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_EXCLUDE_INDIVIDUALS_FROM_SAMPLE_DATA);
            } else if (args[i].startsWith(CNVFilter.COMMAND_BREAK_UP_CENTROMERES)) {
                filter.setBreakupCentromeres(ext.parseBooleanArg(args[i]));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_BREAK_UP_CENTROMERES);
            } else if (args[i].startsWith(CNVFilter.COMMAND_BREAK_UP_CENTROMERES_SOURCE_FILE)) {
                filter.setCentromereBoundariesFromFile(ext.parseStringArg(args[i], null));
                filter.addCommandLineFilter(args[i], CNVFilter.COMMAND_BREAK_UP_CENTROMERES_SOURCE_FILE);
            }
        }
        return filter;

    }

    public static LocusSet<CNVariant> filterCNVFile(Project proj, String cnvFile, String out, CNVFilter cnvFilter, boolean mergePrior, boolean freqFilterPrior, FreqFilter freqFilter, boolean ifMergeAndFreqMergeFirst) {
        CNVariant[] cnvs = CNVariant.loadPlinkFile(proj.PROJECT_DIRECTORY.getValue() + cnvFile, false);
        return filterCNVFile(proj, cnvs, out, cnvFilter, mergePrior, freqFilterPrior, freqFilter, ifMergeAndFreqMergeFirst, false);
    }
    /**
     * 
     * @param cnvsAsPositions
     *            get the start and stop positions from the {@link CNVariant}s themselves so that a {@link MarkerSet} is not needed
     */
    public static LocusSet<CNVariant> filterCNVFile(Project proj, CNVariant[] cnvs, String out, CNVFilter cnvFilter, boolean mergePrior, boolean freqFilterPrior, FreqFilter freqFilter, boolean ifMergeAndFreqMergeFirst, boolean cnvsAsPositions) {
        Logger log = proj.getLog();
        if (mergePrior && freqFilterPrior) {
            if (ifMergeAndFreqMergeFirst) {
                int numPrior = cnvs.length;
                int iter = 1;
                do {
                    numPrior = cnvs.length;
                    log.report("Merging CNVs, iteration " + iter);
                    cnvs = FilterCalls.mergeCNVsInMemory(proj, cnvs, FilterCalls.DEFAULT_CLEAN_FACTOR, cnvsAsPositions);
                    log.report("CNV merging iteration " + iter++ + " + complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");
                } while (numPrior > cnvs.length);

                numPrior = cnvs.length;
                cnvs = FilterCalls.filterBasedOnNumberOfCNVsAtLocusInMemory(proj, cnvs, freqFilter.totalRequired, freqFilter.delRequired, freqFilter.dupRequired, freqFilter.totalLimitedTo, freqFilter.delLimitedTo, freqFilter.dupLimitedTo, freqFilter.proportionOfProbesThatNeedToPassForFinalInclusion, cnvsAsPositions);
                log.report("CNV filtering by frequency complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");
            } else {
                int numPrior = cnvs.length;
                cnvs = FilterCalls.filterBasedOnNumberOfCNVsAtLocusInMemory(proj, cnvs, freqFilter.totalRequired, freqFilter.delRequired, freqFilter.dupRequired, freqFilter.totalLimitedTo, freqFilter.delLimitedTo, freqFilter.dupLimitedTo, freqFilter.proportionOfProbesThatNeedToPassForFinalInclusion, cnvsAsPositions);
                log.report("CNV filtering by frequency complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");

                numPrior = cnvs.length;
                int iter = 1;
                do {
                    numPrior = cnvs.length;
                    log.report("Merging CNVs, iteration " + iter);
                    cnvs = FilterCalls.mergeCNVsInMemory(proj, cnvs, FilterCalls.DEFAULT_CLEAN_FACTOR, cnvsAsPositions);
                    log.report("CNV merging iteration " + iter++ + " + complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");
                } while (numPrior > cnvs.length);
            }
        } else if (mergePrior) {
            int numPrior = cnvs.length;
            int iter = 1;
            do {
                numPrior = cnvs.length;
                log.report("Merging CNVs, iteration " + iter);
                cnvs = FilterCalls.mergeCNVsInMemory(proj, cnvs, FilterCalls.DEFAULT_CLEAN_FACTOR, cnvsAsPositions);
                log.report("CNV merging iteration " + iter++ + " + complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");
            } while (numPrior > cnvs.length);
        } else if (freqFilterPrior) {
            int numPrior = cnvs.length;
            cnvs = FilterCalls.filterBasedOnNumberOfCNVsAtLocusInMemory(proj, cnvs, freqFilter.totalRequired, freqFilter.delRequired, freqFilter.dupRequired, freqFilter.totalLimitedTo, freqFilter.delLimitedTo, freqFilter.dupLimitedTo, freqFilter.proportionOfProbesThatNeedToPassForFinalInclusion, cnvsAsPositions);
            log.report("CNV filtering by frequency complete: started with " + numPrior + " CNVs, now have " + cnvs.length + " CNVs remaining.");
        }
        
        String outFile = Files.isRelativePath(out) ? proj.PROJECT_DIRECTORY.getValue() + out : out;
        return CNVFilter.filterCNVs(cnvs, outFile, cnvFilter, log);
    }
    
    public static void main(String[] args) {
        int numArgs = args.length;
        String filename = null;
        String cnvFile = "Genvisis.cnv";
        String out = "Genvisis.filt.cnv";
        boolean merge = false;
        boolean freq = false;
        boolean mergeFirst = false;
        String logfile = null;
        Project proj;

        String usage = "\n" + "cnv.qc.CNVFilter requires 1 argument\n";
        usage += "   (1) project file name (i.e. " + CNVFilter.COMMAND_PROJECT + filename + " (no default))\n" + "";
        usage += "   (2) cnv file name (relative to the project directory) (i.e. " + CNVFilter.COMMAND_CNV_FILE + cnvFile + " ( default))\n" + "";
        usage += "   (3) output (relative to the project directory) (i.e. " + CNVFilter.COMMAND_CNV_FILE_OUT + out + " ( default))\n" + "";

        if (ext.indexOfStr(CNVFilter.COMMAND_PROJECT, args, true, false) >= 0) {
            proj = new Project(ext.parseStringArg(args[ext.indexOfStr(CNVFilter.COMMAND_PROJECT, args, true, false)], ""), logfile, false);
        } else {
            proj = new Project(filename, null, false);
        }
        CNVFilter cnvFilter = setupCNVFilterFromArgs(proj, args, null, true, proj.getLog());
        FreqFilter freqFilter = CNVFilter.setupFreqFilterFromArgs(args, cnvFilter, proj.getLog());
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
                System.err.println(usage);
                System.exit(1);
            } else if (args[i].startsWith("file=")) {
                filename = ext.parseStringArg(args[i], "");
                numArgs--;
            } else if (args[i].startsWith("log=")) {
                logfile = ext.parseStringArg(args[i], "");
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_CNV_FILE_OUT)) {
                out = ext.parseStringArg(args[i], "");
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_CNV_FILE)) {
                cnvFile = ext.parseStringArg(args[i], "");
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_PROJECT)) {
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_MERGE)) {
                merge = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_FREQ_FILTER)) {
                freq = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (args[i].startsWith(CNVFilter.COMMAND_MERGE_FIRST)) {
                mergeFirst = ext.parseBooleanArg(args[i]);
                numArgs--;
            } else if (cnvFilter.isCommandLineFilterInEffect(args[i])) {
                numArgs--;
            } else {
                System.err.println("Error - invalid argument: " + args[i]);
            }
        }
        if (numArgs != 0) {
            System.err.println(usage);
            System.exit(1);
        }
        try {
            filterCNVFile(proj, cnvFile, out, cnvFilter, merge, freq, freqFilter, mergeFirst);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
}
