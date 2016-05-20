package cnv.plots;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.TreeMap;

import javax.swing.SwingUtilities;

import cnv.filesys.Project;
import common.Aliases;
import common.Files;
import common.Logger;
import common.ext;

class MetaStudy {
	private final ArrayList<StudyData> studies;
	private ArrayList<StudyData> sorted;
	private final HashMap<String, StudyData> nameMap;
	private final float metaBeta;
	private final float metaStderr;
	private final float[] metaConf = new float[2];
    private ArrayList<String> sortOrder = null;
    private boolean shouldSort;
    private boolean currentSortIsNaturalSort = false;
	
	public MetaStudy(float metaBeta, float metaStderr) {
		studies = new ArrayList<StudyData>();
		nameMap = new HashMap<String, StudyData>();
		this.metaBeta = metaBeta;
		this.metaStderr = metaStderr;
		this.metaConf[0] = (float) (metaBeta - 1.96 * metaStderr);
		this.metaConf[1] = (float) (metaBeta + 1.96 * metaStderr);
	}
	
	public void addStudy(StudyData studyData) {
		studies.add(studyData);
		nameMap.put(studyData.getLabel(), studyData);
	}

	String findLongestStudyName() {
		String longest = "";
		for(StudyData ft : getStudies()){
			longest = longest.length() < ft.getDisplayLabel().length() ? ft.getDisplayLabel() : longest;
		}
		return longest;
	}
	
	float calcSumZScore() {
		float sum = 0;
		for	(StudyData ft : studies){
			sum += ft.getZScore();
		}
		return sum;
	}
	
	float findMaxZScore() {
		float max = Float.MIN_VALUE;
		for (StudyData data: studies) {
			max = Math.max(max, data.getZScore());
		}
		return max;
	}
	
	private ArrayList<StudyData> getSorted(ArrayList<String> order) {
	    if (this.sorted == null || this.sorted.isEmpty() || currentSortIsNaturalSort) {
	        this.sorted = new ArrayList<StudyData>();
	        for (int i = order.size() - 1; i >= 0; i--) {
	            String name = ext.replaceAllWith(order.get(i), ForestPlot.REPLACEMENTS_FOOLISHLY_HARD_CODED);
	            String repl = null;
	            if (name.equals("")) {
	                this.sorted.add(new StudyBreak());
	            } else {
	                if (name.split("\t").length > 1) {
	                    repl = name.split("\t")[1];
	                    name = name.split("\t")[0];
	                }
	                StudyData sd = nameMap.get(name);
	                if (sd == null) {
	                    sd = new StudyBreak();
	                } else {
    	                if (repl != null) {
    	                    sd.setReplacementLabel(repl);
    	                }
	                }
	                sorted.add(sd);
	            }
	        }
	    }
        currentSortIsNaturalSort = false;
	    return this.sorted;
	}
	
	private ArrayList<StudyData> getSorted() {
		if (this.sorted == null || !currentSortIsNaturalSort) {
			this.sorted = new ArrayList<StudyData>();
			
			TreeMap<String, String> zeroStudyMap = new TreeMap<String, String>();
			TreeMap<Float, String> betaStudyMap = new TreeMap<Float, String>();
			for (StudyData study : studies) {
				if (study.getBeta(false) == 0.0f) {
					zeroStudyMap.put(study.getLabel(), study.getLabel());
				} else {
					betaStudyMap.put(study.getBeta(false), study.getLabel());
				}
			}
			ArrayList<StudyData> desc = new ArrayList<StudyData>();
			for (java.util.Map.Entry<String, String> entry : zeroStudyMap.entrySet()) {
				desc.add(nameMap.get(entry.getValue()));
			}
			for (java.util.Map.Entry<Float, String> entry : betaStudyMap.entrySet()) {
				desc.add(nameMap.get(entry.getValue()));
			}
			for (int i = desc.size() - 1; i >= 0; i--) {
				sorted.add(desc.get(i));
			}
		}
		currentSortIsNaturalSort = true;
		return this.sorted;
	}

	public ArrayList<StudyData> getStudies() {
	    return this.shouldSort ? (this.sortOrder == null || this.sortOrder.isEmpty() ? getSorted() : getSorted(this.sortOrder)) : this.studies;
	}

	public float[] getMetaConf(boolean odds) {
		return odds ? new float[]{(float) Math.exp(metaConf[0]), (float) Math.exp(metaConf[1])} : metaConf;
	}

//	public float getMetaBeta() {
//		return metaBeta;
//	}
//
//	public float getMetaStderr() {
//	    return metaStderr;
//	}
	public float getMetaBeta(boolean odds) {
	    return (float) (odds ? Math.exp(metaBeta) : metaBeta);
	}
	
	public float getMetaStderr(boolean odds) {
	    return (float) (odds ? Math.exp(metaStderr) : metaStderr);
	}

    public void setSort(boolean sortedDisplay, ArrayList<String> sortOrder) {
        this.shouldSort = sortedDisplay;
        this.sortOrder = sortOrder;
        this.sorted = null;
    }
	
}

class StudyBreak extends StudyData {
    // placeholder class for visual breaks
    public StudyBreak() {
        this ("", 0f, 0f, 0, (byte) 0);
    }
    private StudyBreak(String label, float beta, float stderr, int color, byte shape) {
        super(label, beta, stderr, color, shape);
    }
}

class StudyData {
	private final String label;
	private String replLabel = null;
	private final float beta;
	private final float stderr;
	private int color;
	private byte shape;
	private final float[] confInterval;
	private final float zScore;

	public StudyData(String label, float beta, float stderr, int color, byte shape) {
		this.label = label;
		this.beta = beta;
		this.stderr = stderr;
		this.color = color;
		this.shape = shape;
		this.confInterval = new float[2];
		this.confInterval[0] = (float) (beta - 1.96 * stderr);
		this.confInterval[1] = (float) (beta + 1.96 * stderr);
		this.zScore = stderr == 0.0f ? 0.0f : Math.abs(beta / stderr);
	}

	public void setReplacementLabel(String repl) {
       this.replLabel = repl;
    }
	
	public String getDisplayLabel() {
	    if (replLabel != null) {
	        return replLabel;
	    }
	    return label;
	}
	
    public String getLabel() {
		return label;
	}

	public float getBeta(boolean odds) {
		return (float) (odds ? Math.exp(beta) : beta);
	}

	public float getStderr(boolean odds) {
		return (float) (odds ? Math.exp(stderr) : stderr);
	}

	public int getColor() {
		return color;
	}

	public byte getShape() {
		return shape;
	}

	public float[] getConfInterval(boolean odds) {
		return odds ? new float[]{(float) Math.exp(confInterval[0]), (float) Math.exp(confInterval[1])} : confInterval;
	}

	public float getZScore() {
		return zScore;
	}
}

class ForestInput {
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((comment == null) ? 0 : comment.hashCode());
		result = prime * result + ((file == null) ? 0 : file.hashCode());
		result = prime * result + ((marker == null) ? 0 : marker.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ForestInput other = (ForestInput) obj;
		if (comment == null) {
			if (other.comment != null)
				return false;
		} else if (!comment.equals(other.comment))
			return false;
		if (file == null) {
			if (other.file != null)
				return false;
		} else if (!file.equals(other.file))
			return false;
		if (marker == null) {
			if (other.marker != null)
				return false;
		} else if (!marker.equals(other.marker))
			return false;
		return true;
	}
	
	final String marker;
	final String file;
	final String comment;
	int[] metaIndicies;
	HashMap<String, Integer> studyToColIndexMap;
	ArrayList<String> studyList;
	MetaStudy ms;
	
	public ForestInput(String marker, String file, String comment) {
		this.marker = marker;
		this.file = file;
		this.comment = comment;
		this.studyToColIndexMap = new HashMap<String, Integer>();
		this.studyList = new ArrayList<String>();
	}

	public void addStudy(String string, int i) {
		studyList.add(string);
		studyToColIndexMap.put(string, i);
	}
	public MetaStudy getMetaStudy() { return ms; }
	public void setMetaStudy(MetaStudy ms) { this.ms = ms; }
}

public class ForestPlot {
	public static final  String[][] REPLACEMENTS_FOOLISHLY_HARD_CODED = new String[][] {{"_WBC_TOTAL", ""},
																						{"_WBC_NEUTRO", ""},
																						{"_", " "},
	};
	
	private static final String[] BETA_META_HEADERS = {"beta", "effect"};
	private static final String[] SE_META_HEADERS = {"se", "stderr"};
	private static final String BETA_PREFIX = "beta.";
	private static final String SE_PREFIX = "se.";
	
	
	private Project proj;
	private String markerFileName;
	private Logger log;
	
	private ForestPanel forestPanel;
	
	private volatile boolean loadingFile;
	private Thread loadingThread;
	private boolean atleastOneStudy;
	private int currentDataIndex;
	private MetaStudy currMetaStudy;
	private float maxZScore;
	private float sumZScore;
	private String longestStudyName;
	private String plotLabel;
	
	private boolean sortedDisplay;
	private ArrayList<String> sortOrder;

	private ArrayList<ForestInput> dataIndices = new ArrayList<ForestInput>();


	/**
	 * @param proj
	 */
	public ForestPlot(Project proj) {
		super();
		this.proj = proj;
		this.log = proj.getLog();
		setup();
	}
	/**
	 * @param markerFileName
	 * @param log
	 */
	public ForestPlot(String markerFileName, Logger log) {
		super();
		this.markerFileName = markerFileName;
		this.log = log;
		setup();
	}
	
	private void setup() {
		forestPanel = new ForestPanel(this, log);
		forestPanel.setLayout(new BorderLayout());
	}
	
	private LinkedHashSet<ForestInput> readMarkerFile(String markerFile) {
	    String file;
		LinkedHashSet<ForestInput> markerNames = new LinkedHashSet<ForestInput>();
		BufferedReader markerReader = Files.getReader(markerFile, false, true, false);
		
		if (markerReader != null) {
    		try {
    			while (markerReader.ready() && !Thread.interrupted()) {
    				String[] line = markerReader.readLine().trim().split("\\t");
				    if (line.length >= 2) {
				        file = line[1];
				        if (!file.contains(":") && !file.startsWith("/") && !Files.exists(file)) {
				            if (Files.exists(ext.verifyDirFormat(ext.parseDirectoryOfFile(markerFile)) + file)) {
				                file = ext.verifyDirFormat(ext.parseDirectoryOfFile(markerFile)) + file;
				            } else {
				                if (log != null) {
				                    log.reportError("Error - file " + file + " not found!");
				                } else {
				                    System.err.println("Error - file " + file + " not found!");
				                }
				            }
				        }
    					markerNames.add(new ForestInput(line[0], file, line.length > 2 ? line[2] : ""));
    				} else if (line.length == 1) {
    					markerNames.add(new ForestInput(line[0], "", ""));
    				}
    			}
    		} catch (IOException e) {
    		    if (log != null) {
    		        log.reportException(e);
    		    } else {
    		        e.printStackTrace();
    		    }
    		}
		}
		
		return markerReader == null || Thread.interrupted() ? new LinkedHashSet<ForestInput>() : markerNames;
	}
	
	private void loadStudyData() throws InterruptedException, RuntimeException {
		HashMap<String, ArrayList<ForestInput>> files = new HashMap<String, ArrayList<ForestInput>>();
		for (ForestInput fi : dataIndices) {
			ArrayList<ForestInput> inputList = files.get(fi.file);
			if (inputList == null) {
				inputList = new ArrayList<ForestInput>();
				files.put(fi.file, inputList);
			}
			inputList.add(fi);
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
		}
		final HashMap<String, Integer> progSteps = new HashMap<String, Integer>();
		for (String file : files.keySet()) {
			int sz = Files.getSize(file, false);
			progSteps.put(file, sz);
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
		}
		
		
		for (final java.util.Map.Entry<String, ArrayList<ForestInput>> fileMap : files.entrySet()) {
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			try {
				LineNumberReader lnr = new LineNumberReader(new java.io.FileReader(fileMap.getKey()));
				lnr.skip(Long.MAX_VALUE);
				lnr.close();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
			}
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			
			BufferedReader dataReader;
			dataReader = Files.getReader(fileMap.getKey(), 
											false, // not a jar file
											true, // verbose mode on 
											log, 
											false // don't kill the whole process, esp. if we're running a GUI
											);
			if (dataReader == null) {
				continue;
			}
			if (Thread.interrupted()) {
				interruptLoading();
				return;
			}
			String delimiter = Files.determineDelimiter(fileMap.getKey(), log);
			String header;
			try {
				header = dataReader.readLine();	// skip header
				String[] hdr;
				if (delimiter.startsWith(",")) {
					hdr = ext.splitCommasIntelligently(header, true, log);
				} else {
					hdr = header.trim().split(delimiter);
				}
				int idIndex = ext.indexFactors(new String[][]{Aliases.MARKER_NAMES}, hdr, false, true, false, false)[0];
				
				for (ForestInput inputData : fileMap.getValue()) {
					if (Thread.interrupted()) {
						interruptLoading();
						return;
					}
					try {
						mapMarkersToCol(inputData, header);
					} catch (InterruptedException e) {
						inputData.metaIndicies = null;
						inputData.studyList.clear();
						inputData.studyToColIndexMap.clear();
						interruptLoading();
						return;
					}
				}
				while (dataReader.ready() && !Thread.interrupted()) {
					String readLine = dataReader.readLine();
					String readData[] = delimiter.equals(",") ? ext.splitCommasIntelligently(readLine, true, log) : readLine.split(delimiter);
					String markerName = readData[idIndex];
					for (ForestInput inputData : fileMap.getValue()) {
						if (Thread.interrupted()) {
							dataReader.close();
							interruptLoading();
							return;
						}
						if(inputData.marker.equals(markerName)) {
//							dataToMetaMap.put(inputData, getMetaStudy(inputData, readData));
							getMetaStudy(inputData, readData);
							atleastOneStudy = true;
						}
					}
				}
				dataReader.close();
				
			} catch (IOException e) {
				log.reportException(e);
			}
	
			if(!atleastOneStudy) {
				log.reportError("Not able to find data for file '"+fileMap.getKey()+"'. Please make sure the given markers are correct and included in data file.");
			}
			

		}
		
		
	}
	
	private void mapMarkersToCol(ForestInput data, String hdr) throws RuntimeException, InterruptedException {
		String delim = data.file.toLowerCase().endsWith(".csv") ? ",!" : ext.determineDelimiter(hdr);
		String[] dataFileHeaders = delim.startsWith(",") ? ext.splitCommasIntelligently(hdr, delim.endsWith("!"), log) : hdr.trim().split(delim);
		for (int i = 0; i < dataFileHeaders.length; i++) {
			for (int j = 0; j < BETA_META_HEADERS.length; j++) {
				if (dataFileHeaders[i].toLowerCase().equals(BETA_META_HEADERS[j])) {
					if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_META_HEADERS[j])) {
						data.metaIndicies = new int[]{i, i+1};
					}
				}
			}
			if (dataFileHeaders[i].toLowerCase().startsWith(BETA_PREFIX)) {
				if (dataFileHeaders[i + 1].toLowerCase().startsWith(SE_PREFIX)) {
					if (data.studyToColIndexMap.containsKey(dataFileHeaders[i].split("\\.")[1])) {
						throw new RuntimeException("Malformed data file: Duplicate study name found in file");
					} else {
						data.addStudy(dataFileHeaders[i].split("\\.")[1], i);
					}
				} else {
					throw new RuntimeException("Malformed data file: SE is not present after Beta for: " + dataFileHeaders[i]);
				}
			}
		}
		if (data.metaIndicies == null) {
			log.reportError("Error - no overall beta/se pairing or effect/stderr pairing was found in file "+data.file);
		}
		dataFileHeaders = null;
	}

	private void interruptLoading() throws InterruptedException {
		atleastOneStudy = false;
		markerFileName = "";
		dataIndices.clear();
//		this.dataToMetaMap.clear();
		clearCurrentData();
		throw new InterruptedException();
	}
	
	protected void setCurrentData(int index) {
		if (this.dataIndices.size() == 0 || index < 0 || index > this.dataIndices.size()) return;
		setCurrentDataIndex(index);
//		setCurrentMetaStudy(dataToMetaMap.get(dataIndices.get(index)));
		setCurrentMetaStudy(dataIndices.get(index).getMetaStudy());
		getCurrentMetaStudy().setSort(isSortedDisplay(), getSortOrder());
		if (getCurrentMetaStudy() == null) {
		    String msg = "Error - could not set index to "+index+" since the data did not load properly; check to see if any results files are missing";
		    if (log != null) {
		        log.reportError(msg);
		    } else {
		        System.err.println(msg);
		    }
			return;
		}
		maxZScore = getCurrentMetaStudy().findMaxZScore();
		sumZScore = getCurrentMetaStudy().calcSumZScore();
		longestStudyName = getCurrentMetaStudy().findLongestStudyName();
		setPlotLabel(dataIndices.get(index).marker);
	}
	
	private void clearCurrentData() {
		setCurrentDataIndex(-1);
		setCurrentMetaStudy(null);
		maxZScore = 0;
		maxZScore = 0;
		sumZScore = 0;
		longestStudyName = "";
		setPlotLabel("");
	}
	
	private void getMetaStudy(ForestInput input, String[] readData) {
		String metaB, metaS;
		float metaBeta, metaStderr;
		
		metaB = readData[input.metaIndicies[0]];
		metaS = readData[input.metaIndicies[1]];
		metaBeta = ext.isValidDouble(metaB) ? Float.parseFloat(metaB) : 0.0f;
		metaStderr = ext.isValidDouble(metaS) ? Float.parseFloat(metaS) : 0.0f;
		
		MetaStudy ms = new MetaStudy(metaBeta, metaStderr);
		ArrayList<StudyData> studies = getStudyEntries(input, readData);
		for (int i = 0; i < studies.size(); i++) {
			ms.addStudy(studies.get(i));
		}
		
		input.setMetaStudy(ms);
//		return ms;
	}
	
	private ArrayList<StudyData> getStudyEntries(ForestInput input, String[] readData) {
		ArrayList<StudyData> studies = new ArrayList<StudyData>();
		String betaVal, seVal;
		for (int i = input.studyList.size() - 1; i >= 0; i--) {
			String studyName = input.studyList.get(i);
			betaVal = readData[input.studyToColIndexMap.get(studyName)];
			seVal = readData[input.studyToColIndexMap.get(studyName) + 1];
			float beta = ext.isValidDouble(betaVal) ? Float.parseFloat(betaVal) : 0.0f;
			float stderr = ext.isValidDouble(seVal) ? Float.parseFloat(seVal) : 0.0f;
			studies.add(new StudyData(ext.replaceAllWith(studyName, REPLACEMENTS_FOOLISHLY_HARD_CODED), beta, stderr, 0, PlotPoint.FILLED_CIRCLE));
		}
		return studies;
	}
	
	
	public void loadMarkerFile() {
		loadMarkerFile(new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						reloadData();
					} catch (InterruptedException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					forestPanel.setPointsGeneratable(markerFileName != null);
					forestPanel.setRectangleGeneratable(markerFileName != null);
					forestPanel.setExtraLayersVisible(new byte[] { 99 });
					
					try {
						SwingUtilities.invokeAndWait(new Runnable() {
							@Override
							public void run() {
								updateGUI();
							}
						});
					} catch (InvocationTargetException e) {
						// TODO Auto-generated catch block
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
					}
					loadingFile = false;
				}
			}, "ForestPlot_loadMarkerFile"));
	}
	
	protected void loadMarkerFile(Thread loadingThread) {
		if (!this.loadingFile) {
			this.loadingFile = true;
			this.loadingThread = loadingThread;
			this.loadingThread.start();
		}
	}
	
	protected void reloadData() throws InterruptedException {
		atleastOneStudy = false;
		
		dataIndices = new ArrayList<ForestInput>();
		if (markerFileName != null) {
			dataIndices.addAll(readMarkerFile(markerFileName));
		} 
		
		if (!dataIndices.isEmpty() && !Thread.interrupted()) {
			loadStudyData();
			setCurrentData(0);
		} else {
			clearCurrentData();
		}
		
	}
	
    public boolean waitForLoad() {
        while(isLoadingFile()) {
            Thread.yield();
        }
        return true;
    }
    
    /**
	 * NOTE: be sure to call 'waitForLoad' before calling screenCap, otherwise data might not be loaded in time.
	 * @param subdir
	 * @param odds
	 */
	public void screenCapAll(String subdir, boolean odds, boolean versionIfExists) {
	    setOddsRatioDisplay(odds);
        ArrayList<ForestInput> data = getDataIndices();
//        for (int i = 0; i < data.size(); i++) {
        for (int i = 0; i < 2; i++) {
            setCurrentData(i);
            getForestPanel().addNotify();
            getForestPanel().setSize(new Dimension(1000,720));
            getForestPanel().createImage();
            getForestPanel().validate();
            String marker = "", filename = "", dataFile = "";
            int count = 1;
            String root = (proj == null ? ext.parseDirectoryOfFile(getMarkerFileName()) : proj.PROJECT_DIRECTORY.getValue());
            root = ext.verifyDirFormat(root);
            if (subdir != null && !subdir.equals("")) {
                root += subdir;
                root = ext.verifyDirFormat(root);
            }
            marker = getDataIndices().get(getCurrentDataIndex()).marker;
            dataFile = ext.rootOf(getDataIndices().get(getCurrentDataIndex()).file, true);
            filename = marker + "_" + dataFile;
            filename = ext.replaceWithLinuxSafeCharacters(filename, true);
            if (new File(root + filename + ".png").exists()) {
                if (versionIfExists) {
                    while (new File(root+filename+".png").exists()) {
                        filename = marker + "_" + dataFile + "_v" + count;
                        filename = ext.replaceWithLinuxSafeCharacters(filename, true);
                        count++;
                    }
                }
            }
            if (getLog() != null) {
            	getLog().report("Writing screenshot to file " + root + filename + ".png");
            } else {
                System.out.println("Writing screenshot to file " + root + filename + ".png");
            }
            getForestPanel().screenCapture(root+filename+".png");
        }
	}
	
	public void updateGUI() {
		forestPanel.paintAgain();
	}
	
	public void loadOrderFile(String filename, boolean shouldSort) {
	    if (!Files.exists(filename)) {
	        String msg = "Error - study order file (" + filename + ") not found!";
	        if (log != null) {
	            log.reportError(msg);
	        } else {
	            System.err.println(msg);
	        }
	        return;
	    }
	    ArrayList<String> order = new ArrayList<String>();
	    try {
    	    BufferedReader reader = Files.getAppropriateReader(filename);
    	    String line = null;
            while((line = reader.readLine()) != null) {
                order.add(line.trim());
            }
            reader.close();
        } catch (IOException e) {
            if (proj != null) {
                proj.message("Error occurred while loading study order file: " + e.getMessage());
            }
            if (log != null) {
                log.reportException(e);
            } else {
                e.printStackTrace();
            }
        }
	    if (log != null) {
	        log.report("Loaded Study Order File: " + filename);
	    } else {
	        System.out.println("Loaded Study Order File: " + filename);
	    }
	    setSortOrder(order);
	    setSortedDisplay(shouldSort);
        if (currMetaStudy != null) {
            currMetaStudy.setSort(this.isSortedDisplay(), this.getSortOrder());
        }
	}
	
	public static void generateNaiveOrderFile(String[] studies, String outFile) {
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		for (String study : studies) {
			String safeStudy = ext.replaceWithLinuxSafeCharacters(study, true);
			if (study.equals(safeStudy)) writer.println(safeStudy);
			else writer.println(safeStudy + "\t" + study);
		}
		writer.flush();
		writer.close();
	}
	
	public ForestPanel getForestPanel() {
		return forestPanel;
	}
	
	public String getMarkerFileName() {
		return markerFileName;
	}
	
	public void setMarkerFileName(String file) {
		markerFileName = file;
		
	}
	
	public int getCurrentDataIndex() {
		return currentDataIndex;
	}

	private void setCurrentDataIndex(int currentDataIndex) {
		this.currentDataIndex = currentDataIndex;
	}
	
	public MetaStudy getCurrentMetaStudy() {
	    if (currMetaStudy != null) {
	        currMetaStudy.setSort(this.isSortedDisplay(), this.getSortOrder());
	    }
		return currMetaStudy;
	}
	
	private void setCurrentMetaStudy(MetaStudy currMetaStudy) {
		this.currMetaStudy = currMetaStudy;
	}
	

    public boolean isSortedDisplay() {
        return sortedDisplay;
    }

    public ArrayList<String> getSortOrder() {
        return sortOrder;
    }

    public void setSortOrder(ArrayList<String> sortOrder) {
        this.sortOrder = sortOrder;
    }
    
	public float getMaxZScore() {
		return maxZScore;
	}

	public float getSumZScore() {
		return sumZScore;
	}
	
	public String getPlotLabel() {
		return plotLabel;
	}

	private void setPlotLabel(String plotLabel) {
		this.plotLabel = plotLabel;
	}
	
	public ArrayList<ForestInput> getDataIndices() {
		return dataIndices;
	}
	
	public String getLongestStudyName() {
		return longestStudyName;
	}
	public Logger getLog() {
		return log;
	}

    public void setSortedDisplay(boolean sorted) {
        this.sortedDisplay = sorted;
    }
    
	public void setOddsRatioDisplay(boolean selected) {
        forestPanel.oddsDisplay = selected;
    }
	
	protected boolean isLoadingFile() {
		return loadingFile;
	}
	
	protected Thread getLoadingThread() {
		return loadingThread;
	}
	
	protected void setLoadingFile(boolean loadingFile) {
		this.loadingFile = loadingFile;
	}

	
}
