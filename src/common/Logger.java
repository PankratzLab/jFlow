package common;

import java.io.*;

import javax.swing.JTextArea;

public class Logger implements Serializable {
	private static final long serialVersionUID = 1L;
	public static final int LEVEL_ONE = 1;
	
	private String filename;
	private boolean logging;
	private int level;
	private JTextArea textArea;
	
	public Logger() {
		this(null, false);
	}
	
	public Logger(String filename) {
		this(filename, false);
	}
	
	public Logger(String filename, boolean append) {
		this(filename, append, 10); 
	}
	
	public Logger(String filename, boolean append, int level) {
		this.filename = filename;
		logging = filename!=null;
		if (logging && !append) {
			new File(filename).delete();
		}
		this.level = level;
		this.textArea = null;
	}
	
	public void linkTextArea(JTextArea text) {
		textArea = text;
	}
	
	public String getFilename() {
		return filename;
	}

	public int getLevel() {
		return level;
	}

	public void setLevel(int level) {
		this.level = level;
	}

	/**
	 * @param str
	 *            report this string with a time stamp and info message
	 */
	public void reportTimeInfo(String str) {
		report(ext.getTime() + " Info - " + str, true, true);
	}

	/**
	 * @param str
	 *            report this string with a time stamp and info message
	 */
	public void reportTimeElapsed( long time) {
		reportTimeInfo("Time elapsed: " + ext.getTimeElapsed(time));
	}
	
	/**
	 * @param file
	 *            report this file already exists with a time stamp and info message
	 */
	public void reportFileExists(String file) {
		reportTimeInfo("File "+ file+" already exists...");
	}

	/**
	 * @param str
	 *            report this string with a time stamp and warning message
	 */
	public void reportTimeWarning(String str) {
		report(ext.getTime() + " Warning - " + str, true, true);
	}

	/**
	 * @param str
	 *            report this string with a time stamp and error message
	 */
	public void reportTimeError(String str) {
		reportError(ext.getTime() + " Error - " + str, true, true);
	}
	
	/**
	 * @param filename report that this file was not found with a time stamp and error message
	 */
	public void reportFileNotFound(String filename) {
		String str = "file \"" + filename + "\" not found in current directory";
		reportTimeError(str);
	}

	/**
	 * @param filenames
	 *            report the files that do not exist in this string array
	 * 
	 */
	public void reportFilesNotFound(String[] filenames) {
		if (filenames == null) {
			reportTimeError("No file handles provided");
		} else {
			boolean haveMissing =false;
			for (int i = 0; i < filenames.length; i++) {
				if (!Files.exists(filenames[i])) {
					reportFileNotFound(filenames[i]);
					haveMissing =true;
				}
			}
			if (!haveMissing) {
				reportTimeInfo("All " + filenames.length + " files exist");
			}
		}
	}

	/**
	 * @param filename
	 *            report that this file had an IO err with time stamp and error message
	 */
	public void reportIOException(String filename) {
		String str = "could not read file \"" + filename + "\"";
		reportTimeError(str);
	}
	
	public void report(String str) {
		report(str, true, true);
	}
	
	public void report(String str, boolean line, boolean reportToScreen) {
		report(str, line, reportToScreen, 0);
	}
	
	public void report(String str, boolean line, boolean reportToScreen, int levelRequiredToReport) {
		PrintWriter writer;
	
		if (level >= levelRequiredToReport && reportToScreen) {
			if (line) {
				System.out.println(str);
			} else {
				System.out.print(str);
			}
		}
		
		if (level >= levelRequiredToReport && textArea != null) {
			textArea.setText(textArea.getText() + str+ (line?"\r\n":""));
			textArea.setCaretPosition(textArea.getDocument().getLength());
		}

		if (level >= levelRequiredToReport && logging) {
			try {
		        writer = new PrintWriter(new FileWriter(filename, true));
		        if (line) {
		        	writer.println(str);
		        } else {
		        	writer.print(str);
		        }
		        writer.close();
	        } catch (Exception e) {
		        System.err.println("Error writing to "+filename);
		        e.printStackTrace();
	        }
		}
	}

	public void reportError(String err) {
		reportError(err, true, true);
	}
	
	public void reportError(String err, boolean line, boolean reportToScreen) {
		reportError(err, line, reportToScreen, 0);
	}
	
	public void reportError(String err, boolean line, boolean reportToScreen, int levelRequiredToReport) {
		PrintWriter writer;
		
		if (level >= levelRequiredToReport && reportToScreen) {
			if (line) {
				System.err.println(err);
			} else {
				System.err.print(err);
			}
		}

		if (level >= levelRequiredToReport && textArea != null) {
			textArea.setText(textArea.getText() + err+ (line?"\r\n":""));
			textArea.setCaretPosition(textArea.getDocument().getLength());
		}

		if (level >= levelRequiredToReport && logging) {
			try {
		        writer = new PrintWriter(new FileWriter(filename, true));
		        if (line) {
		        	writer.println(err);
		        } else {
		        	writer.print(err);
		        }
		        writer.close();
	        } catch (Exception e) {
		        System.err.println("Error writing to "+filename);
		        e.printStackTrace();
	        }
		}
	}

	public void reportException(Exception e) {
		reportException(e, 0);
	}
	
	public void reportException(Exception e, int levelRequiredToReport) {
		PrintWriter writer;
		
		e.printStackTrace();
		if (level >= levelRequiredToReport && logging) {
			e.printStackTrace();
			try {
		        writer = new PrintWriter(new FileWriter(filename, true));
				e.printStackTrace(writer);
		        writer.close();
	        } catch (Exception e2) {
		        System.err.println("Error writing to "+filename);
		        e2.printStackTrace();
	        }
		}
		
		if (level >= 11) {
			System.exit(1);
		}
	}
	
	public void timestamp() {
		report(ext.getDate()+"\t"+ext.getTime());
	}
	
	public long memoryTotal(){
		long memory;
		
		report("Total heap size is: " +ext.prettyUpSize(memory = Runtime.getRuntime().totalMemory(), 1));
		
		return memory;
	}

	public long memoryFree(){
		long memory;
		
		report("Free heap size is: " +ext.prettyUpSize(memory = Runtime.getRuntime().freeMemory(), 1));
		
		return memory;
	}

	public double memoryPercentFree(){
		double percentFree;
		
		percentFree = ((float) 100*Runtime.getRuntime().freeMemory() / Runtime.getRuntime().totalMemory());
		
		report("Percent free heap size is: " +ext.formDeci(percentFree,1) + "%");
		
		return percentFree;
	}

	public long memoryUsed(){
		long memory;
		
		report("Used heap size is: " +ext.prettyUpSize(memory = (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()), 1));
		
		return memory;
	}

	public long memoryMax(){
		long memory;
		
		report("Used heap size is: " +ext.prettyUpSize(memory = Runtime.getRuntime().maxMemory(), 1));
		
		return memory;
	}	
}
