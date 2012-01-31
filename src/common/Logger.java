package common;

import java.io.*;

public class Logger {
	public static final int LEVEL_ONE = 1;
	
	private String filename;
	private boolean logging;
	private int level;
	
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
	}
	
	public String getFilename() {
		return filename;
	}

	public int getLevel() {
		return level;
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
			try {
		        writer = new PrintWriter(new FileWriter(filename, true));
				e.printStackTrace(writer);
		        writer.close();
	        } catch (Exception e2) {
		        System.err.println("Error writing to "+filename);
		        e2.printStackTrace();
	        }
		}
	}
	
	public void timestamp() {
		report(ext.getDate()+"\t"+ext.getTime());
	}
}
