package org.genvisis.qsub;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;

public final class QueueProperties {
	
	private static final String PROPERTIES_FILE = "queues.properties";
	
	private static Logger log = new Logger();
	private static List<JobQueue> qList;
	
	private enum QueueKeys {
		Q_NAME() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getName(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setName(rawValue); }
		},
		Q_WALL_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinWalltime(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMinWalltime(parseIntOrNeg1(rawValue)); }
		},
		Q_WALL_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxWalltime(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMaxWalltime(parseIntOrNeg1(rawValue)); }
		},
		Q_WALL_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultWalltime(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setDefaultWalltime(parseIntOrNeg1(rawValue)); }
		},
		Q_MEM_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinMem(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMinMem(parseLongOrNeg1(rawValue)); }
		},
		Q_MEM_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxMem(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMaxMem(parseLongOrNeg1(rawValue)); }
		},
		Q_MEM_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultMem(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setDefaultMem(parseLongOrNeg1(rawValue)); }
		},
		Q_PROC_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinProc(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMinProc(parseIntOrNeg1(rawValue)); }
		},
		Q_PROC_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxProc(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMaxProc(parseIntOrNeg1(rawValue)); }
		},
		Q_PROC_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultProc(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setDefaultProcCnt(parseIntOrNeg1(rawValue)); }
		},
		Q_NODE_MIN() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMinNodeCnt(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMinNodeCnt(parseIntOrNeg1(rawValue)); }
		},
		Q_NODE_MAX() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getMaxNodeCnt(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setMaxNodeCnt(parseIntOrNeg1(rawValue)); }
		},
		Q_NODE_DEF() {
			@Override
			public String getValue(JobQueue q) { return this.toString() + "=" + q.getDefaultNodeCnt(); }
			@Override
			public void setValue(JobQueue q, String rawValue) { q.setDefaultNodeCnt(parseIntOrNeg1(rawValue)); }
		};
		
		public abstract String getValue(JobQueue q);
		public abstract void setValue(JobQueue q, String rawValue);
		
		protected long parseLongOrNeg1(String v) {
			long val = -1;
			try {
				val = Long.parseLong(v);
			} catch (NumberFormatException e) {}
			return val;
		}
		private static int parseIntOrNeg1(String v) {
			int val = -1;
			try {
				val = Integer.parseInt(v);
			} catch (NumberFormatException e) {}
			return val;
		}
	}
	
	private QueueProperties() {}
	
	public static List<JobQueue> getJobQueues() {
		if (qList == null) {
			load();
		}
		return qList == null ? new ArrayList<JobQueue>() : qList;
	}
	
	public static synchronized void load() {
		if (Files.exists(PROPERTIES_FILE)) {
			loadFile();
		} else {
			log.report("No queue properties file found; Initializing programmatically.  Please edit the generated queue.properties file by hand with relevant information.");
			init();
		}
	}
	
	private static void loadFile() {
		BufferedReader reader;
		try {
			reader = Files.getAppropriateReader(PROPERTIES_FILE);
		} catch (FileNotFoundException e) {
			log.reportError("Couldn't read properties file: " + PROPERTIES_FILE);
			log.reportException(e);
			return;
		}
		String line;
		JobQueue q = null;
		qList = new ArrayList<JobQueue>();
		try {
			while ((line = reader.readLine()) != null) {
				if (line.startsWith(QueueKeys.Q_NAME.toString())) {
					if (q != null) {
						qList.add(q);
					}
					q = new JobQueue();
				}
				parseLine(line, q);
			}
		} catch (IOException e) {
			log.reportError("Encountered error when reading queue properties file.");
			log.reportException(e);
			qList = null;
			return;
		}
		qList.add(q);
	}
	
	private static void parseLine(String line, JobQueue q) {
		String[] pts = line.split("=");
		QueueKeys.valueOf(pts[0]).setValue(q, pts[1]);
	}
	
	
	public static synchronized void voidQueues() {
		qList = null;
		if (Files.exists(PROPERTIES_FILE)) {
			boolean del = (new File(PROPERTIES_FILE)).delete();
			if (!del) {
				log.reportError("Could not delete PBS queue properties file.  Please delete this file manually and re-initialize.");
			}
		}
	}
	
	private static synchronized void init() {
		if (qList == null) {
			List<JobQueue> allQs = QueuesParser.parseAllowedQueues(log);
			qList = new ArrayList<JobQueue>();
			for (JobQueue q : allQs) {
				if (q.isAllowed()) {
					qList.add(q);
				}
			}
			save();
		}
	}
	
	
	/**
	 * Write the current properties to disk
	 */
	private static synchronized void save() {
		PrintWriter out;
		
		if (qList != null) {
			try {
				out = Files.getAppropriateWriter(PROPERTIES_FILE);
				QueueKeys[] keys = QueueKeys.values();
				
				for (JobQueue q : qList) {
					out.println(QueueKeys.Q_NAME.getValue(q));
					for (int i = 1; i < keys.length; i++) {
						out.println(keys[i].getValue(q));
					}
					out.println();
				}
			  out.flush();
				out.close();
			} catch (Exception e) {
				System.err.println("Failed to save PBS queue properties: " + PROPERTIES_FILE);
			}
		}
	}

	
}
