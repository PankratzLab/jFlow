package org.genvisis.cnv;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.QueueControl.JobQueue;

public class QueueProperties {
	
	private static final String PROPERTIES_FILE = "queues.properties";
	private static final String INIT_FAILURE = "Failed to create PBS queue properties: " + PROPERTIES_FILE
      + ". Genvisis will continue, but generated PBS files should be reviewed before submission."	;
	
	private static Logger log;
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
	
	
	private void parseLine(String line, JobQueue q) {
		String[] pts = line.split("=");
		QueueKeys.valueOf(pts[0]).setValue(q, pts[1]);
	}
	
	
	public static synchronized void voidQueues() {
		// clear queue container object
		
		if (Files.exists(PROPERTIES_FILE)) {
			boolean del = (new File(PROPERTIES_FILE)).delete();
			if (!del) {
				log.reportError("Could not delete PBS queue properties file.  Please delete this file manually and re-initialize.");
			}
		}
	}
	
	private static synchronized void init() {
//		if (props == null) {
//			props = new Properties();
//			try {
//				File propFile = new File(PROPERTIES_FILE);
//				if (!propFile.exists() && !propFile.createNewFile()) {
//					System.err.println(INIT_FAILURE);
//					System.exit(-1);
//				}
//				InputStream is = new FileInputStream(propFile);
//				props.load(is);
//				is.close();
//				
//				List<JobQueue> qList = QueueControl.parseAllowedQueues(log);
//				
//				
//				save();
//			} catch (Exception e) {
//				System.err.println(INIT_FAILURE);
//				System.exit(-1);
//			}
//		}
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
