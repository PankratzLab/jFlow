package org.genvisis.common;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class QueueControl {

	public static class JobQueue {
		private String name;
		private boolean allowed;
		private int minWalltime = -1; // may not exist
		private int maxWalltime = -1; // hours
		private long defaultMem = -1; // bytes
		private long minMem = -1; // bytes
		private long maxMem = -1; // bytes
		private int jobsInQueue = -1; // VERY VOLATILE
		private int minProc = -1; // may not exist
		private int maxProc = -1;
		private int minNodeCnt = -1;
		private int maxNodeCnt = -1;
		private QueueType type;
		private ArrayList<String> routeDests = new ArrayList<String>();
		private ArrayList<String> allowedGroups = new ArrayList<String>();
		private boolean userAccessControlEnabled;
		private boolean groupAccessControlEnabled;
		private boolean groupControlSloppy = false;
		private ArrayList<String> allowedUsers = new ArrayList<String>();

		public String getName() {
			return name;
		}

		public void setName(String name) {
			this.name = name;
		}

		public boolean isAllowed() {
			return allowed;
		}

		public void setAllowed(boolean allowed) {
			this.allowed = allowed;
		}

		public int getMaxWalltime() {
			return maxWalltime;
		}

		public void setMaxWalltime(int maxWalltime) {
			this.maxWalltime = maxWalltime;
		}

		public long getDefaultMem() {
			return defaultMem;
		}

		public void setDefaultMem(long defaultMem) {
			this.defaultMem = defaultMem;
		}

		public long getMinMem() {
			return minMem;
		}
		
		public void setMinMem(long minMem) {
			this.minMem = minMem;
		}
		
		public long getMaxMem() {
			return maxMem;
		}
		
		public void setMaxMem(long maxMem) {
			this.maxMem = maxMem;
		}

		public int getJobsInQueue() {
			return jobsInQueue;
		}

		public void setJobsInQueue(int jobsInQueue) {
			this.jobsInQueue = jobsInQueue;
		}

		public int getMinProc() {
			return minProc;
		}

		public void setMinProc(int minProc) {
			this.minProc = minProc;
		}

		public QueueType getType() {
			return type;
		}

		public void setType(QueueType type) {
			this.type = type;
		}

		public ArrayList<String> getRouteDests() {
			return routeDests;
		}

		public void setRouteDests(ArrayList<String> routeDests) {
			this.routeDests = routeDests;
		}

		public int getMaxProc() {
			return maxProc;
		}

		public void setMaxProc(int maxProc) {
			this.maxProc = maxProc;
		}

		public int getMinWalltime() {
			return minWalltime;
		}

		public void setMinWalltime(int minWalltime) {
			this.minWalltime = minWalltime;
		}

		public int getMinNodeCnt() {
			return minNodeCnt;
		}

		public void setMinNodeCnt(int minNodeCnt) {
			this.minNodeCnt = minNodeCnt;
		}

		public int getMaxNodeCnt() {
			return maxNodeCnt;
		}

		public void setMaxNodeCnt(int maxNodeCnt) {
			this.maxNodeCnt = maxNodeCnt;
		}

		public ArrayList<String> getAllowedGroups() {
			return allowedGroups;
		}

		public ArrayList<String> getAllowedUsers() {
			return allowedUsers;
		}

		public void setAllowedUsers(ArrayList<String> allowedUsers) {
			this.allowedUsers = allowedUsers;
		}

		public boolean isGroupControlSloppy() {
			return groupControlSloppy;
		}

		public void setGroupControlSloppy(boolean groupControlSloppy) {
			this.groupControlSloppy = groupControlSloppy;
		}

		public void setUserAccessControlEnabled(boolean control) {
			this.userAccessControlEnabled = control;
		}

		public void setGroupAccessControlEnabled(boolean control) {
			this.groupAccessControlEnabled = control;
		}

		public boolean isGroupAccessControlEnabled() {
			return this.groupAccessControlEnabled;
		}

		public boolean isUserAccessControlEnabled() {
			return this.userAccessControlEnabled;
		}

	}

	enum QueueType {
		EXEC,
		ROUTE;
	}

	private static final String TAG_WALLTIME_MIN = "resources_min.walltime = ";
	private static final String TAG_WALLTIME_MAX = "resources_max.walltime = ";
	private static final String TAG_PROC_MIN = "resources_min.procct = ";
	private static final String TAG_PROC_MAX = "resources_max.procct = ";
	private static final String TAG_NODE_MIN = "resources_min.nodect = ";
	private static final String TAG_NODE_MAX = "resources_max.nodect = ";
	private static final String TAG_TYPE = "queue_type = ";
	private static final String TAG_TOTAL_JOBS = "total_jobs = ";
	private static final String TAG_USERS_LIST = "acl_users = ";
	private static final String TAG_USER_CTRL = "acl_user_enable = ";
	private static final String TAG_GROUP_CTRL = "acl_group_enable = ";
	private static final String TAG_GROUPS = "acl_groups = ";
	private static final String TAG_GROUP_CTRL_SLOPPY = "acl_group_sloppy = ";
	private static final String TAG_ROUTE_DEST = "route_destinations = ";
	private static final String TAG_MEM = "resources_default.mem = ";

	private static String[] loadIDInfo() throws IOException {
		Runtime rt = Runtime.getRuntime();
		String[] commands = {"id"};
		Process proc = rt.exec(commands);

		BufferedReader stdInput = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		String s = stdInput.readLine();
		stdInput.close();
		return s.split("\\s");
	}

	/**
	 * Runs 'id' command and parses current user name. 
	 * 
	 * @return
	 */
	public static String getUserName() {
		String uid = "";
		try {
			uid = loadIDInfo()[0];
			uid = uid.substring(uid.indexOf('(') + 1, uid.length() - 1);
		} catch (IOException e) {
			// error from running 'id' or from reading input; either way, unable to determine username
		}
		return uid;
	}

	/**
	 * Runs 'id' command and parses user's current group. 
	 * 
	 * @return
	 * @throws IOException Thrown by Runtime.exec and BufferedReader when reading command output
	 */
	public static String getCurrentGroup() {
		String gid = "";
		try {
			gid = loadIDInfo()[1];
			gid = gid.substring(gid.indexOf('(') + 1, gid.length() - 1);
		} catch (IOException e) {
			// error from running 'id' or from reading input; either way, unable to determine group
		}
		return gid;
	}

	/**
	 * Runs 'id' command and parses all groups for current user. 
	 * 
	 * @return
	 * @throws IOException Thrown by Runtime.exec and BufferedReader when reading command output
	 */
	public static ArrayList<String> getUserGroups() {
		String[] grps = new String[0];
		try {
			loadIDInfo()[2].split(",");
		} catch (IOException e) {
			// error from running 'id' or from reading input; either way, unable to determine groups
		}
		ArrayList<String> groups = new ArrayList<String>();
		for (String g : grps) { 
			groups.add(g.substring(g.indexOf('(') + 1, g.length() - 1));
		}
		return groups;
	}
	
	
	public static ArrayList<JobQueue> parseAllowedQueues(Logger log) {
		return parseAllowedQueues(getUserName(), getCurrentGroup(), getUserGroups().toArray(new String[0]), log);
	}

	public static ArrayList<JobQueue> parseAllowedQueues(String username, String currgroup, String[] allGroups, Logger log) {
		return filterAllowed(parseQueues(username, currgroup, allGroups, log));
	}

	/**
	 * Runs 'qstat -Qf' command and parses queue information. 
	 * 
	 * @return
	 */
	public static ArrayList<JobQueue> parseQueues(String username, String currGroup, String[] allGroups, Logger log) {
		Runtime rt = Runtime.getRuntime();
		String[] commands = {"qstat", "-Qf"};
		BufferedReader stdInput = null;
		
		try {
			Process proc = rt.exec(commands);
			stdInput = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		} catch (IOException e) {
			log.reportError("Failed to parse available queues: " + e.getMessage());
		}

		ArrayList<JobQueue> queues = new ArrayList<JobQueue>();
		
		if (stdInput == null) {
			return queues;
		}
		
		JobQueue curr = null;
		String s;
		try {
			while ((s = stdInput.readLine()) != null) {
				s = s.trim();
				if (s.startsWith("Queue: ")) {
					if (curr != null) {
						checkIfQueueAllowed(curr, username, currGroup, allGroups);
						queues.add(curr);
					}
					curr = new JobQueue();
					curr.setName(s.substring("Queue: ".length()).trim());
				}
				if (s.startsWith(TAG_TYPE)) {
					if (s.endsWith("Route")) {
						curr.setType(QueueType.ROUTE);
					} else {
						curr.setType(QueueType.EXEC);
					}
				}
				if (s.startsWith(TAG_TOTAL_JOBS)) {
					try {
						curr.setJobsInQueue(Integer.parseInt(s.substring(TAG_TOTAL_JOBS.length()).trim()));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'total jobs' tag but couldn't parse value: " + s.substring(13).trim());
					}
				}
				if (s.startsWith(TAG_WALLTIME_MIN)) {
					String s1 = s.substring(TAG_WALLTIME_MIN.length()).trim();
					try {
						curr.setMinWalltime(Integer.parseInt(s1.substring(0, s1.indexOf(':'))));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'minimum walltime' tag but couldn't parse value: " + s1.substring(0, s1.indexOf(':')).trim());
					}
				}
				if (s.startsWith(TAG_WALLTIME_MAX)) {
					String s1 = s.substring(TAG_WALLTIME_MAX.length()).trim();
					try {
						curr.setMaxWalltime(Integer.parseInt(s1.substring(0, s1.indexOf(':'))));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'maximum walltime' tag but couldn't parse value: " + s1.substring(0, s1.indexOf(':')).trim());
					}
				}
				if (s.startsWith(TAG_MEM)) {
					String b = s.substring(TAG_MEM.length());
					try {
						curr.setDefaultMem(Long.parseLong(b.substring(0, b.length() - 1)));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'assigned memory' tag but couldn't parse value: " + b.substring(0, b.indexOf(':')).trim());
					}
				}
				if (s.startsWith(TAG_PROC_MAX)) {
					try {
						curr.setMaxProc(Integer.parseInt(s.substring(TAG_PROC_MAX.length())));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'maximum processors' tag but couldn't parse value: " + s.substring(TAG_PROC_MAX.length()));
					}
				}
				if (s.startsWith(TAG_PROC_MIN)) {
					try {
						curr.setMinProc(Integer.parseInt(s.substring(TAG_PROC_MIN.length())));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'minimum processors' tag but couldn't parse value: " + s.substring(TAG_PROC_MIN.length()));
					}
				}
				if (s.startsWith(TAG_NODE_MAX)) {
					try {
						curr.setMaxNodeCnt(Integer.parseInt(s.substring(TAG_NODE_MAX.length())));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'maximum nodes' tag but couldn't parse value: " + s.substring(TAG_NODE_MAX.length()));
					}
				}
				if (s.startsWith(TAG_NODE_MIN)) {
					try {
						curr.setMinNodeCnt(Integer.parseInt(s.substring(TAG_NODE_MIN.length())));
					} catch (NumberFormatException e1) {
						log.reportError("Found 'minimum processors' tag but couldn't parse value: " + s.substring(TAG_NODE_MIN.length()));
					}
				}
				if (s.startsWith(TAG_USERS_LIST)) {
					String[] allowedUsers = s.substring(TAG_USERS_LIST.length()).split(",");
					for (String u : allowedUsers) {
						curr.getAllowedUsers().add(u);
					}
				}
				if (s.startsWith(TAG_USER_CTRL)) {
					curr.setUserAccessControlEnabled(s.toUpperCase().endsWith("TRUE"));
				}
				if (s.startsWith(TAG_GROUP_CTRL)) {
					curr.setGroupAccessControlEnabled(s.toUpperCase().endsWith("TRUE"));
				}
				if (s.startsWith(TAG_GROUP_CTRL_SLOPPY)) {
					curr.setGroupControlSloppy(s.toUpperCase().endsWith("TRUE"));
				}
				if (s.startsWith(TAG_ROUTE_DEST)) {
					String[] dests = s.substring(TAG_ROUTE_DEST.length()).trim().split(",");
					for (String d : dests) {
						curr.getRouteDests().add(d);
					}
				}
				if (s.startsWith(TAG_GROUPS)) {
					String[] groups = s.substring(TAG_ROUTE_DEST.length()).trim().split(",");
					for (String g : groups) {
						curr.getAllowedGroups().add(g);
					}
				}
			}
		} catch (IOException e) {
			log.reportError("Error reading data from qstat command: " + e.getMessage());
		}
		if (curr != null) {
			checkIfQueueAllowed(curr, username, currGroup, allGroups);
			queues.add(curr);
		}
		try {
			stdInput.close();
		} catch (IOException e) {
			// ignore - probably a bad idea
		}

		return queues;
	}


	private static void checkIfQueueAllowed(JobQueue curr, String username, String currGroup, String[] allGroups) {
		if (!curr.isGroupAccessControlEnabled() && !curr.isUserAccessControlEnabled()) {
			curr.setAllowed(true);
		} else {
			if (curr.isUserAccessControlEnabled()) {
				if (username != null && curr.getAllowedUsers().contains(username)) {
					curr.setAllowed(true);
				}
			}
			if (curr.isGroupAccessControlEnabled()) {
				if (curr.isGroupControlSloppy()) {
					if (allGroups != null && allGroups.length > 0) {
						for (String s : allGroups) {
							if (curr.getAllowedGroups().contains(s)) {
								curr.setAllowed(true);
								break;
							}
						}
					} else if (currGroup != null) {
						if (curr.getAllowedGroups().contains(currGroup)) {
							curr.setAllowed(true);
						}
					}
				} else {
					if (currGroup != null && curr.getAllowedGroups().contains(currGroup)) {
						curr.setAllowed(true);
					}
				}
			}
		}
	}

	private static ArrayList<JobQueue> filterAllowed(ArrayList<JobQueue> all) {
		ArrayList<JobQueue> filt = new ArrayList<JobQueue>();
		for (JobQueue f : all) {
			if (f.isAllowed()) {
				filt.add(f);
			}
		}
		return filt;
	}

	/**
	 * <p>
	 * Takes a list of PBS queues, filters to usable queues, and returns a queue based on (in
	 * descending order) three metrics:<br />
	 * <ol>
	 * <li>An access-controlled queue (any; if multiple, the queue with the fewest currently-running
	 * jobs).</li>
	 * <li>A routing queue with the largest number of sub-queues allowed for this user.</li>
	 * <li>The queue with the most number of jobs in it (assuming that this queue is intended by
	 * system administrators to be the default queue).</li>
	 * </ol>
	 * 
	 * @param queues List of PBS Queues
	 * @return
	 */
	public static JobQueue findSensibleDefault(ArrayList<JobQueue> queues) {
		ArrayList<JobQueue> allowed = filterAllowed(queues); // maybe already filtered, but just in case

		JobQueue routeQueueMost = null;
		int svRt = 0;
		JobQueue exclusive = null;
		JobQueue mostUsed = null;
		for (JobQueue q : allowed) {
			if (mostUsed == null || mostUsed.getJobsInQueue() < q.getJobsInQueue()) {
				mostUsed = q;
			}

			if (q.getType() == QueueType.ROUTE) {
				int found = 0;
				for (String r : q.getRouteDests()) {
					for (JobQueue q2 : allowed) {
						if (q2.getName().equals(r)) {
							found++;
							break;
						}
					}
				}
				if (routeQueueMost == null || svRt < found) {
					routeQueueMost = q;
					svRt = found;
				}
			}

			if ((q.isUserAccessControlEnabled() || q.isGroupAccessControlEnabled())
					&& (exclusive == null || exclusive.getJobsInQueue() > q.getJobsInQueue())) {
				exclusive = q;
			}
		}

		return exclusive != null ? exclusive : routeQueueMost != null ? routeQueueMost : mostUsed;
	}
	
}
