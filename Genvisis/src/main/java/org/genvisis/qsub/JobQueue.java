package org.genvisis.qsub;

import java.util.ArrayList;

public final class JobQueue {
	
	public enum QueueType {
		EXEC, ROUTE;
	}

	private String name;
	private boolean allowed;
	private int minWalltime = -1; // may not exist
	private int maxWalltime = -1; // hours
	private int defaultWalltime = -1; // hours
	private long defaultMem = -1; // bytes
	private long minMem = -1; // bytes
	private long maxMem = -1; // bytes
	private int jobsInQueue = -1; // VERY VOLATILE
	private int defaultProc = -1;
	private int minProc = -1; // may not exist
	private int maxProc = -1;
	private int defaultNodeCnt = -1;
	private int minNodeCnt = -1;
	private int maxNodeCnt = -1;
	private JobQueue.QueueType type;
	private ArrayList<String> routeDests = new ArrayList<String>();
	private ArrayList<String> allowedGroups = new ArrayList<String>();
	private boolean userAccessControlEnabled;
	private boolean groupAccessControlEnabled;
	private boolean groupControlSloppy = false;
	private ArrayList<String> allowedUsers = new ArrayList<String>();
	private boolean isDefaultQueue = false;

	public boolean isDefaultQueue() {
		return isDefaultQueue;
	}

	public void setDefaultQueue(boolean isDefaultQueue) {
		this.isDefaultQueue = isDefaultQueue;
	}

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

	public JobQueue.QueueType getType() {
		return type;
	}

	public void setType(JobQueue.QueueType type) {
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

	public int getDefaultWalltime() {
		return this.defaultWalltime;
	}

	public int getDefaultProc() {
		return this.defaultProc;
	}

	public int getDefaultNodeCnt() {
		return this.defaultNodeCnt;
	}

	public void setDefaultNodeCnt(int nodeCnt) {
		this.defaultNodeCnt = nodeCnt;
	}

	public void setDefaultProcCnt(int procCnt) {
		this.defaultProc = procCnt;
	}

	public void setDefaultWalltime(int defWall) {
		this.defaultWalltime = defWall;
	}

}