package org.genvisis.one.george;

public class Node {
	private int id = -1;
	private int id2 = -1;
	private Node[] children = null;
	
	public Node(int id, Node[] children) {
		this.id = id;
		this.children = children;
	}
	public Node(int id, int id2, Node[] children) {
		this.id = id;
		this.id2 = id2;
		this.children = children;
	}
	
	public int[] getIds() {
		if (this.id2 == -1) // uninitialized
			return new int[]{this.id};
		else
			return new int[]{this.id, this.id2};
	}
	public Node[] getChildren() { // Returns pointer to children
		return this.children;
	}
}
