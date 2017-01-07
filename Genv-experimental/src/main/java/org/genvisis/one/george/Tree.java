package org.genvisis.one.george;

public class Tree {
	Node root;
	int maxLeafId;
	
	public Integer[] shallowcopy(Integer[] ls) {
		Integer[] ls2 = new Integer[ls.length];
		for (int i = 0; i < ls.length; i++)
			ls2[i] = ls[i];
		return ls2;
	}
	
	public Tree(Node root, int maxLeafId) {
		this.root = root;
		this.maxLeafId = maxLeafId;
	}
	// Returns an array s.t. 1 indicates leafId is a member of that parent (e.g. T cell is a lymphocyte)
	// leafId is one indexed (as is the input)
	// leafId = 0 means not in any category and not in tree. Returns array of all 0's
	private Tuple<Integer[], Boolean> getArrayHelp(int leafId, Node n, Integer[] storeArr) 
			throws Exception {
		if (leafId == 0)
			return new Tuple<Integer[], Boolean>(storeArr, true);
		
		Integer[] newStoreArr = shallowcopy(storeArr);
		
		for (int i = 0; i < n.getIds().length; i++) {
			newStoreArr[n.getIds()[i]-1] = 1; // zero indexed
		}
		if (newStoreArr[leafId-1] == 1) {
			System.out.println("true");
			return new Tuple<Integer[], Boolean>(newStoreArr, true);
		}
		if (n.getChildren() != null) {
			for (int i = 0; i < n.getChildren().length; i++) {
				Node child = n.getChildren()[i];
				Tuple<Integer[], Boolean> tup = getArrayHelp(leafId, child, newStoreArr);
				//Integer[] arr = tup.x;
				boolean success = tup.y;
				if (success)
					return tup;
			}
		}
		return new Tuple<Integer[], Boolean>(null, false);
	}
	public Integer[] getArray(int leafId) throws Exception {
		Integer[] storeArr = new Integer[maxLeafId]; // leaf ids from 1 to maxLeafId
		for (int i = 0; i < storeArr.length; i++)
			storeArr[i] = 0;
		Tuple<Integer[], Boolean> tup = getArrayHelp(leafId, root, storeArr);
		return tup.x;
	}
}
