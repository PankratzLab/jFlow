/**
 * 
 */
package org.genvisis.seq.analysis.mtdna;

/**
 * Adapted from http://www.geeksforgeeks.org/longest-prefix-matching-a-trie-based-solution-in-java/
 * 
 * Tries to return the best haplogroup match, given that there is a unique sample available
 *
 */
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

// HaploNode, which stores a character of a Haplogroup and the children in a HashMap
class HaploNode {
	private final char value;
	private HashMap<Character, HaploNode> children;
	private HashSet<String> inds;// holds available samples
	private boolean bIsEnd;

	public HaploNode(char ch, HashSet<String> inds) {
		value = ch;
		children = new HashMap<Character, HaploNode>();
		this.inds = inds;
		bIsEnd = false;
	}

	public Map<Character, HaploNode> getChildren() {
		return children;
	}

	public void setIsEnd(boolean val) {
		bIsEnd = val;
	}

	public boolean isEnd() {
		return bIsEnd;
	}

	public HashSet<String> getInds() {
		return inds;
	}
}


/**
 * Implements the actual Trie
 *
 */
class HaploTrie {
	private HaploNode root;
	private HashSet<String> used;

	public HaploTrie() {
		root = new HaploNode((char) 0, new HashSet<String>());
		used = new HashSet<String>();
	}

	/**
	 * @param haplogroup the haplogroup to insert
	 * @param inds the sample ids with this haplogroup
	 */
	public void insert(String haplogroup, HashSet<String> inds) {

		// Find length of the given word
		int length = haplogroup.length();
		HaploNode crawl = root;
		crawl.getInds().addAll(inds);

		// Traverse through all characters of given word
		for (int level = 0; level < length; level++) {
			Map<Character, HaploNode> child = crawl.getChildren();
			char ch = haplogroup.charAt(level);
			// If there is already a child for current character of given word
			if (child.containsKey(ch)) {
				crawl = child.get(ch);

			} else // Else create a child
			{
				HaploNode temp = new HaploNode(ch, new HashSet<String>());
				child.put(ch, temp);
				crawl = temp;
			}
			crawl.getInds().addAll(inds);


		}

		// Set bIsEnd true for last character
		crawl.setIsEnd(true);
		crawl.getInds().addAll(inds);

	}

	/**
	 * Holds results returned from {@link HaploTrie#getBestMatchHaplotypeSamples(String)}
	 *
	 */
	public static class HaplogroupMatchResult {
		private String sample;
		private String closestMatch;

		private HaplogroupMatchResult(String sample, String closestMatch) {
			super();
			this.sample = sample;
			this.closestMatch = closestMatch;
		}

		public String getSample() {
			return sample;
		}

		public String getClosestMatch() {
			return closestMatch;
		}

	}

	/**
	 * @param haplotype get the best remaining haplogropu match
	 * @return
	 */
	public HaplogroupMatchResult getBestMatchHaplotypeSamples(String haplotype) {
		StringBuilder result = new StringBuilder(""); // Initialize resultant string
		int length = haplotype.length(); // Find length of the input string
		// Initialize reference to traverse through Trie
		HaploNode crawl = root;

		// Iterate through all characters of input string 'str' and traverse
		// down the Trie
		int level;
		int prevMatch = 0;
		for (level = 0; level < length; level++) {
			// Find current character of str
			char ch = haplotype.charAt(level);

			// HashMap of current Trie node to traverse down
			Map<Character, HaploNode> child = crawl.getChildren();
			if (child.containsKey(ch)) {
				crawl.getInds().removeAll(used);
				child.get(ch).getInds().removeAll(used);
			}
			if (child.containsKey(ch)&& !child.get(ch).getInds().isEmpty()
					&& !crawl.getInds().isEmpty()) {
				result.append(ch); // Update result

				crawl = child.get(ch); // Update crawl to move down in Trie


				// If this is end of a word or no samples remain, then update prevMatch
				if (crawl.isEnd())
					prevMatch = level + 1;
			} else
				break;

		}

		// If the last processed character did not match end of a word,
		// return the previously matching prefix
		String samp = crawl.getInds().iterator().next();
		used.add(samp);
		crawl.getInds().removeAll(used);


		if (!crawl.isEnd())
			return new HaplogroupMatchResult(samp, result.substring(0, prevMatch));
		else
			return new HaplogroupMatchResult(samp, result.toString());
	}
}


